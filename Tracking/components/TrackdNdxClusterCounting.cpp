#include "TrackdNdxClusterCounting.h"

// Gaudi
#include <GaudiKernel/MsgStream.h>

// EDM4hep
#include "edm4hep/Quantity.h"
#include "edm4hep/utils/vector_utils.h"

// dd4hep
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"

// STL
#include <limits>
#include <random>

DECLARE_COMPONENT(TrackdNdxClusterCounting)
// Old name of this algorithm, from when it was linking against Delphes
DECLARE_COMPONENT_WITH_ID(TrackdNdxClusterCounting, "TrackdNdxDelphesBased")

TrackdNdxClusterCounting::TrackdNdxClusterCounting(const std::string& name, ISvcLocator* svcLoc)
    : Transformer(
          name, svcLoc,
          {KeyValues("InputLinkCollection", {"TrackMCParticleLinks"}), KeyValues("HeaderName", {"EventHeader"})},
          {KeyValues("OutputCollection", {"RecDqdxCollection"})}) {}

StatusCode TrackdNdxClusterCounting::initialize() {
  m_uniqueIDSvc = service("UniqueIDGenSvc");
  if (!m_uniqueIDSvc) {
    error() << "Unable to locate the UniqueIDGenSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate the GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  // Load the geometry parameters from the XML file
  // Cast to mm as this algorithm tries to use mm or 1/mm everywhere
  double Rmin = m_geoSvc->getDetector()->constantAsDouble(m_Rmin_parameter_name.value()) / dd4hep::mm;
  double Rmax = m_geoSvc->getDetector()->constantAsDouble(m_Rmax_parameter_name.value()) / dd4hep::mm;
  double Zmax = m_geoSvc->getDetector()->constantAsDouble(m_Zmax_parameter_name.value()) / dd4hep::mm;
  double Zmin;
  // For fwd-bwd symmetric detectors, use negative Zmax as Zmin
  if (m_Zmin_parameter_name.value() == m_Zmax_parameter_name.value()) {
    Zmin = -Zmax;
  } else {
    Zmin = m_geoSvc->getDetector()->constantAsDouble(m_Zmin_parameter_name.value()) / dd4hep::mm;
  }

  debug() << "Geometry parameters:" << endmsg;
  debug() << "Rmin: " << Rmin << " mm" << endmsg;
  debug() << "Rmax: " << Rmax << " mm" << endmsg;
  debug() << "Zmin: " << Zmin << " mm" << endmsg;
  debug() << "Zmax: " << Zmax << " mm" << endmsg;

  m_driftVolume = {Rmin, Rmax, Zmin, Zmax};

  const auto gas = k4RecTracker::ClusterCounting::toGasMixture(m_GasSel.value());
  if (!gas) {
    error() << "Unknown gas selection " << m_GasSel.value() << endmsg;
    return StatusCode::FAILURE;
  }
  m_clusterParametrisation = k4RecTracker::ClusterCounting::Parametrisation::forGas(*gas);

  // Make sure fill factor is between 0 and 1
  if (m_fill_factor.value() < 0.0 || m_fill_factor.value() > 1.0) {
    warning() << "Fill factor of " << m_fill_factor.value() << " is not between 0 and 1, setting to 1.0" << endmsg;
    m_fill_factor.set(1.0);
  }

  return StatusCode::SUCCESS;
}

edm4hep::RecDqdxCollection TrackdNdxClusterCounting::operator()(const edm4hep::TrackMCParticleLinkCollection& input,
                                                                const edm4hep::EventHeaderCollection& header) const {
  edm4hep::RecDqdxCollection outputCollection;

  std::mt19937_64 random_engine;
  auto engine_seed = m_uniqueIDSvc->getUniqueID(header, this->name());
  random_engine.seed(engine_seed);

  // Dummy value if dN/dx calculation somehow fails
  const double dummy_value = -999.0;

  debug() << "Processing new Event" << endmsg;

  unsigned int track_count = 0;
  for (const auto& link : input) {
    debug() << "Processing track " << track_count++ << endmsg;

    // Initialise dN/dx value to dummy value
    double dNdx_value = dummy_value;

    // Get the track and corresponding MCParticle
    const auto& mc_particle = link.getTo();
    const auto& track = link.getFrom();
    // Check for validity
    if (!mc_particle.isAvailable() || !track.isAvailable()) {
      warning() << "Invalid link found, skipping." << endmsg;
      continue;
    }

    const auto store_value = [&]() {
      debug() << "dNdx value: " << dNdx_value << " (clusters/mm)" << endmsg;

      auto dqdx = outputCollection.create();
      edm4hep::Quantity q;
      q.value = static_cast<float>(dNdx_value);

      dqdx.setDQdx(q);
      dqdx.setTrack(track);
    };

    //////////////////////////
    // Particle Information //
    //////////////////////////
    double momentum = edm4hep::utils::magnitude(mc_particle.getMomentum());
    double mass = mc_particle.getMass();

    double betagamma = momentum / mass;
    debug() << "MCParticle betagamma: " << betagamma << endmsg;
    // Check if betagamma is in valid range of the parametrisation
    if (betagamma < m_clusterParametrisation->minBetaGamma()) {
      debug() << "beta*gamma value below lower limit of range of the parametrisation ("
              << m_clusterParametrisation->minBetaGamma() << "-" << m_clusterParametrisation->maxBetaGamma()
              << "), dN/dx will be set to dummy value: " << dummy_value << " clusters/mm" << endmsg;
      store_value();
      continue;
    } else if (betagamma > m_clusterParametrisation->maxBetaGamma()) {
      debug() << "beta*gamma value above upper limit of range of the parametrisation ("
              << m_clusterParametrisation->minBetaGamma() << "-" << m_clusterParametrisation->maxBetaGamma()
              << "), the value at the upper limit (Fermi plateau) will be used as approximation." << endmsg;
    }

    // Get number of clusters per length (clamped to the range of the parametrisation)
    double nclusters_per_mm = m_clusterParametrisation->clustersPerMM(betagamma);
    debug() << "Number of clusters per mm: " << nclusters_per_mm << endmsg;

    ///////////////////////
    // Track Information //
    ///////////////////////
    // Use track state at IP, since energy loss in tracking system is negligible
    const auto track_state = track.getTrackState(edm4hep::TrackState::AtIP).value();

    // Note: track length will already be in mm, since this is what the track parametrisation uses, so no need to
    // cast it to dd4hep::mm
    double track_length = k4RecTracker::ClusterCounting::trackLengthInCylinder(
        track_state.D0, track_state.omega, track_state.Z0, track_state.tanLambda, m_driftVolume);
    // Check if track length calculation was successful
    if (track_length < std::numeric_limits<double>::epsilon()) {
      warning() << "Track does not cross the drift volume, dN/dx will be set to dummy value: " << dummy_value
                << " clusters/mm" << endmsg;
      store_value();
      continue;
    }
    debug() << "Track length inside full chamber: " << track_length << " mm" << endmsg;
    // Apply fill factor to track length
    track_length *= m_fill_factor.value();
    debug() << "Track length after applying fill factor: " << track_length << " mm" << endmsg;

    /////////////////////////////////////////////
    // Draw Number of Clusters from Poissonian //
    /////////////////////////////////////////////

    double nclusters_mean = nclusters_per_mm * track_length;
    std::poisson_distribution<int> poisson_dist(nclusters_mean);

    int n_cluster = poisson_dist(random_engine);
    debug() << "Track has " << n_cluster << " clusters." << endmsg;

    dNdx_value = n_cluster / track_length;

    store_value();
  }

  return outputCollection;
}
