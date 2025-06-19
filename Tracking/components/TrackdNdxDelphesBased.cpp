#include "TrackdNdxDelphesBased.h"

// Gaudi
#include <GaudiKernel/MsgStream.h>

// EDM4hep
#include "edm4hep/Quantity.h"
#include "edm4hep/utils/vector_utils.h"

// dd4hep
#include "DD4hep/Detector.h"

// ROOT
#include <TVectorD.h>

// STL
#include <random>
#include <limits>

DECLARE_COMPONENT(TrackdNdxDelphesBased)

TrackdNdxDelphesBased::TrackdNdxDelphesBased(const std::string& name, ISvcLocator* svcLoc)
    : Transformer(name, svcLoc,
        {KeyValues("InputLinkCollection", {"TrackMCParticleLinks"}),
         KeyValues("HeaderName", {"EventHeader"})},
        {KeyValues("OutputCollection", {"RecDqdxCollection"})}) {}

StatusCode TrackdNdxDelphesBased::initialize() {
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

    // Initialise the delphes track util
    m_delphesTrkUtil = new TrkUtil();
    if (!m_delphesTrkUtil) {
        error() << "Failed to create delphes TrkUtil instance" << endmsg;
        return StatusCode::FAILURE;
    }

    // Load the geometry parameters from the XML file
    // Cast to mm as this algorithm tries to use mm or 1/mm everywhere
    double Rmin = m_geoSvc->getDetector()->constantAsDouble(m_Rmin_parameter_name.value())/dd4hep::mm;
    double Rmax = m_geoSvc->getDetector()->constantAsDouble(m_Rmax_parameter_name.value())/dd4hep::mm;
    double Zmax = m_geoSvc->getDetector()->constantAsDouble(m_Zmax_parameter_name.value())/dd4hep::mm;
    double Zmin;
    // For fwd-bwd symmetric detectors, use negative Zmax as Zmin
    if (m_Zmin_parameter_name.value() == m_Zmax_parameter_name.value()) {
        Zmin = -Zmax;
    } else {
        Zmin = m_geoSvc->getDetector()->constantAsDouble(m_Zmin_parameter_name.value())/dd4hep::mm;
    }

    debug() << "Geometry parameters:" << endmsg;
    debug() << "Rmin: " << Rmin << " mm" << endmsg;
    debug() << "Rmax: " << Rmax << " mm" << endmsg;
    debug() << "Zmin: " << Zmin << " mm" << endmsg;
    debug() << "Zmax: " << Zmax << " mm" << endmsg;

    m_delphesTrkUtil->SetDchBoundaries(Rmin, Rmax, Zmin, Zmax);
    m_delphesTrkUtil->SetGasMix(m_GasSel.value());

    // Make sure fill factor is between 0 and 1
    if (m_fill_factor.value() < 0.0 || m_fill_factor.value() > 1.0) {
        warning() << "Fill factor of "<< m_fill_factor.value() << " is not between 0 and 1, setting to 1.0" << endmsg;
        m_fill_factor.set(1.0);
    }

    return StatusCode::SUCCESS;
}

edm4hep::RecDqdxCollection TrackdNdxDelphesBased::operator()(const edm4hep::TrackMCParticleLinkCollection& input,
                                                            const edm4hep::EventHeaderCollection& header) const {
    edm4hep::RecDqdxCollection outputCollection;

    std::mt19937_64 random_engine;
    auto engine_seed = m_uniqueIDSvc->getUniqueID(header, this->name());
    random_engine.seed(engine_seed);

    // Dummy value if dN/dx calculation somehow fails
    const double dummy_value = -1.0;
    bool success = true;

    debug() << "Processing new Event" << endmsg;

    unsigned int i = 0;
    for (const auto& link : input) {
        debug() << "Processing track " << i++ << endmsg;

        // Get the track and corresponding MCParticle
        const auto& mc_particle = link.getTo();
        const auto& track = link.getFrom();
        // Check for validity
        if (!mc_particle.isAvailable() || !track.isAvailable()) {
            warning() << "Invalid link found, skipping." << endmsg;
            continue;
        }

        //////////////////////////
        // Particle Information //
        //////////////////////////
        double momentum = edm4hep::utils::magnitude(mc_particle.getMomentum());
        double mass     = mc_particle.getMass();

        double betagamma = momentum/mass;
        debug() << "MCParticle betagamma: " << betagamma << endmsg;
        // Check if betagamma is in valid range of delphes parametrisation (status: 16 June 2025)
        if (betagamma < 0.5 || betagamma > 20000.0) {
            warning() << "beta*gamma value outside of \"good\" range of delphes parametrisation (0.5-20000), dN/dx will be set to dummy value: " 
                      << dummy_value << " clusters/mm" << endmsg;
            success = false;
        }

        // Get number of clusters per length from delphes
        // Output from delphes function is in 1/m, so to convert to 1/mm we need to scale accordingly
        double nclusters_per_mm = m_delphesTrkUtil->Nclusters(betagamma, m_GasSel.value()) / 1000.0;
        debug() << "Number of clusters per mm: " << nclusters_per_mm  << endmsg;

        ///////////////////////
        // Track Information //
        ///////////////////////
        // Use track state at IP, since this corresponds to delphes and energy loss in tracking system is negligible
        const auto& track_state = track.getTrackStates(edm4hep::TrackState::AtIP);

        // Convert edm4hep::TrackState to delphes parameters
        // Inverse conversion from https://github.com/key4hep/k4SimDelphes/blob/main/converter/src/DelphesEDM4HepConverter.cc#L532
        // Note the same order of parameters as in delphes: D0, phi, C, Z0, cot(theta)
        // See: https://github.com/delphes/delphes/blob/98f15add056e657e39bda9e32ccd97ef427ce04c/external/TrackCovariance/TrkUtil.cc#L961
        TVectorD delphes_track(5);
        delphes_track[0] = track_state.D0;
        delphes_track[1] = track_state.phi;
        const double scale = -2.0;            // delphes uses C instead of omega, scale is used to convert
        delphes_track[2] = track_state.omega / scale;
        delphes_track[3] = track_state.Z0;
        delphes_track[4] = track_state.tanLambda; // tanLambda and cot(theta) are the same thing (see DelphesEDM4HepConverter.cc)

        // Note: track length will already be in mm, since this is what delphes and the track parametrisation use
        // so no need to cast it to dd4hep::mm
        double track_length = m_delphesTrkUtil->TrkLen(delphes_track);
        // Check if track length calculation was successful
        if (track_length < std::numeric_limits<double>::epsilon()) {
            warning() << "Delphes track length calculation returned 0.0, dN/dx will be set to dummy value: " 
                      << dummy_value << " clusters/mm" << endmsg;
            success = false;
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

        double dNdx_value = (success)? (n_cluster/track_length) : dummy_value;
        debug() << "dNdx value: " << dNdx_value << " (clusters/mm)" << endmsg;

        auto dqdx = outputCollection.create();
        edm4hep::Quantity q;
        q.value = static_cast<float>(dNdx_value);

        dqdx.setDQdx(q);
        dqdx.setTrack(track);
    }


    return outputCollection;
}
