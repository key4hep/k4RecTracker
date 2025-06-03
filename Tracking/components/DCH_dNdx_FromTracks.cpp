#include "DCH_dNdx_FromTracks.h"

// Gaudi
#include <GaudiKernel/MsgStream.h>

// EDM4hep
#include "edm4hep/Quantity.h"
#include "edm4hep/utils/vector_utils.h"

// dd4hep
#include "DD4hep/Detector.h"

// ROOT
#include <TVectorD.h>

DECLARE_COMPONENT(DCH_dNdx_FromTracks)

DCH_dNdx_FromTracks::DCH_dNdx_FromTracks(const std::string& name, ISvcLocator* svcLoc)
    : Transformer(name, svcLoc,
        {KeyValues("InputLinkCollection", {"TrackMCParticleLinks"}),
         KeyValues("HeaderName", {"EventHeader"})},
        {KeyValues("OutputCollection", {"MCParticles"})}) {}

StatusCode DCH_dNdx_FromTracks::initialize() {
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
    double Rmin = m_geoSvc->getDetector()->constantAsDouble(m_Rmin_parameter_name.value());
    double Rmax = m_geoSvc->getDetector()->constantAsDouble(m_Rmax_parameter_name.value());
    double Zmax = m_geoSvc->getDetector()->constantAsDouble(m_Zmax_parameter_name.value());
    double Zmin;
    // For fwd-bwd symmetric detectors, use negative Zmax as Zmin
    if (m_Zmin_parameter_name.value() == m_Zmax_parameter_name.value()) {
        Zmin = -Zmax;
    } else {
        Zmin = m_geoSvc->getDetector()->constantAsDouble(m_Zmin_parameter_name.value());
    }

    debug() << "Geometry parameters:" << endmsg;
    debug() << "Rmin: " << Rmin << " cm" << endmsg;
    debug() << "Rmax: " << Rmax << " cm" << endmsg;
    debug() << "Zmin: " << Zmin << " cm" << endmsg;
    debug() << "Zmax: " << Zmax << " cm" << endmsg;
    // dd4hep uses cm, while delphes (or more importantly the track parametrisation) uses mm, so we need to convert
    m_delphesTrkUtil->SetDchBoundaries(Rmin*10, Rmax*10, Zmin*10, Zmax*10);
    m_delphesTrkUtil->SetGasMix(m_GasSel.value());

    return StatusCode::SUCCESS;
}

edm4hep::RecDqdxCollection DCH_dNdx_FromTracks::operator()(const edm4hep::TrackMCParticleLinkCollection& input,
                                                            const edm4hep::EventHeaderCollection& header) const {
    edm4hep::RecDqdxCollection outputCollection;

    auto engine_seed = m_uniqueIDSvc->getUniqueID(header, this->name());
    m_engine.seed(engine_seed);

    std::uniform_real_distribution<double> mean_dist(5.0, 20.0);
    double random_mean = mean_dist(m_engine);

    std::poisson_distribution<int> poisson_dist(random_mean);

    debug() << "Processing new Event" << endmsg;
    debug() << "Random mean for this event: " << random_mean << endmsg;


    unsigned int i = 0;
    for (auto link : input) {
        debug() << "Processing track " << i++ << endmsg;

        // Get the track and corresponding MCParticle
        const auto& mc_particle = link.getTo();
        const auto& track = link.getFrom();
        // Check for validity
        if (!mc_particle.isAvailable() || !track.isAvailable()) {
            warning() << "Invalid link found, skipping." << endmsg;
            continue;
        }

        /////////////////////////
        // Cluster Information //
        /////////////////////////
        double momentum = edm4hep::utils::magnitude(mc_particle.getMomentum());
        double mass     = mc_particle.getMass();
        // debug() << "MCParticle momentum: " << momentum << endmsg;
        // debug() << "MCParticle mass: " << mass << endmsg;
        // TODO: safeguards for zero mass and no charge

        double betagamma = momentum/mass;

        ///////////////////////
        // Track Information //
        ///////////////////////
        // Use track state at IP, since this corresponds to delphes and energy loss in tracking system is negligible
        const auto& track_state = track.getTrackStates(1);

        // Convert edm4hep::TrackState to delphes parameters
        // Inverse conversion from https://github.com/key4hep/k4SimDelphes/blob/main/converter/src/DelphesEDM4HepConverter.cc#L532
        // Note the same order of parameters as in delphes: D0, phi, C, Z0, cot(theta)
        // See: https://github.com/delphes/delphes/blob/98f15add056e657e39bda9e32ccd97ef427ce04c/external/TrackCovariance/TrkUtil.cc#L961
        TVectorD delphes_track(5);
        delphes_track[0] = track_state.D0;
        delphes_track[1] = track_state.phi;
        double scale = -2.0;            // delphes uses C instead of omega, scale is used to convert
        delphes_track[2] = track_state.omega / scale;
        delphes_track[3] = track_state.Z0;
        delphes_track[4] = track_state.tanLambda; // tanLambda and cot(theta) are the same thing (see DelphesEDM4HepConverter.cc)

        // Note: track length will be in mm, since this is what delphes uses
        double track_length = m_delphesTrkUtil->TrkLen(delphes_track);
        debug() << "Track length inside chamber: " << track_length << " mm" << endmsg;

        /////////////////////////////////////////////
        // Draw Number of Clusters from Poissonian //
        /////////////////////////////////////////////

        int random_hits = poisson_dist(m_engine);
        // debug() << "RANDOM number of clusters: " << random_hits << endmsg;
        info() << "Track has " << random_hits << " clusters." << endmsg;

        auto dqdx = outputCollection.create();
        edm4hep::Quantity q;
        q.value = static_cast<float>(random_hits);        

        dqdx.setDQdx(q);
        dqdx.setTrack(track);
    }


    return outputCollection;
}
