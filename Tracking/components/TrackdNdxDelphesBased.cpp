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
    debug() << "Rmin: " << Rmin/dd4hep::cm << " cm" << endmsg;
    debug() << "Rmax: " << Rmax/dd4hep::cm << " cm" << endmsg;
    debug() << "Zmin: " << Zmin/dd4hep::cm << " cm" << endmsg;
    debug() << "Zmax: " << Zmax/dd4hep::cm << " cm" << endmsg;

    // dd4hep uses cm, while delphes (or more importantly the track parametrisation) uses mm, so we need to convert
    m_delphesTrkUtil->SetDchBoundaries(Rmin/dd4hep::mm, Rmax/dd4hep::mm, Zmin/dd4hep::mm, Zmax/dd4hep::mm);

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
    const double dummy_value = -999 * (1/dd4hep::m);
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
                      << dummy_value/(1/dd4hep::m) << " clusters/m" << endmsg;
            success = false;
        }

        // Get number of clusters per length from delphes (output is in m^-1)
        double nclusters_per_meter = m_delphesTrkUtil->Nclusters(betagamma, m_GasSel.value()) * (1/dd4hep::m);
        debug() << "Number of clusters per meter: " << nclusters_per_meter / (1/dd4hep::m) << endmsg;

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

        // Note: track length will be in mm, since this is what delphes and the track parametrisation use
        double track_length = m_delphesTrkUtil->TrkLen(delphes_track)*dd4hep::mm;
        // Check if track length calculation was successful
        if (track_length < std::numeric_limits<double>::epsilon()) {
            warning() << "Delphes track length calculation returned 0.0, dN/dx will be set to dummy value: " 
                      << dummy_value/(1/dd4hep::m) << " clusters/m" << endmsg;
            success = false;
        }
        debug() << "Track length inside full chamber: " << track_length/dd4hep::mm << " mm" << endmsg;
        // Apply fill factor to track length
        track_length *= m_fill_factor.value();
        debug() << "Track length after applying fill factor: " << track_length/dd4hep::mm << " mm" << endmsg;

        /////////////////////////////////////////////
        // Draw Number of Clusters from Poissonian //
        /////////////////////////////////////////////

        double nclusters_mean = nclusters_per_meter * track_length;
        std::poisson_distribution<int> poisson_dist(nclusters_mean);

        int n_cluster = poisson_dist(random_engine);
        debug() << "Track has " << n_cluster << " clusters." << endmsg;

        double dNdx_value = (success)? (n_cluster/track_length) : dummy_value;
        debug() << "dNdx value: " << dNdx_value / (1/dd4hep::m) << " (clusters/m)" << endmsg;

        auto dqdx = outputCollection.create();
        edm4hep::Quantity q;
        q.value = static_cast<float>(dNdx_value/(1/dd4hep::m));

        dqdx.setDQdx(q);
        dqdx.setTrack(track);
    }


    return outputCollection;
}
