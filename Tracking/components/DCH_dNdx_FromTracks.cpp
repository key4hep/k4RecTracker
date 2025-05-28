#include "DCH_dNdx_FromTracks.h"

// Gaudi
#include <GaudiKernel/MsgStream.h>

// EDM4hep
#include "edm4hep/Quantity.h"
#include "edm4hep/utils/vector_utils.h"

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
        // Basic debug info
        debug() << "Linked MCParticle PDGID: " << mc_particle.getPDG() << endmsg;
        debug() << "Linking Track ID: " << track.id() << endmsg;

        double momentum = edm4hep::utils::magnitude(mc_particle.getMomentum());
        double mass     = mc_particle.getMass();
        debug() << "MCParticle momentum: " << momentum << endmsg;
        debug() << "MCParticle mass: " << mass << endmsg;
        // TODO: safeguards for zero mass and no charge

        double betagamma = momentum/mass;



        int random_hits = poisson_dist(m_engine);
        debug() << "Random number of hits: " << random_hits << endmsg;
        info() << "Track has " << random_hits << " hits." << endmsg;

        auto dqdx = outputCollection.create();
        edm4hep::Quantity q;
        q.type = 0;
        q.value = static_cast<float>(random_hits);
        q.error = std::sqrt(q.value);

        

        dqdx.setDQdx(q);
        dqdx.setTrack(track);
    }


    return outputCollection;
}
