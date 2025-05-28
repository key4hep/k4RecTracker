#include "DCH_dNdx_FromTracks.h"

// Gaudi
#include <GaudiKernel/MsgStream.h>

// EDM4hep
#include "edm4hep/Quantity.h"

DECLARE_COMPONENT(DCH_dNdx_FromTracks)

DCH_dNdx_FromTracks::DCH_dNdx_FromTracks(const std::string& name, ISvcLocator* svcLoc)
    : Transformer(name, svcLoc,
        {KeyValues("InputLinkCollection", {"TrackMCParticleLinks"}),
         KeyValues("HeaderName", {"EventHeader"}),
         KeyValues("InputTrackCollection", {"Tracks"})},
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
                                                            const edm4hep::EventHeaderCollection& header,
                                                            const edm4hep::TrackCollection& tracks) const {
    edm4hep::RecDqdxCollection outputCollection;

    auto engine_seed = m_uniqueIDSvc->getUniqueID(header, this->name());
    m_engine.seed(engine_seed);

    std::uniform_real_distribution<double> mean_dist(5.0, 20.0);
    double random_mean = mean_dist(m_engine);

    std::poisson_distribution<int> poisson_dist(random_mean);

    debug() << "Processing new Event" << endmsg;
    debug() << "Random mean for this event: " << random_mean << endmsg;

    if (input.size() != tracks.size()) {
        warning() << "Mismatch between TrackMCParticleLinkCollection size (" << input.size()
                  << ") and TrackCollection size (" << tracks.size() << "). Processing only the minimum of the two." << endmsg;
    }
    size_t n = std::min(input.size(), tracks.size());

    for (size_t i = 0; i < n; ++i) {
        debug() << "Processing track " << i << endmsg;

        int random_hits = poisson_dist(m_engine);
        debug() << "Random number of hits: " << random_hits << endmsg;
        info() << "Track " << i << " has " << random_hits << " hits." << endmsg;

        auto dqdx = outputCollection.create();
        edm4hep::Quantity q;
        q.type = 0;
        q.value = static_cast<float>(random_hits);
        q.error = std::sqrt(q.value);

        const auto& track = tracks[i];
        debug() << "Using Track with ID " << track.id() << endmsg;

        const auto& link = input[i];
        const auto& mc_particle = link.getTo();
        const auto& link_track = link.getFrom();
        debug() << "Linked MCParticle PDGID: " << mc_particle.getPDG() << endmsg;
        debug() << "Linking Track ID: " << link_track.id() << endmsg;

        if (track!= link_track) {
            warning() << "Not the same tracks" << endmsg;
        } else {
            debug() << "Track and link_track are the same." << endmsg;
        }


        dqdx.setDQdx(q);
        dqdx.setTrack(track);
    }


    return outputCollection;
}
