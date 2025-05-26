// Gaudi

// Gaudi services
#include <k4Interface/IUniqueIDGenSvc.h>

// k4FWCore
#include "k4FWCore/Transformer.h"

// edm4hep
#include "edm4hep/TrackCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/RecDqdxCollection.h"

// STL
#include <random>


struct DCH_dNdx_FromTracks final
        : k4FWCore::Transformer<edm4hep::RecDqdxCollection(const edm4hep::TrackCollection&, const edm4hep::EventHeaderCollection&)> {
    DCH_dNdx_FromTracks(const std::string& name, ISvcLocator* svcLoc)
        : Transformer(name, svcLoc, 
            {KeyValues("InputCollection", {"Tracks"}),
             KeyValues("HeaderName", {"EventHeader"})},
            {KeyValues("OutputCollection", {"MCParticles"})}) {}

    edm4hep::RecDqdxCollection operator()(const edm4hep::TrackCollection& input, const edm4hep::EventHeaderCollection& header) const override {
        edm4hep::RecDqdxCollection outputCollection;

        auto engine_seed = m_uniqueIDSvc->getUniqueID(header, this->name());
        m_engine.seed(engine_seed);

        // Generate a random mean in range [5.0, 20.0]
        std::uniform_real_distribution<double> mean_dist(5.0, 20.0);
        double random_mean = mean_dist(m_engine);

        std::poisson_distribution<int> poisson_dist(random_mean);

        debug() << "Processing new Event" << endmsg;
        debug() << "Random mean for this event: " << random_mean << endmsg;
        for (size_t i = 0; i < input.size(); ++i) {
            debug() << "Processing track " << i << endmsg;
            int random_hits = poisson_dist(m_engine);
            debug() << "Random number of hits: " << random_hits << endmsg;
            auto dqdx = outputCollection.create();

            edm4hep::Quantity q;
            q.type = 0;
            q.value = static_cast<float>(random_hits);
            q.error = std::sqrt(q.value);
            dqdx.setDQdx(q);
            dqdx.setTrack(input[i]);

        }
        return outputCollection;
    }

    StatusCode initialize() final {
        m_uniqueIDSvc = service("UniqueIDGenSvc");
        if (!m_uniqueIDSvc) {
            error() << "Unable to locate the UniqueIDGenSvc" << endmsg;
            return StatusCode::FAILURE;
        }
        return StatusCode::SUCCESS;
    }

    private:
        SmartIF<IUniqueIDGenSvc> m_uniqueIDSvc{nullptr};
        inline static thread_local std::mt19937_64 m_engine;
};

DECLARE_COMPONENT(DCH_dNdx_FromTracks)
