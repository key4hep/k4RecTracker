#ifndef DCH_DNDX_FROM_TRACKS_H
#define DCH_DNDX_FROM_TRACKS_H

// Gaudi
#include <GaudiKernel/SmartIF.h>
#include <GaudiKernel/ISvcLocator.h>

// k4Interface
#include <k4Interface/IUniqueIDGenSvc.h>

// k4FWCore
#include "k4FWCore/Transformer.h"

// edm4hep
#include "edm4hep/TrackMCParticleLinkCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/RecDqdxCollection.h"

// delphes
#include "TrackCovariance/TrkUtil.h"

// STL
#include <random>

class DCH_dNdx_FromTracks final
    : public k4FWCore::Transformer<edm4hep::RecDqdxCollection(const edm4hep::TrackMCParticleLinkCollection&, const edm4hep::EventHeaderCollection&)> {

public:
    DCH_dNdx_FromTracks(const std::string& name, ISvcLocator* svcLoc);

    edm4hep::RecDqdxCollection operator()(const edm4hep::TrackMCParticleLinkCollection& input,
                                          const edm4hep::EventHeaderCollection& header) const override;

    StatusCode initialize() final;

private:
    SmartIF<IUniqueIDGenSvc> m_uniqueIDSvc{nullptr};
    inline static thread_local std::mt19937_64 m_engine;

    TrkUtil m_delphesTrkUtil;
};

#endif // DCH_DNDX_FROM_TRACKS_H