#pragma once

// k4FWCore
#include "k4FWCore/Transformer.h"

// edm4hep
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"

// edm4hep extension
#include "extension/SenseWireHitCollection.h"

class DCHdigi_v02 final
    : public k4FWCore::Transformer<extension::SenseWireHitCollection(const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&)> {

public:
    DCHdigi_v02(const std::string& name, ISvcLocator* svcLoc);

    extension::SenseWireHitCollection operator()(const edm4hep::SimTrackerHitCollection& input,
                                          const edm4hep::EventHeaderCollection& header) const override;

    StatusCode initialize() override final;

private:
    mutable unsigned int m_event_counter;

};


DECLARE_COMPONENT(DCHdigi_v02);
