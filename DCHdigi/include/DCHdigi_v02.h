#pragma once

// k4FWCore
#include "k4FWCore/Transformer.h"

// edm4hep
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"

// edm4hep extension
#include "extension/SenseWireHitCollection.h"

// DD4hep
#include "DDSegmentation/BitFieldCoder.h"

// DDRec
#include "DDRec/DCH_info.h"

class DCHdigi_v02 final
    : public k4FWCore::Transformer<extension::SenseWireHitCollection(const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&)> {

public:
    DCHdigi_v02(const std::string& name, ISvcLocator* svcLoc);

    extension::SenseWireHitCollection operator()(const edm4hep::SimTrackerHitCollection& input,
                                          const edm4hep::EventHeaderCollection& header) const override;

    StatusCode initialize() override final;

private:
    mutable unsigned int m_event_counter;

    dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

    // Drift chamber info extension for geometry calculations
    dd4hep::rec::DCH_info* m_dch_info{nullptr};

    /// Convert EDM4hep Vector3d to TVector3
    TVector3 Convert_EDM4hepVector_to_TVector3(const edm4hep::Vector3d& v, double scale) const {
        return {v[0] * scale, v[1] * scale, v[2] * scale};
    };

};


DECLARE_COMPONENT(DCHdigi_v02);
