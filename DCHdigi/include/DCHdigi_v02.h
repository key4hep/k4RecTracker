#pragma once

// Gaudi
#include <GaudiKernel/SmartIF.h>
#include <GaudiKernel/ISvcLocator.h>

// k4FWCore
#include "k4FWCore/Transformer.h"

// k4Interface
#include <k4Interface/IUniqueIDGenSvc.h>
#include <k4Interface/IGeoSvc.h>

// edm4hep
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"

// edm4hep extension
#include "extension/SenseWireHitCollection.h"

// DD4hep
#include "DDSegmentation/BitFieldCoder.h"

// DDRec
#include "DDRec/DCH_info.h"

// delphes
#include "TrackCovariance/TrkUtil.h"

class DCHdigi_v02 final
    : public k4FWCore::Transformer<extension::SenseWireHitCollection(const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&)> {

public:
    DCHdigi_v02(const std::string& name, ISvcLocator* svcLoc);

    extension::SenseWireHitCollection operator()(const edm4hep::SimTrackerHitCollection& input,
                                          const edm4hep::EventHeaderCollection& header) const override;

    StatusCode initialize() override final;

private:
    SmartIF<IUniqueIDGenSvc> m_uniqueIDSvc{nullptr};
    Gaudi::Property<std::string> m_uidSvcName{this, "uidSvcName", "UniqueIDGenSvc", "The name of the UniqueIDGenSvc instance"};

    SmartIF<IGeoSvc> m_geoSvc{nullptr};
    Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};

    mutable unsigned int m_event_counter;

    TrkUtil m_delphesTrkUtil;

    dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

    // Detector name
    Gaudi::Property<std::string> m_dch_name{this, "DCH_name", "DCH_v2", "Name of the Drift Chamber detector"};

    // Drift chamber info extension for geometry calculations
    dd4hep::rec::DCH_info* m_dch_info{nullptr};

    // z resolution in mm
    Gaudi::Property<double> m_z_resolution{
        this, 
        "zResolution_mm", 
        1.0,
        "Spatial resolution in the z direction (along the wire) in mm."};
    // xy resolution in mm
    Gaudi::Property<double> m_xy_resolution{
        this, 
        "xyResolution_mm", 
        0.1,
        "Spatial resolution in the xy direction in mm."};

    // Gas mixture in the chamber
    Gaudi::Property<int> m_GasSel{
        this,
        "GasSel",
        {0},
        "Gas selection: 0: He(90%)-Isobutane(10%), 1: pure He, 2: Ar(50%)-Ethane(50%), 3: pure Ar."};


    /// Convert EDM4hep Vector3d to TVector3
    TVector3 Convert_EDM4hepVector_to_TVector3(const edm4hep::Vector3d& v, double scale) const {
        return {v[0] * scale, v[1] * scale, v[2] * scale};
    };
    /// Convert TVector3 to EDM4hep Vector3d
    edm4hep::Vector3d Convert_TVector3_to_EDM4hepVector(const TVector3& v, double scale) const {
        return {v.x() * scale, v.y() * scale, v.z() * scale};
    };

    // /// Function to calculate the drift time from the distance to the wire
    // double get_drift_time(double distance_to_wire_mm) const;

    // /// Function to calculate the time it takes for the signal to travel from the wire to the readout electronics
    // double get_signal_travel_time(double distance_to_readout_mm) const;

};


DECLARE_COMPONENT(DCHdigi_v02);
