#pragma once

// Gaudi
#include <GaudiKernel/SmartIF.h>
#include <GaudiKernel/ISvcLocator.h>

// k4Interface
#include <k4Interface/IUniqueIDGenSvc.h>
#include <k4Interface/IGeoSvc.h>

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

class TrackdNdxDelphesBased final
    : public k4FWCore::Transformer<edm4hep::RecDqdxCollection(const edm4hep::TrackMCParticleLinkCollection&, const edm4hep::EventHeaderCollection&)> {

public:
    TrackdNdxDelphesBased(const std::string& name, ISvcLocator* svcLoc);

    edm4hep::RecDqdxCollection operator()(const edm4hep::TrackMCParticleLinkCollection& input,
                                          const edm4hep::EventHeaderCollection& header) const override;

    StatusCode initialize() final;

private:
    SmartIF<IUniqueIDGenSvc> m_uniqueIDSvc{nullptr};
    SmartIF<IGeoSvc> m_geoSvc{nullptr};
    inline static thread_local std::mt19937_64 m_engine;

    TrkUtil* m_delphesTrkUtil;

    // Drift Chamber (geometry) parameters
    Gaudi::Property<std::string> m_Zmax_parameter_name{
        this,
        "ZmaxParameterName",
        {"DCH_gas_Lhalf"},
        "Name of XML file parameter describing the z extent of the active volume (max value in +z direction)"};
    Gaudi::Property<std::string> m_Zmin_parameter_name{
        this,
        "ZminParameterName",
        {"DCH_gas_Lhalf"},
        "Name of XML file parameter describing the z extent of the active volume (min value in -z direction). For forward-backward symmetric detectors, use same name as for ZmaxParameterName. This value is then automatically converted into the negative value."};
    Gaudi::Property<std::string> m_Rmin_parameter_name{
        this,
        "RminParameterName",
        {"DCH_gas_inner_cyl_R"},
        "Name of XML file parameter describing the inner radius of the active volume."};
    Gaudi::Property<std::string> m_Rmax_parameter_name{
        this,
        "RmaxParameterName",
        {"DCH_gas_outer_cyl_R"},
        "Name of XML file parameter describing the outer radius of the active volume."};
    Gaudi::Property<int> m_GasSel{
        this,
        "GasSel",
        {0},
        "Gas selection: 0: He(90%)-Isobutane(10%), 1: pure He, 2: Ar(50%)-Ethane(50%), 3: pure Ar."};
    Gaudi::Property<double> m_fill_factor{
        this,
        "FillFactor",
        {1.0},
        "Factor (between 0 and 1) describing the fraction of the detector volume that is active (eg. for Straw Tube Tracker the factor is significantly below 1 due to gaps between tubes). The factor is used to scale the calculated track length and thus the corresponding number of clusters."};
};