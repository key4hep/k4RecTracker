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

    TrkUtil* m_delphesTrkUtil;

    // Drift Chamber (geometry) parameters
    Gaudi::Property<double> m_Zmin{
        this,
        "Zmin",
        {-2250.0},
        "smallest Z value inside drift chamber in mm.)"};
    Gaudi::Property<double> m_Zmax{
        this,
        "Zmax",
        {2250.0},
        "largest Z value inside drift chamber in mm.)"};
    Gaudi::Property<double> m_Rmin{
        this,
        "Rmin",
        {349.8},
        "inner radius of the drift chamber in mm.)"};
    Gaudi::Property<double> m_Rmax{
        this,
        "Rmax",
        {2015.0},
        "outer radius of the drift chamber in mm.)"};
    Gaudi::Property<int> m_GasSel{
        this,
        "GasSel",
        {0},
        "Gas selection: 0: He(90%)-Isobutane(10%), 1: pure He, 2: Ar(50%)-Ethane(50%), 3: pure Ar."};
};

#endif // DCH_DNDX_FROM_TRACKS_H