#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/SmartIF.h"

// K4FWCORE & podio
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "podio/UserDataCollection.h"

// EDM4HEP
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "edm4hep/Vector3d.h"

// DD4HEP
#include "DD4hep/Detector.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Vector3D.h"
#include "DDSegmentation/BitFieldCoder.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

/** @class MuRWELLDigitizer
 *
 *  Cluster-based digitizer for the IDEA muon system (uRWELL chambers).
 *
 *  Each Geant4 SimTrackerHit is first assigned to the (y,z) strip-cell
 *  given by its cellID. Cells that share a chamber/slice and fall within
 *  the same time window are then grouped into clusters:
 *
 *    1. Pick the unused fired cell with the highest summed energy (seed).
 *    2. Greedily attach the highest-energy fired cell that is (N°2 -1)-neighbour-
 *       adjacent to any cell already in the cluster, where N is the cluster size,
 *        *as long as* doing so
 *       keeps the number of distinct Y-strip indices and the number of
 *       distinct Z-strip indices both within MaxClusterSize. This matches
 *       the test-beam cluster-size definition (strips per view).
 *    3. Emit one digi hit at the energy-weighted centroid plus optional
 *       Gaussian smearing, apply per-cluster efficiency, and repeat.
 *
 *  Inputs:  inputSimHits (SimTrackerHitCollection)
 *  Outputs: outputDigiHits, outputSimDigiLink, and five validation collections
 *           (simDigiDifferenceX/Y/Z, clusterSizeY/Z), all returned by operator().
 *
 *  @author Mahmoud Althakeel
 *  @date   2026-06
 */

// All collections produced per event, in the order they are declared in the constructor.
using MuRWELLDigitizerOutputs = std::tuple<edm4hep::TrackerHitPlaneCollection,             // outputDigiHits
                                           edm4hep::TrackerHitSimTrackerHitLinkCollection, // outputSimDigiLink
                                           podio::UserDataCollection<double>,              // simDigiDifferenceX
                                           podio::UserDataCollection<double>,              // simDigiDifferenceY
                                           podio::UserDataCollection<double>,              // simDigiDifferenceZ
                                           podio::UserDataCollection<int>,                 // clusterSizeY
                                           podio::UserDataCollection<int>>;                // clusterSizeZ

class MuRWELLDigitizer final
    : public k4FWCore::MultiTransformer<MuRWELLDigitizerOutputs(const edm4hep::SimTrackerHitCollection&)> {
public:
  explicit MuRWELLDigitizer(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  MuRWELLDigitizerOutputs operator()(const edm4hep::SimTrackerHitCollection& simHits) const override;

private:
  // ---- Geometry --------------------------------------------------------
  Gaudi::Property<std::string> m_detectorName{this, "detectorName", "Muon-System", "DD4hep detector name"};
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "MuonSystemCollection",
                                             "Name of the detector readout"};
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
  SmartIF<IGeoSvc> m_geoSvc;
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder = nullptr;
  dd4hep::VolumeManager m_volman;
  const dd4hep::rec::SurfaceMap* m_surfMap = nullptr; // populated only when forceHitsOntoSurface is true

  // ---- Smearing / efficiency / clustering knobs ------------------------
  // uRWELL is a planar 2D readout: the chamber-local x is the gas-gap normal
  // (not measured) and only the two
  // in-plane surface axes (u = local y, v = local z) are smeared.
  FloatProperty m_u_resolution{this, "uResolution", 0.4,
                               "Spatial resolution along surface u (= local y, in-plane) [mm]"};
  FloatProperty m_v_resolution{this, "vResolution", 0.4,
                               "Spatial resolution along surface v (= local z, in-plane) [mm]"};
  FloatProperty m_t_resolution{
      this, "tResolution", -1.0,
      "Time resolution [ns]; <=0 disables time smearing (digi time = earliest contributing SimHit time)"};
  FloatProperty m_efficiency{this, "efficiency", 0.98, "Per-cluster detection efficiency"};
  FloatProperty m_timeWindow{this, "timeWindow", 25.0, "Time window for grouping SimHits into the same cluster [ns]"};
  Gaudi::Property<unsigned int> m_maxClusterSize{
      this, "maxClusterSize", 3, "Maximum number of distinct strip indices per view (Y and Z) inside one cluster"};
  BooleanProperty m_forceHitsOntoSurface{
      this, "forceHitsOntoSurface", false,
      "If true, project smeared digi positions that fall outside the chamber sensitive surface back onto it"};
  Gaudi::Property<unsigned int> m_printFrequency{
      this, "printFrequency", 100, "Print one INFO line every N events with the in/out hit count (0 = never)"};

  // ---- Random services -------------------------------------------------
  SmartIF<IRndmGenSvc> m_randSvc;
  Rndm::Numbers m_gauss_u; // smear along surface u (local y)
  Rndm::Numbers m_gauss_v; // smear along surface v (local z)
  Rndm::Numbers m_gauss_t; // smear time if tResolution > 0
  Rndm::Numbers m_flat;
};
