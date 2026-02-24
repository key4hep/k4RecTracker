#pragma once

#include "IChargeCollector.h" 
#include "VTXdigi_tools.h"

#include "Gaudi/Property.h"

#include "Gaudi/Accumulators/RootHistogram.h" // added by Jona
#include "Gaudi/Accumulators/Histogram.h" // added by Jona

#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/EventHeaderCollection.h" // added by Jona
#include "edm4hep/TrackerHitPlaneCollection.h" // added by Jona
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h" // added by Jona

// K4FWCORE
#include "k4Interface/IGeoSvc.h"
#include "k4FWCore/Transformer.h" // added by Jona
#include "k4Interface/IUniqueIDGenSvc.h" // added by Jona

// DD4HEP
#include "DDRec/SurfaceManager.h"

// ROOT
#include "TRandom2.h"

// C++ std
#include <string> // added by Jona
#include <vector>
#include <cmath> // for std::fmod
#include <typeinfo>

// for debugging-csv
#include <iostream>
#include <fstream>  // for std::ofstream
#include <sstream>

#include "GAUDI_VERSION.h"

/** @class VTXdigi_Modular
 *
 * Creates trackerHits from simHits. Produces clusters from simHits, outputs either the cluster centre or all hits in the cluster as digitized hits.
 *
 *  @author Jona Dilg, Armin Ilg
 *  @date   2025-09
 */

/* -- Forward declarations -- */
namespace VTXdigi_tools {
  class IChargeCollector;
}

struct VTXdigi_Modular final : k4FWCore::MultiTransformer <std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> (const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&)> {

  VTXdigi_Modular(const std::string& name, ISvcLocator* svcLoc);
  
  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> operator() (const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const override;

private: 

  /* ---- Initialization & finalization functions ---- */

  void InitServicesAndGeometry();
  void InitLayersAndSensors();
  void InitHistograms();

  /* ---- Core algorithm functions ---- */

  bool CheckEventSetup(const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const;

  bool CheckSimhitLayer(const edm4hep::SimTrackerHit& simHit) const;

  std::vector<VTXdigi_tools::Cluster> Clusterize(const VTXdigi_tools::HitMap& hitMap) const;

  void CreateDigiHits(edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitLinks, const dd4hep::DDSegmentation::CellID& cellID, const TGeoHMatrix& trafoMatrix, const std::vector<VTXdigi_tools::Cluster>& clusters) const;
  
  void FillHistograms_perSimHit(const VTXdigi_tools::SimHitWrapper& hit) const;
  void FillHistograms_perPixel(const dd4hep::DDSegmentation::CellID& cellID, const VTXdigi_tools::Pixel& pix, const std::pair<float, float> clusterPos_local) const;
  void FillHistograms_perDigiHit(const std::unordered_set<std::shared_ptr<const edm4hep::SimTrackerHit>>& simTrackerHits, const edm4hep::TrackerHitPlane& digiHit, const TGeoHMatrix& trafoMatrix, const int clusterSize) const;
  
  /* -- Properties -- */

  const std::string m_undefinedString = "UNDEFINED";

  Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", m_undefinedString, "Name of the subdetector (eg. \"Vertex\")"};
  Gaudi::Property<std::string> m_subDetChildName{this, "SubDetectorChildName", m_undefinedString, "Name of the subdetector child (eg. \"VertexBarrel\"), if applicable. If undefined, the subdetector itself is assumed to contain layers as children."};

  Gaudi::Property<std::string> m_geoServiceName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
  Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
  Gaudi::Property<bool> m_clusterize{this, "Clusterize", true, "Whether to clusterize hits or not. If false, each pixel hit is output as a separate trackerHit. If true, the cluster centre and total charge are output."};


  /* -- Properties mainlyrelated to the main event loop -- */

  Gaudi::Property<std::vector<int>> m_layers{this, "Layers", {}, "Which layers to run on (0-indexed). If empty, all layers are run."};

  /* -- Properties and members related to the various charge collection algorithms-- */

  Gaudi::Property<std::string> m_chargeCollectionMethod{this, "ChargeCollectionMethod", "Drift", "Method used for charge collection: \"Fast\", \"Drift\", \"LookupTable\", etc."};
  Gaudi::Property<float> m_depletedRegionDepthCenter{this, "DepletedRegionDepthCenter", 0.0f, "Depth of the depleted region center for charge collection (in mm), wrt to the pixel center at 0 mm. Used for the digitised hit position."};
  
  
  Gaudi::Property<bool> m_debugHistograms{this, "DebugHistograms", false, "Flag to create and fill debug histograms. Not recommended for multithreading, might lead to crashes. Default is false."};
  Gaudi::Property<int> m_infoPrintInterval{this, "InfoPrintInterval", 100, "Interval for printing information during processing."};

  /* LUT */
  Gaudi::Property<std::string> m_LUT_FileName{this, "LookupTableFile", "", "File to load the lookup table from. Must be given if ChargeCollectionMethod is set to \"LookupTable\"."};
  Gaudi::Property<float> m_LUT_stepLength{this, "LookupTableSegmentStepLength", 0.0001f, "Length of the segments that a particle path through the sensor is split into. The deposited charge is distributed evenly over the segments, and each segments charge is distributed according to the in-pixel bin the segment center falls into. In mm. Defaults to 0.0004 mm."};
  
  /* -- Services, geometry variables -- */
  
  SmartIF<IUniqueIDGenSvc> m_uidService;
  SmartIF<IGeoSvc> m_geoService;
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_cellIdDecoder;
  const dd4hep::Detector* m_detector = nullptr;
  const dd4hep::rec::SurfaceMap* m_surfaceMap; // map from cellID (unsigned long, without segmentation bits) to simSurface (dd4hep::rec::ISurface*)
  dd4hep::VolumeManager m_volumeManager; // volume manager to get the physical cell sensitive volume
  dd4hep::DetElement m_subDetector; // subdetector DetElement. contains layers as children
  
  /* -- Constants -- */
  
  const float m_chargePerkeV = 273.97f; // number of electron-hole pairs created per keV of deposited energy in silicon. eh-pair ~ 3.65 eV
  
  /* -- Member variables -- */

  std::unique_ptr<VTXdigi_tools::IChargeCollector> m_chargeCollector = nullptr;

  std::pair<size_t, size_t> m_pixelCount = {0, 0}; 
  std::pair<float, float> m_pixelPitch = {0.0f, 0.0f};
  float m_sensorThickness = 0.0f;
  std::pair<float, float> m_sensorLength = {0.0f, 0.0f};
  TGeoRotation m_sensorNormalRotation = TGeoRotation("sensorNormalRotation"); // rotation to rotate the sensor local coordinate system. Initialised to unit matrix. The same information was previously stored in m_localNormalVectorDir.

  /* -- Counters -- */

  mutable Gaudi::Accumulators::Counter<> m_counter_eventsRead{this, "Events read"};
  mutable Gaudi::Accumulators::Counter<> m_counter_eventsRejected_noSimHits{this, "Events rejected (no simHits)"};
  mutable Gaudi::Accumulators::Counter<> m_counter_eventsAccepted{this, "Events accepted"};
  mutable Gaudi::Accumulators::Counter<> m_counter_simHitsRead{this, "SimTrackerHits read"};
  mutable Gaudi::Accumulators::Counter<> m_counter_simHitsRejected_LayerNotToBeDigitized{this, "SimTrackerHits rejected (layer not to be digitized)"};
  mutable Gaudi::Accumulators::Counter<> m_counter_simHitsAccepted{this, "SimTrackerHits accepted"};

  /* -- Histograms -- */

  enum { 
    hist1d_simHitE, 
    hist1d_simHitCharge,
    hist1d_simHitMomentum_keV,
    hist1d_simHitMomentum_MeV,
    hist1d_simHitMomentum_GeV,
    hist1d_digiHitCharge,
    hist1d_simHitPDG,
    hist1d_digiHitsPerSimHit,
    hist1d_clusterSize,
    hist1d_clusterSize_createdInGenerator,
    hist1d_clusterSize_createdInSimulation,
    hist1d_residualU, 
    hist1d_residualV, 
    hist1d_residualW, 
    hist1d_residualR, 
    hist1d_pixelDistToClusterCenterU,
    hist1d_pixelDistToClusterCenterV,
    hist1dArrayLen
  }; 
  mutable std::unordered_map<
    int, // layer number
    std::array<
      std::unique_ptr<
        Gaudi::Accumulators::StaticHistogram<
          1, 
          Gaudi::Accumulators::atomicity::full,
          float
        >
      >,
      hist1dArrayLen
    >
  > m_hist1d;

  enum { 
    hist1dglobal_pathLength,
    hist1dglobal_pathLength_Geant4,
    hist1dglobal_pathLength_ratio,
    hist1dglobalArrayLen
  };
  std::array<
    std::unique_ptr<
      Gaudi::Accumulators::StaticHistogram<
        1, 
        Gaudi::Accumulators::atomicity::full,
        float
      >
    >,
    hist1dglobalArrayLen
  > m_hist1dglobal;

  // enum{
  //   histProfile1d_clusterSize_vs_hit_z,
  //   histProfile1d_clusterSize_vs_module_z,
  //   histProfile1d_clusterSize_vs_hit_cosTheta,
  //   histProfile1dArrayLen };
  // mutable std::vector<
  //   std::array<
  //     std::unique_ptr<
  //       Gaudi::Accumulators::StaticProfileHistogram<
  //         1,
  //         Gaudi::Accumulators::atomicity::full,
  //         float
  //       >
  //     >,
  //     histProfile1dArrayLen
  //   >
  // > m_histProfile1d;  

  enum { 
    hist2d_hitMap_simHits,
    hist2d_hitMap_pixelHits, 
    hist2d_clusterSize_vs_hit_z,
    hist2d_clusterSize_vs_hit_z_createdInGenerator,
    hist2d_clusterSize_vs_hit_z_createdInSim,
    hist2dArrayLen
  };
  mutable std::unordered_map<
    int, // layer number
    std::array<
      std::unique_ptr<
        Gaudi::Accumulators::StaticHistogram<
          2,
          Gaudi::Accumulators::atomicity::full,
          float
        >
      >,
      hist2dArrayLen
    >
  > m_hist2d;

public:

  inline std::array<float, 3> SensorDimensions() const { return {m_sensorLength.first, m_sensorLength.second, m_sensorThickness}; }

  inline std::pair<float, float> PixelPitch() const { return m_pixelPitch; }

  inline std::pair<size_t, size_t> PixelCount() const { return m_pixelCount; }

  /** @brief Increment histograms. To be called from the charge collector, once per simHit */
  void FillHistograms_fromChargeCollector_perSimHit(const float pathLength, const float pathLength_Geant4) const;

  inline std::string LutFileName() const { return m_LUT_FileName; }
  inline float LutStepLength() const { return m_LUT_stepLength; }

}; // class VTXdigi_Modular

