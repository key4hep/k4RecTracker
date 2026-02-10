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
  class ChargeCollector_SinglePixel;
}

struct VTXdigi_Modular final : k4FWCore::MultiTransformer <std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> (const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&)> {

  VTXdigi_Modular(const std::string& name, ISvcLocator* svcLoc);
  
  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> operator() (const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const override;

private: 

  // friend class VTXdigi_tools::IChargeCollector;
  friend class VTXdigi_tools::ChargeCollector_SinglePixel;

  /* ---- Initialization & finalization functions ---- */

  void InitServicesAndGeometry();
  void InitLayersAndSensors();
  void InitHistograms();

  /* ---- Core algorithm functions ---- */

  bool CheckEventSetup(const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const;

  bool CheckSimhitLayer(const edm4hep::SimTrackerHit& simHit) const;
  
  /* -- Properties -- */

  const std::string m_undefinedString = "UNDEFINED";

  Gaudi::Property<std::string> m_subDetName{this, "SubDetectorName", m_undefinedString, "Name of the subdetector (eg. \"Vertex\")"};
  Gaudi::Property<std::string> m_subDetChildName{this, "SubDetectorChildName", m_undefinedString, "Name of the subdetector child (eg. \"VertexBarrel\"), if applicable. If undefined, the subdetector itself is assumed to contain layers as children."};

  Gaudi::Property<std::string> m_geoServiceName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
  Gaudi::Property<std::string> m_encodingStringVariable{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};

  /* -- Properties mainlyrelated to the main event loop -- */

  Gaudi::Property<std::vector<int>> m_layersToDigitize{this, "LayersToDigitize", {}, "Which layers to digitize (0-indexed). If empty, all layers are digitized."};

  /* -- Properties and members related to the various charge collection algorithms-- */

  Gaudi::Property<std::string> m_chargeCollectionMethod{this, "ChargeCollectionMethod", "Drift", "Method used for charge collection: \"Fast\", \"Drift\", \"LookupTable\", etc."};
  
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

  int m_layerCount = 0; 
  std::array<size_t, 2> m_pixelCount = {0, 0}; 
  std::array<float, 2> m_pixelPitch = {0.0f, 0.0f};
  float m_sensorThickness = 0.0f;
  std::array<float, 2> m_sensorLength = {0.0f, 0.0f};
  TGeoRotation m_sensorNormalRotation = TGeoRotation("sensorNormalRotation"); // rotation to rotate the sensor local coordinate system. Initialised to unit matrix. The same information was previously stored in m_localNormalVectorDir.

  /* -- Counters -- */

  mutable Gaudi::Accumulators::Counter<> m_counter_eventsRead{this, "Events read"};
  mutable Gaudi::Accumulators::Counter<> m_counter_eventsRejected_noSimHits{this, "Events rejected (no simHits)"};
  mutable Gaudi::Accumulators::Counter<> m_counter_eventsAccepted{this, "Events accepted"};
  mutable Gaudi::Accumulators::Counter<> m_counter_simHitsRead{this, "SimTrackerHits read"};
  mutable Gaudi::Accumulators::Counter<> m_counter_simHitsRejected_LayerNotToBeDigitized{this, "SimTrackerHits rejected (layer not to be digitized)"};
  mutable Gaudi::Accumulators::Counter<> m_counter_simHitsAccepted{this, "SimTrackerHits accepted"};

}; // class VTXdigi_Modular





