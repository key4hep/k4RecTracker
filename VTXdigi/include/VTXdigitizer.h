#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

// K4FWCORE
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"

// EDM4HEP
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#if __has_include("edm4hep/TrackerHit3DCollection.h")
#include "edm4hep/TrackerHit3DCollection.h"
#else
#include "edm4hep/TrackerHitCollection.h"
namespace edm4hep {
  using TrackerHit3DCollection = edm4hep::TrackerHitCollection;
}  // namespace edm4hep
#endif

// DD4HEP
#include "DD4hep/Detector.h"  // for dd4hep::VolumeManager
#include "DDRec/Vector3D.h"
#include "DDRec/SurfaceManager.h"

#include "DDSegmentation/BitFieldCoder.h"

/** @class VTXdigitizer
 *
 *  Algorithm for creating digitized (meaning 'reconstructed' for now) vertex detector hits (edm4hep::TrackerHit3D) from Geant4 hits (edm4hep::SimTrackerHit).
 *  
 *  @author Brieuc Francois
 *  @date   2023-03
 *
 */

class VTXdigitizer : public Gaudi::Algorithm {
public:
  explicit VTXdigitizer(const std::string&, ISvcLocator*);
  virtual ~VTXdigitizer();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Execute.
   *   @return status code
   */
  virtual StatusCode execute(const EventContext&) const final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:
  // Input sim vertex hit collection name
  mutable DataHandle<edm4hep::SimTrackerHitCollection> m_input_sim_hits{"inputSimHits", Gaudi::DataHandle::Reader, this};
  // Output digitized vertex hit collection name
  mutable DataHandle<edm4hep::TrackerHit3DCollection> m_output_digi_hits{"outputDigiHits", Gaudi::DataHandle::Writer, this};

  // Detector name
  Gaudi::Property<std::string> m_detectorName{this, "detectorName", "Vertex", "Name of the detector (default: Vertex)"};
  // Detector readout names
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "VertexBarrelCollection", "Name of the detector readout"};
  // Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  // Decoder for the cellID
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  // Volume manager to get the physical cell sensitive volume
  dd4hep::VolumeManager m_volman;

  // x resolution in um
  FloatProperty m_x_resolution{this, "xResolution", 0.1, "Spatial resolution in the x direction [um] (r-phi direction in barrel, z direction in disks)"};

  // y resolution in um
  FloatProperty m_y_resolution{this, "yResolution", 0.1, "Spatial resolution in the y direction [um] (r direction in barrel, r-phi direction in disks)"};

  // t resolution in ns
  FloatProperty m_t_resolution{this, "tResolution", 0.1, "Time resolution [ns]"};

  // Surface manager used to project hits onto sensitive surface with forceHitsOntoSurface argument
  const dd4hep::rec::SurfaceMap* _map;

  // Option to force hits onto sensitive surface
  BooleanProperty m_forceHitsOntoSurface{this, "forceHitsOntoSurface", false, "Project hits onto the surface in case they are not yet on the surface (default: false"};

  // Random Number Service
  IRndmGenSvc* m_randSvc;

  // Gaussian random number generator used for smearing
  Rndm::Numbers m_gauss_x;
  Rndm::Numbers m_gauss_y;
  Rndm::Numbers m_gauss_time;
};
