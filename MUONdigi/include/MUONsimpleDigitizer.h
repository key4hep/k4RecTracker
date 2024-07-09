#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
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

/** @class MUONsimpleDigitizer
 *
 *  Algorithm for creating digitized Muon system hits (still based on edm4hep::TrackerHit3D) from edm4hep::SimTrackerHit.
 *  You have to specify the expected resolution in z and in xy.
 *  
 *  @author Mahmoud Ali
 *  @date   2023-09
 *
 */

class MUONsimpleDigitizer : public GaudiAlgorithm {
public:
  explicit MUONsimpleDigitizer(const std::string&, ISvcLocator*);
  virtual ~MUONsimpleDigitizer();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Execute.
   *   @return status code
   */
  virtual StatusCode execute() final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:
  // Input sim tracker hit collection name
  DataHandle<edm4hep::SimTrackerHitCollection> m_input_sim_hits{"inputSimHits", Gaudi::DataHandle::Reader, this};
  // Output digitized tracker hit collection name
  DataHandle<edm4hep::TrackerHit3DCollection> m_output_digi_hits{"outputDigiHits", Gaudi::DataHandle::Writer, this};

  // Detector readout name
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "MuonSystemCollection", "Name of the detector readout"};
  // Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  // Decoder for the cellID
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  // Volume manager to get the physical cell sensitive volume
  dd4hep::VolumeManager m_volman;

  // z position resolution in mm
  FloatProperty m_x_resolution{this, "xResolution", 1.0,
                               "Spatial resolution in the x direction [mm]"};
  // xy resolution in mm
  FloatProperty m_y_resolution{this, "yResolution", 1.0, "Spatial resolution in the y direction [mm]"};

  // Random Number Service
  IRndmGenSvc* m_randSvc;
  // Gaussian random number generator used for the smearing of the x position
  Rndm::Numbers m_gauss_x;
  // Gaussian random number generator used for the smearing of the y position
  Rndm::Numbers m_gauss_y;
};
