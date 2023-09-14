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
#include "edm4hep/TrackerHitCollection.h"

// DD4HEP
#include "DD4hep/Detector.h"  // for dd4hep::VolumeManager
#include "DDSegmentation/BitFieldCoder.h"

/** @class VTXdigitizer
 *
 *  Algorithm for creating digitized (meaning 'reconstructed' for now) vertex detector hits (edm4hep::TrackerHit) from Geant4 hits (edm4hep::SimTrackerHit).
 *  
 *  @author Brieuc Francois
 *  @date   2023-03
 *
 */

class VTXdigitizer : public GaudiAlgorithm {
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
  virtual StatusCode execute() final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:
  // Input sim vertex hit collection name
  DataHandle<edm4hep::SimTrackerHitCollection> m_input_sim_hits{"inputSimHits", Gaudi::DataHandle::Reader, this};
  // Output digitized vertex hit collection name
  DataHandle<edm4hep::TrackerHitCollection> m_output_digi_hits{"outputDigiHits", Gaudi::DataHandle::Writer, this};

  // Detector readout names
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "VTXDCollection", "Name of the detector readout"};
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

  // Random Number Service
  IRndmGenSvc* m_randSvc;
  // Gaussian random number generator used for smearing
  Rndm::Numbers m_gauss_x;
  Rndm::Numbers m_gauss_y;
  Rndm::Numbers m_gauss_t;
};
