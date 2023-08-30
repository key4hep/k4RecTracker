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

/** @class DCHdigitizer
 *
 *  Algorithm for creating digitized drift chamber hits (still based on edm4hep::TrackerHit) from edm4hep::SimTrackerHit.
 *  You have to specify the expected resolution in z and in xy (distance to the wire). The smearing is applied in the wire reference frame.
 *  
 *  @author Brieuc Francois
 *  @date   2023-03
 *
 */

class DCHdigitizer : public GaudiAlgorithm {
public:
  explicit DCHdigitizer(const std::string&, ISvcLocator*);
  virtual ~DCHdigitizer();
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
  DataHandle<edm4hep::TrackerHitCollection> m_output_digi_hits{"outputDigiHits", Gaudi::DataHandle::Writer, this};

  // Detector readout name
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "CDCHHits", "Name of the detector readout"};
  // Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  // Decoder for the cellID
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  // Volume manager to get the physical cell sensitive volume
  dd4hep::VolumeManager m_volman;

  // z position resolution in mm
  FloatProperty m_z_resolution{this, "zResolution", 1.0,
                               "Spatial resolution in the z direction (from reading out the wires at both sides) [mm]"};
  // xy resolution in mm
  FloatProperty m_xy_resolution{this, "xyResolution", 0.1, "Spatial resolution in the xy direction [mm]"};

  // Random Number Service
  IRndmGenSvc* m_randSvc;
  // Gaussian random number generator used for the smearing of the z position
  Rndm::Numbers m_gauss_z;
  // Gaussian random number generator used for the smearing of the xy position
  Rndm::Numbers m_gauss_xy;
};
