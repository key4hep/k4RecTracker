#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

// K4FWCORE
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"

// EDM4HEP & PODIO
#include "edm4hep/SimTrackerHitCollection.h"
#include "podio/UserDataCollection.h"

// EDM4HEP extension
#include "extension/DriftChamberDigiCollection.h"
#include "extension/DriftChamberDigiLocalCollection.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"

// DD4HEP
#include "DD4hep/Detector.h"  // for dd4hep::VolumeManager
#include "DDSegmentation/BitFieldCoder.h"

/** @class DCHsimpleDigitizerExtendedEdm
 *
 *  Algorithm for creating digitized drift chamber hits (extension::DriftChamberDigi) from edm4hep::SimTrackerHit.
 *  You have to specify the expected resolution in z and in xy (distance to the wire). The smearing is applied in the wire reference frame.
 *  
 *  @author Brieuc Francois
 *  @date   2023-03
 *
 */

class DCHsimpleDigitizerExtendedEdm : public GaudiAlgorithm {
public:
  explicit DCHsimpleDigitizerExtendedEdm(const std::string&, ISvcLocator*);
  virtual ~DCHsimpleDigitizerExtendedEdm();
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
  DataHandle<extension::DriftChamberDigiCollection> m_output_digi_hits{"outputDigiHits", Gaudi::DataHandle::Writer, this};
  // Output association between digitized and simulated hit collections
  DataHandle<extension::MCRecoDriftChamberDigiAssociationCollection> m_output_sim_digi_association{"outputSimDigiAssociation", Gaudi::DataHandle::Writer, this};
  // Output digitized tracker hit in local coordinates collection name. Only filled in debug mode
  DataHandle<extension::DriftChamberDigiLocalCollection> m_output_digi_local_hits{"outputDigiLocalHits", Gaudi::DataHandle::Writer, this};

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
  // Flag to produce debugging distributions
  Gaudi::Property<bool> m_debugMode{this, "debugMode", false, "Flag to produce debugging distributions"};
  // Declaration of debugging distributions
  DataHandle<podio::UserDataCollection<double>> m_leftHitSimHitDeltaDistToWire{"leftHitSimHitDeltaDistToWire", Gaudi::DataHandle::Writer, this}; // mm
  DataHandle<podio::UserDataCollection<double>> m_leftHitSimHitDeltaLocalZ{"leftHitSimHitDeltaLocalZ", Gaudi::DataHandle::Writer, this}; // mm
  DataHandle<podio::UserDataCollection<double>> m_rightHitSimHitDeltaDistToWire{"rightHitSimHitDeltaDistToWire", Gaudi::DataHandle::Writer, this}; // mm
  DataHandle<podio::UserDataCollection<double>> m_rightHitSimHitDeltaLocalZ{"rightHitSimHitDeltaLocalZ", Gaudi::DataHandle::Writer, this}; // mm

  // Random Number Service
  IRndmGenSvc* m_randSvc;
  // Gaussian random number generator used for the smearing of the z position
  Rndm::Numbers m_gauss_z;
  // Gaussian random number generator used for the smearing of the xy position
  Rndm::Numbers m_gauss_xy;
};
