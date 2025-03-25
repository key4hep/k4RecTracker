#pragma once

// GAUDI
#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

// K4FWCORE
#include "k4FWCore/DataHandle.h"

// EDM4HEP
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHit3DCollection.h"

// DD4HEP
#include "DD4hep/Detector.h"

/** @class ARCdigitizer
 *
 *  Algorithm for creating digitized (meaning 'reconstructed' for now) ARC hits (edm4hep::TrackerHit3D) from Geant4 hits
 * (edm4hep::SimTrackerHit).
 *
 *  @author Brieuc Francois, Matthew Basso
 *  @date   2023-03
 *
 */

class ARCdigitizer : public Gaudi::Algorithm {
public:
  explicit ARCdigitizer(const std::string&, ISvcLocator*);
  virtual ~ARCdigitizer();
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
  // Input sim tracker hit collection name
  mutable DataHandle<edm4hep::SimTrackerHitCollection> m_input_sim_hits{"inputSimHits", Gaudi::DataHandle::Reader,
                                                                        this};
  // Output digitized tracker hit collection name
  mutable DataHandle<edm4hep::TrackerHit3DCollection> m_output_digi_hits{"outputDigiHits", Gaudi::DataHandle::Writer,
                                                                         this};
  // Flat value for SiPM efficiency
  FloatProperty m_flat_SiPM_effi{this, "flatSiPMEfficiency", -1.0,
                                 "Flat value for SiPM quantum efficiency (<0 := disabled)"};
  // Apply the SiPM efficiency to digitized hits instead of simulated hits
  BooleanProperty m_apply_SiPM_effi_to_digi{this, "applySiPMEffiToDigiHits", false,
                                            "Apply the SiPM efficiency to digitized hits instead of simulated hits"};

  // Detector geometry
  dd4hep::Detector* m_detector;
  // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  // Uniform random number generator used for the SiPM quantum efficiency
  Rndm::Numbers m_uniform;
};
