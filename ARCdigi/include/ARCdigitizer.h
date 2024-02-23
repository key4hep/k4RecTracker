#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

// K4FWCORE
#include "k4FWCore/DataHandle.h"

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
#include "DD4hep/Detector.h"

/** @class ARCdigitizer
 *
 *  Algorithm for creating digitized (meaning 'reconstructed' for now) ARC hits (edm4hep::TrackerHit3D) from Geant4 hits (edm4hep::SimTrackerHit).
 *  
 *  @author Brieuc Francois, Matthew Basso
 *  @date   2023-03
 *
 */

class ARCdigitizer : public GaudiAlgorithm {
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
  // Flat value for SiPM efficiency
  FloatProperty m_flat_SiPM_effi{this, "flatSiPMEfficiency", -1.0, "Flat value for SiPM quantum efficiency (<0 := disabled)"};
  // Apply the SiPM efficiency to digitized hits instead of simulated hits
  BooleanProperty m_apply_SiPM_effi_to_digi{this, "applySiPMEffiToDigiHits", false, "Apply the SiPM efficiency to digitized hits instead of simulated hits"};

  // Detector geometry
  dd4hep::Detector* m_detector;
  // Random Number Service
  IRndmGenSvc* m_randSvc;
  // Uniform random number generator used for the SiPM quantum efficiency
  Rndm::Numbers m_uniform;
};
