#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"

// K4FWCORE
#include "k4FWCore/DataHandle.h"

// EDM4HEP
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"

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
};
