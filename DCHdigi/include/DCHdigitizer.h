#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"

// K4FWCORE
#include "k4FWCore/DataHandle.h"

// EDM4HEP
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

/** @class DCHdigitizer
 *
 *  Algorithm for creating digitized (meaning 'reconstructed' for now) drift chamber hits (edm4hep::TrackerHit) from Geant4 hits (edm4hep::SimTrackerHit).
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
};
