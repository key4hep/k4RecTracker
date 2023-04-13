#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"

// K4FWCORE
#include "k4FWCore/DataHandle.h"

// EDM4HEP
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"

// GENFIT
#include "WireMeasurement.h"

/** @class GenFitter
 *
 *  
 *  
 *  @author Maria Dolores Garcia, Brieuc Francois
 *  @date   2023-03
 *
 */

class GenFitter : public GaudiAlgorithm {
public:
  explicit GenFitter(const std::string&, ISvcLocator*);
  virtual ~GenFitter();
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
  // Input tracker hit collection name
  DataHandle<edm4hep::TrackerHitCollection> m_input_hits{"inputHits", Gaudi::DataHandle::Reader, this};
  // Output track collection name
  DataHandle<edm4hep::TrackCollection> m_output_tracks{"outputTracks", Gaudi::DataHandle::Writer, this};
  // Transient genfit measurements used internally by genfit to run the tracking 
  std::vector<genfit::WireMeasurement> m_wire_measurements;


};
