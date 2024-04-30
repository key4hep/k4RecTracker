#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"

// K4FWCORE
#include "k4FWCore/DataHandle.h"

// EDM4HEP
#if __has_include("edm4hep/TrackerHit3DCollection.h")
#include "edm4hep/TrackerHit3DCollection.h"
#else
#include "edm4hep/TrackerHitCollection.h"
namespace edm4hep {
  using TrackerHit3DCollection = edm4hep::TrackerHitCollection;
}  // namespace edm4hep
#endif

#include "edm4hep/TrackCollection.h"
// EDM4HEP extension
#include "extension/DriftChamberDigiCollection.h"
#include "extension/DriftChamberDigiLocalCollection.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"

// GENFIT
//#include "WireMeasurement.h"

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
  DataHandle<extension::DriftChamberDigiCollection> m_input_hits_CDC{"inputHits_CDC", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHit3DCollection> m_input_hits_VTXIB{"inputHits_VTXIB", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHit3DCollection> m_input_hits_VTXD{"inputHits_VTXD", Gaudi::DataHandle::Reader, this};
  DataHandle<edm4hep::TrackerHit3DCollection> m_input_hits_VTXOB{"inputHits_VTXOB", Gaudi::DataHandle::Reader, this};

  // Output track collection name
  DataHandle<edm4hep::TrackCollection> m_output_tracks{"outputTracks", Gaudi::DataHandle::Writer, this};
  // Transient genfit measurements used internally by genfit to run the tracking
  //std::vector<genfit::WireMeasurement> m_wire_measurements;
};
