#pragma once

// GAUDI
#include "Gaudi/Property.h"
// #include "GaudiAlg/GaudiAlgorithm.h"
#include "Gaudi/Algorithm.h"
// K4FWCORE
#include "k4FWCore/DataHandle.h"

// EDM4HEP
#if __has_include("edm4hep/TrackerHit3DCollection.h")
#include "edm4hep/TrackerHit3DCollection.h"
#else
#include "edm4hep/TrackerHit3D.h"
#include "extension/TrackerHit3D.h"
namespace edm4hep {
  using TrackerHit3DCollection = edm4hep::TrackerHitCollection;
}  // namespace edm4hep
#endif

#include "extension/TrackCollection.h"
#include "extension/TrackCollection.h"

// GENFIT
//#include "WireMeasurement.h"

/** @class Evaluation_efficiency
 *
 *  
 *  
 *  @author Andrea De Vita, Maria Dolores Garcia, Brieuc Francois
 *  @date   2024-06
 *
 */

class Evaluation_efficiency : public Gaudi::Algorithm {
public:
  explicit Evaluation_efficiency(const std::string&, ISvcLocator*);
  virtual ~Evaluation_efficiency();
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

  // Output track collection name
  mutable DataHandle<extension::TrackCollection> m_input_tracks{"inputTracks", Gaudi::DataHandle::Writer, this};
 
};
