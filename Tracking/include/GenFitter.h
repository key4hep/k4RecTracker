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
#include "extension/TrackerHit3DCollection.h"
// EDM4HEP extension
#include "extension/DriftChamberDigiCollection.h"
#include "extension/DriftChamberDigiLocalCollection.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
#include "onnxruntime_cxx_api.h"
#include <memory> 
#include <vector> 
#include <string>   
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

class GenFitter : public Gaudi::Algorithm {
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
  virtual StatusCode execute(const EventContext&) const final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:
  /// Pointer to the ONNX enviroment
  std::unique_ptr<Ort::Env> fEnv;
  /// Pointer to the ONNX inference session
  std::unique_ptr<Ort::Session> fSession;
  /// ONNX settings
  Ort::SessionOptions fSessionOptions;
  /// ONNX memory info
  const OrtMemoryInfo* fInfo;
  struct MemoryInfo;
  /// the input names represent the names given to the model
  /// when defining  the model's architecture (if applicable)
  /// they can also be retrieved from model.summary()
  std::vector<const char*> fInames;
  std::vector<const char*> fOnames;

  std::string modelPath={};
  // Input tracker hit collection name
  mutable DataHandle<extension::DriftChamberDigiCollection> m_input_hits_CDC{"inputHits_CDC", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<edm4hep::TrackerHit3DCollection> m_input_hits_VTXIB{"inputHits_VTXIB", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<edm4hep::TrackerHit3DCollection> m_input_hits_VTXD{"inputHits_VTXD", Gaudi::DataHandle::Reader, this};
  mutable DataHandle<edm4hep::TrackerHit3DCollection> m_input_hits_VTXOB{"inputHits_VTXOB", Gaudi::DataHandle::Reader, this};

  // Output track collection name
  mutable DataHandle<extension::TrackCollection> m_output_tracks{"outputTracks", Gaudi::DataHandle::Writer, this};
  mutable DataHandle<extension::TrackerHit3DCollection> m_output_hits{"outputHits", Gaudi::DataHandle::Writer, this};
};
