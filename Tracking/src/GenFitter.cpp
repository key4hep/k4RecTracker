#include "GenFitter.h"

DECLARE_COMPONENT(GenFitter)

GenFitter::GenFitter(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {
  declareProperty("inputHits", m_input_hits, "Input tracker hit collection name");
  declareProperty("outputTracks", m_output_tracks, "Output track collection name");
}

GenFitter::~GenFitter() {}

StatusCode GenFitter::initialize() { return StatusCode::SUCCESS; }

StatusCode GenFitter::execute() {
  // Get the input collection with tracker hits
  const edm4hep::TrackerHit3DCollection* input_hits = m_input_hits.get();
  verbose() << "Input Hit collection size: " << input_hits->size() << endmsg;

  // Convert edm4hep::TrackerHitCollection to Genfit measurements
  //for (const auto& input_hit : *input_hits) {
  //  m_wire_measurements.push_back(genfit::WireMeasurement());
  //}

  // Produce dummy tracks
  edm4hep::TrackCollection* output_tracks = m_output_tracks.createAndPut();
  auto                      output_track  = output_tracks->create();
  output_track.setChi2(1.);
  output_track.setNdf(1);
  return StatusCode::SUCCESS;
}

StatusCode GenFitter::finalize() { return StatusCode::SUCCESS; }
