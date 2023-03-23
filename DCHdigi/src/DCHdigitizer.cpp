#include "DCHdigitizer.h"

DECLARE_COMPONENT(DCHdigitizer)

DCHdigitizer::DCHdigitizer(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {
  declareProperty("inputSimHits", m_input_sim_hits, "Input sim tracker hit collection name");
  declareProperty("outputRecHits", m_output_rec_hits, "Output reconstructed tracker hit collection name");
}

DCHdigitizer::~DCHdigitizer() {}

StatusCode DCHdigitizer::initialize() { return StatusCode::SUCCESS; }

StatusCode DCHdigitizer::execute() { 
  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  verbose() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Digitize the sim hits
  edm4hep::TrackerHitCollection* output_rec_hits = m_output_rec_hits.createAndPut();
  for (const auto& input_sim_hit : *input_sim_hits) {
    auto output_rec_hit = output_rec_hits->create();
    output_rec_hit.setEDep(input_sim_hit.getEDep());
    output_rec_hit.setPosition(input_sim_hit.getPosition());
  }
  return StatusCode::SUCCESS;
}

StatusCode DCHdigitizer::finalize() { return StatusCode::SUCCESS; }
