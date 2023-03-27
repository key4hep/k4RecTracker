#include "VTXdigitizer.h"

DECLARE_COMPONENT(VTXdigitizer)

VTXdigitizer::VTXdigitizer(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {
  declareProperty("inputSimHits", m_input_sim_hits, "Input sim vertex hit collection name");
  declareProperty("outputDigiHits", m_output_digi_hits, "Output digitized vertex hit collection name");
}

VTXdigitizer::~VTXdigitizer() {}

StatusCode VTXdigitizer::initialize() { return StatusCode::SUCCESS; }

StatusCode VTXdigitizer::execute() { 
  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  verbose() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Digitize the sim hits
  edm4hep::TrackerHitCollection* output_digi_hits = m_output_digi_hits.createAndPut();
  for (const auto& input_sim_hit : *input_sim_hits) {
    auto output_digi_hit = output_digi_hits->create();
    output_digi_hit.setEDep(input_sim_hit.getEDep());
    output_digi_hit.setPosition(input_sim_hit.getPosition());
  }
  return StatusCode::SUCCESS;
}

StatusCode VTXdigitizer::finalize() { return StatusCode::SUCCESS; }
