#include "Evaluation_efficiency.h"
#include <iostream>
#include <sstream>
#include <vector>




DECLARE_COMPONENT(Evaluation_efficiency)

GenFitter::Evaluation_efficiency(const std::string& aName, ISvcLocator* aSvcLoc) : Gaudi::Algorithm(aName, aSvcLoc) {
  declareProperty("inputTracks", m_input_tracks, "Input track collection name");
}

Evaluation_efficiency::~Evaluation_efficiency() {}

StatusCode Evaluation_efficiency::initialize() { 
  return StatusCode::SUCCESS;}


StatusCode Evaluation_efficiency::execute(const EventContext&) const {
  // Get the input collection with tracker hits
  

  return StatusCode::SUCCESS;
}

StatusCode GenFitter::finalize() { return StatusCode::SUCCESS; }

