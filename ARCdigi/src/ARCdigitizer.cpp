#include "ARCdigitizer.h"

// EDM4HEP
#include "edm4hep/Vector3d.h"

// DD4HEP
#include "DDRec/CellIDPositionConverter.h"

// STL
#include <unordered_map>

DECLARE_COMPONENT(ARCdigitizer)

ARCdigitizer::ARCdigitizer(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {
  declareProperty("inputSimHits", m_input_sim_hits, "Input sim tracker hit collection name");
  declareProperty("outputDigiHits", m_output_digi_hits, "Output digitized tracker hit collection name");
  declareProperty("detectorCompact", m_det_compact, "Path to detector compact");
}

ARCdigitizer::~ARCdigitizer() {}

StatusCode ARCdigitizer::initialize() {
  m_detector = dd4hep::Detector::make_unique(this->name() + "_detector");
  if (m_det_compact.value().empty()) {
    error() << "Detector compact must be provided!" << endmsg;
    return StatusCode::FAILURE;
  }
  m_detector->fromCompact(m_det_compact);
  return StatusCode::SUCCESS;
}

StatusCode ARCdigitizer::execute() {

  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  verbose() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Dictionary to keep track of cell IDs and deposited energies
  std::unordered_map<uint64_t, float> merged_digi_hits;

  // Digitize the sim hits
  for (const auto& input_sim_hit : *input_sim_hits) {
    auto cell = input_sim_hit.getCellID();
    if (merged_digi_hits.find(cell) == merged_digi_hits.end()) merged_digi_hits[cell] = 0.;
    merged_digi_hits[cell] += input_sim_hit.getEDep();
  }

  // Get our cell ID -> position converter
  dd4hep::rec::CellIDPositionConverter converter(*m_detector);

  // Write the digitized hits
  edm4hep::TrackerHitCollection* output_digi_hits = m_output_digi_hits.createAndPut();
  for (auto it = merged_digi_hits.begin(); it != merged_digi_hits.end(); it++) {
    auto output_digi_hit = output_digi_hits->create();
    auto pos = converter.position(it->first);
    output_digi_hit.setPosition(edm4hep::Vector3d(pos.X(), pos.Y(), pos.Z()));
    output_digi_hit.setEDep(it->second);
  }

  return StatusCode::SUCCESS;
}

StatusCode ARCdigitizer::finalize() {
  return StatusCode::SUCCESS;
}
