#include "ARCdigitizer.h"

// EDM4HEP
#include "edm4hep/Vector3d.h"

// DD4HEP
#include "DDRec/CellIDPositionConverter.h"

// STL
#include <unordered_map>
#include <utility>

DECLARE_COMPONENT(ARCdigitizer)

ARCdigitizer::ARCdigitizer(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {
  declareProperty("inputSimHits", m_input_sim_hits, "Input sim tracker hit collection name");
  declareProperty("outputDigiHits", m_output_digi_hits, "Output digitized tracker hit collection name");
  declareProperty("detectorCompact", m_det_compact, "Path to detector compact");
  declareProperty("flatSiPMEfficiency", m_flat_SiPM_effi, "Flat value for SiPM quantum efficiency (<0 := disabled)");
  declareProperty("applySiPMEffiToDigiHits", m_apply_SiPM_effi_to_digi,
                  "Apply the SiPM efficiency to digitized hits instead of simulated hits");
}

ARCdigitizer::~ARCdigitizer() {}

StatusCode ARCdigitizer::initialize() {
  m_detector = dd4hep::Detector::make_unique(this->name() + "_detector");
  if (m_det_compact.value().empty()) {
    error() << "Detector compact must be provided!" << endmsg;
    return StatusCode::FAILURE;
  }
  // Initialize random service
  if (service("RndmGenSvc", m_randSvc).isFailure()) {
    error() << "Couldn't get RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_uniform.initialize(m_randSvc, Rndm::Flat(0.0, 1.0)).isFailure()) {
    error() << "Couldn't initialize RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }
  // Sanity check on efficiency cut
  if (m_flat_SiPM_effi > 1.0) {
    error() << "Flat SiPM efficiency cannot exceed 1!" << endmsg;
    return StatusCode::FAILURE;
  }
  m_detector->fromCompact(m_det_compact);
  return StatusCode::SUCCESS;
}

StatusCode ARCdigitizer::execute() {
  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  verbose() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Dictionary to keep track of cell IDs and (summed) deposited energies / (earliest) arrival times
  std::unordered_map<uint64_t, std::pair<float, float>> merged_digi_hits;

  // Digitize the sim hits
  for (const auto& input_sim_hit : *input_sim_hits) {
    // Throw away simulated hits based on flat SiPM efficiency
    if (!m_apply_SiPM_effi_to_digi && m_flat_SiPM_effi >= 0.0 && m_uniform.shoot() > m_flat_SiPM_effi)
      continue;
    auto cell = input_sim_hit.getCellID();
    if (merged_digi_hits.find(cell) == merged_digi_hits.end())
      merged_digi_hits[cell] = std::pair<float, float>(0.0, -1.0);
    merged_digi_hits[cell].first += input_sim_hit.getEDep();
    if (merged_digi_hits[cell].second < 0.0 || input_sim_hit.getTime() < merged_digi_hits[cell].second)
      merged_digi_hits[cell].second = input_sim_hit.getTime();
  }

  // Get our cell ID -> position converter
  dd4hep::rec::CellIDPositionConverter converter(*m_detector);

  // Write the digitized hits
  edm4hep::TrackerHitCollection* output_digi_hits = m_output_digi_hits.createAndPut();
  for (auto it = merged_digi_hits.begin(); it != merged_digi_hits.end(); it++) {
    // Throw away digitized hits based on flat SiPM efficiency
    if (m_apply_SiPM_effi_to_digi && m_flat_SiPM_effi >= 0.0 && m_uniform.shoot() > m_flat_SiPM_effi)
      continue;
    auto output_digi_hit = output_digi_hits->create();
    auto pos             = converter.position(it->first);
    output_digi_hit.setCellID(it->first);
    output_digi_hit.setPosition(edm4hep::Vector3d(pos.X(), pos.Y(), pos.Z()));
    output_digi_hit.setEDep(it->second.first);
    output_digi_hit.setTime(it->second.second);
  }

  verbose() << "Output Digi Hit collection size: " << output_digi_hits->size() << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode ARCdigitizer::finalize() { return StatusCode::SUCCESS; }
