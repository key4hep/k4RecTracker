#include "DCHdigitizer.h"

// DD4hep
#include "DDRec/Vector3D.h"

DECLARE_COMPONENT(DCHdigitizer)

DCHdigitizer::DCHdigitizer(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {
  declareProperty("inputSimHits", m_input_sim_hits, "Input sim tracker hit collection name");
  declareProperty("outputDigiHits", m_output_digi_hits, "Output digitized tracker hit collection name");
}

DCHdigitizer::~DCHdigitizer() {}

StatusCode DCHdigitizer::initialize() {

  // Initialize random services
  if (service("RndmGenSvc", m_randSvc).isFailure()) {
    error() << "Couldn't get RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_gauss_z.initialize(m_randSvc, Rndm::Gauss(0., m_z_resolution)).isFailure()) {
    error() << "Couldn't initialize RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_gauss_xy.initialize(m_randSvc, Rndm::Gauss(0., m_xy_resolution)).isFailure()) {
    error() << "Couldn't initialize RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

StatusCode DCHdigitizer::execute() {
  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  verbose() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Digitize the sim hits
  edm4hep::TrackerHitCollection* output_digi_hits = m_output_digi_hits.createAndPut();
  for (const auto& input_sim_hit : *input_sim_hits) {
    auto output_digi_hit = output_digi_hits->create();
    output_digi_hit.setEDep(input_sim_hit.getEDep());
    // FIXME temporary solution before to have the distance to the wire
    edm4hep::Vector3d true_global_position_edm = input_sim_hit.getPosition();
    dd4hep::rec::Vector3D true_global_position(true_global_position_edm.x, true_global_position_edm.y, true_global_position_edm.z);
    dd4hep::rec::Vector3D true_global_position_atZ0(true_global_position_edm.x, true_global_position_edm.y, 0.);
    double true_global_radius_atZ0 = true_global_position_atZ0.r();
    std::cout << "Radius " << true_global_radius_atZ0 << std::endl;
    double reco_global_radius_atZ0 = true_global_radius_atZ0 + m_gauss_xy.shoot();
    edm4hep::Vector3d smeared_position();
    output_digi_hit.setPosition(input_sim_hit.getPosition());
  }
  return StatusCode::SUCCESS;
}

StatusCode DCHdigitizer::finalize() { return StatusCode::SUCCESS; }
