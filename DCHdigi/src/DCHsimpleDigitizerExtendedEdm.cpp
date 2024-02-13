#include "DCHsimpleDigitizerExtendedEdm.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"

// ROOT
#include "Math/Cylindrical3D.h"

DECLARE_COMPONENT(DCHsimpleDigitizerExtendedEdm)

DCHsimpleDigitizerExtendedEdm::DCHsimpleDigitizerExtendedEdm(const std::string& aName, ISvcLocator* aSvcLoc)
    : GaudiAlgorithm(aName, aSvcLoc), m_geoSvc("GeoSvc", "DCHsimpleDigitizerExtendedEdm") {
  declareProperty("inputSimHits", m_input_sim_hits, "Input sim tracker hit collection name");
  declareProperty("outputDigiHits", m_output_digi_hits, "Output digitized tracker hit collection name");
}

DCHsimpleDigitizerExtendedEdm::~DCHsimpleDigitizerExtendedEdm() {}

StatusCode DCHsimpleDigitizerExtendedEdm::initialize() {
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

  // check if readout exists
  if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  // set the cellID decoder
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  // retrieve the volume manager
  m_volman = m_geoSvc->getDetector()->volumeManager();

  return StatusCode::SUCCESS;
}

StatusCode DCHsimpleDigitizerExtendedEdm::execute() {
  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  debug() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Digitize the sim hits
  extension::DriftChamberDigiCollection* output_digi_hits = m_output_digi_hits.createAndPut();
  for (const auto& input_sim_hit : *input_sim_hits) {
    auto output_digi_hit = output_digi_hits->create();
    // smear the hit position: need to go in the wire local frame to smear in the direction aligned/perpendicular with the wire for z/distanceToWire, taking e.g. stereo angle into account
    // retrieve the cell detElement
    dd4hep::DDSegmentation::CellID cellID         = input_sim_hit.getCellID();
    auto                           cellDetElement = m_volman.lookupDetElement(cellID);
    // retrieve the wire (in DD4hep 1.23 there is no easy way to access the volume daughters we have to pass by detElements, in later versions volumes can be used)
    const std::string& wireDetElementName =
        Form("superLayer_%ld_layer_%ld_phi_%ld_wire", m_decoder->get(cellID, "superLayer"),
             m_decoder->get(cellID, "layer"), m_decoder->get(cellID, "phi"));
    dd4hep::DetElement wireDetElement = cellDetElement.child(wireDetElementName);
    // get the transformation matrix used to place the wire
    const auto& wireTransformMatrix = wireDetElement.nominal().worldTransformation();
    // Retrieve global position in mm and apply unit transformation (translation matrix is stored in cm)
    double simHitGlobalPosition[3] = {input_sim_hit.getPosition().x * dd4hep::mm,
                                      input_sim_hit.getPosition().y * dd4hep::mm,
                                      input_sim_hit.getPosition().z * dd4hep::mm};
    double simHitLocalPosition[3]  = {0, 0, 0};
    // get the simHit coordinate in cm in the wire reference frame to be able to apply smearing of radius perpendicular to the wire
    wireTransformMatrix.MasterToLocal(simHitGlobalPosition, simHitLocalPosition);
    debug() << "Cell ID string: " << m_decoder->valueString(cellID) << endmsg;
    debug() << "Global simHit x " << simHitGlobalPosition[0] << " --> Local simHit x " << simHitLocalPosition[0]
            << " in cm" << endmsg;
    debug() << "Global simHit y " << simHitGlobalPosition[1] << " --> Local simHit y " << simHitLocalPosition[1]
            << " in cm" << endmsg;
    debug() << "Global simHit z " << simHitGlobalPosition[2] << " --> Local simHit z " << simHitLocalPosition[2]
            << " in cm" << endmsg;
    // build a vector to easily apply smearing of distance to the wire
    dd4hep::rec::Vector3D simHitLocalPositionVector(simHitLocalPosition[0], simHitLocalPosition[1],
                                                    simHitLocalPosition[2]);
    // get the smeared distance to the wire (cylindrical coordinate as the smearing should be perpendicular to the wire)
    debug() << "Original distance to wire: " << simHitLocalPositionVector.rho() << endmsg;
    double smearedDistanceToWire = simHitLocalPositionVector.rho() + m_gauss_xy.shoot() * dd4hep::mm;
    while(smearedDistanceToWire < 0){
      debug() << "Negative smearedDistanceToWire (" << smearedDistanceToWire << ") shooting another random number" << endmsg;
      smearedDistanceToWire = simHitLocalPositionVector.rho() + m_gauss_xy.shoot() * dd4hep::mm;
    }
    debug() << "Smeared distance to wire: " << smearedDistanceToWire << endmsg;
    // smear the z position (in local coordinate the z axis is aligned with the wire i.e. it take the stereo angle into account);
    double smearedZ = simHitLocalPositionVector.z() + m_gauss_z.shoot() * dd4hep::mm;
    double leftHitLocalPosition[3]  = {-1 * smearedDistanceToWire, 0, smearedZ};
    double rightHitLocalPosition[3]  = {smearedDistanceToWire, 0, smearedZ};
    double leftHitGlobalPosition[3]  = {0, 0, 0};
    double rightHitGlobalPosition[3]  = {0, 0, 0};
    wireTransformMatrix.LocalToMaster(leftHitLocalPosition, leftHitGlobalPosition);
    wireTransformMatrix.LocalToMaster(rightHitLocalPosition, rightHitGlobalPosition);
    // fill the output DriftChamberDigi (making sure we are back in mm)
    output_digi_hit.setCellID(cellID);
    output_digi_hit.setLeftPosition(edm4hep::Vector3d({leftHitGlobalPosition[0] / dd4hep::mm, leftHitGlobalPosition[1] / dd4hep::mm, leftHitGlobalPosition[2] / dd4hep::mm}));
    output_digi_hit.setRightPosition(edm4hep::Vector3d({rightHitGlobalPosition[0] / dd4hep::mm, rightHitGlobalPosition[1] / dd4hep::mm, rightHitGlobalPosition[2] / dd4hep::mm}));
  }
  debug() << "Output Digi Hit collection size: " << output_digi_hits->size() << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode DCHsimpleDigitizerExtendedEdm::finalize() { return StatusCode::SUCCESS; }
