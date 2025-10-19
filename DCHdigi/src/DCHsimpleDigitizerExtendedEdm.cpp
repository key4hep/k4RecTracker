#include "DCHsimpleDigitizerExtendedEdm.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"

// ROOT
#include "Math/Cylindrical3D.h"

DECLARE_COMPONENT(DCHsimpleDigitizerExtendedEdm)

DCHsimpleDigitizerExtendedEdm::DCHsimpleDigitizerExtendedEdm(const std::string& aName, ISvcLocator* aSvcLoc)
    : Gaudi::Algorithm(aName, aSvcLoc), m_geoSvc("GeoSvc", "DCHsimpleDigitizerExtendedEdm") {
  declareProperty("inputSimHits", m_input_sim_hits, "Input sim tracker hit collection name");
  declareProperty("outputDigiHits", m_output_digi_hits, "Output digitized tracker hit collection name");
  declareProperty("outputSimDigiAssociation", m_output_sim_digi_association,
                  "Output name for the association between digitized and simulated hit collections");
}

DCHsimpleDigitizerExtendedEdm::~DCHsimpleDigitizerExtendedEdm() {}

StatusCode DCHsimpleDigitizerExtendedEdm::initialize() {
  // Initialize random services
  m_randSvc = service("RndmGenSvc", false);
  if (!m_randSvc) {
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

StatusCode DCHsimpleDigitizerExtendedEdm::execute(const EventContext&) const {
  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  debug() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Prepare a collection for digitized hits in local coordinate (only filled in debug mode)
  extension::DriftChamberDigiCollection* output_digi_hits = m_output_digi_hits.createAndPut();
  extension::DriftChamberDigiLocalCollection* output_digi_local_hits = m_output_digi_local_hits.createAndPut();

  // Prepare a collection for the association between digitized and simulated hit, setting weights to 1
  extension::MCRecoDriftChamberDigiAssociationCollection* digi_sim_associations =
      m_output_sim_digi_association.createAndPut();

  // Prepare collections for the debugging distributions
  auto leftHitSimHitDeltaDistToWire = m_leftHitSimHitDeltaDistToWire.createAndPut();
  auto leftHitSimHitDeltaLocalZ = m_leftHitSimHitDeltaLocalZ.createAndPut();
  auto rightHitSimHitDeltaDistToWire = m_rightHitSimHitDeltaDistToWire.createAndPut();
  auto rightHitSimHitDeltaLocalZ = m_rightHitSimHitDeltaLocalZ.createAndPut();

  // Digitize the sim hits
  for (const auto& input_sim_hit : *input_sim_hits) {
    auto output_digi_hit = output_digi_hits->create();
    // smear the hit position: need to go in the wire local frame to smear in the direction aligned/perpendicular with
    // the wire for z/distanceToWire, taking e.g. stereo angle into account retrieve the cell detElement
    dd4hep::DDSegmentation::CellID cellID = input_sim_hit.getCellID();
    auto cellDetElement = m_volman.lookupDetElement(cellID);
    // retrieve the wire (in DD4hep 1.23 there is no easy way to access the volume daughters we have to pass by
    // detElements, in later versions volumes can be used)
    const std::string& wireDetElementName =
        Form("superLayer_%ld_layer_%ld_phi_%ld_wire", m_decoder->get(cellID, "superLayer"),
             m_decoder->get(cellID, "layer"), m_decoder->get(cellID, "phi"));
    dd4hep::DetElement wireDetElement = cellDetElement.child(wireDetElementName);
    // get the transformation matrix used to place the wire (DD4hep works with cm)
    const auto wireTransformMatrix = wireDetElement.nominal().worldTransformation();
    // Retrieve global position in mm and transform it to cm because the DD4hep translation matrix is stored in cm
    double simHitGlobalPosition[3] = {input_sim_hit.getPosition().x * dd4hep::mm,
                                      input_sim_hit.getPosition().y * dd4hep::mm,
                                      input_sim_hit.getPosition().z * dd4hep::mm};
    double simHitLocalPosition[3] = {0, 0, 0};
    // get the simHit coordinate in cm in the wire reference frame to be able to apply smearing of radius perpendicular
    // to the wire
    wireTransformMatrix.MasterToLocal(simHitGlobalPosition, simHitLocalPosition);
    debug() << "Cell ID string: " << m_decoder->valueString(cellID) << endmsg;
    debug() << "Global simHit x " << simHitGlobalPosition[0] << " --> Local simHit x " << simHitLocalPosition[0]
            << " in cm" << endmsg;
    debug() << "Global simHit y " << simHitGlobalPosition[1] << " --> Local simHit y " << simHitLocalPosition[1]
            << " in cm" << endmsg;
    debug() << "Global simHit z " << simHitGlobalPosition[2] << " --> Local simHit z " << simHitLocalPosition[2]
            << " in cm" << endmsg;
    // build a vector to easily apply smearing of distance to the wire, going back to mm
    dd4hep::rec::Vector3D simHitLocalPositionVector(
        simHitLocalPosition[0] / dd4hep::mm, simHitLocalPosition[1] / dd4hep::mm, simHitLocalPosition[2] / dd4hep::mm);
    // get the smeared distance to the wire (cylindrical coordinate as the smearing should be perpendicular to the wire)
    debug() << "Original distance to wire: " << simHitLocalPositionVector.rho() << " mm" << endmsg;
    double smearedDistanceToWire = simHitLocalPositionVector.rho() + m_gauss_xy.shoot();
    while (smearedDistanceToWire < 0) {
      debug() << "Negative smearedDistanceToWire (" << smearedDistanceToWire << ") shooting another random number"
              << endmsg;
      smearedDistanceToWire = simHitLocalPositionVector.rho() + m_gauss_xy.shoot();
    }
    debug() << "Smeared distance to wire: " << smearedDistanceToWire << " mm " << endmsg;
    // smear the z position (in local coordinate the z axis is aligned with the wire i.e. it take the stereo angle into
    // account);
    double smearedZ = simHitLocalPositionVector.z() + m_gauss_z.shoot();
    // NB: here we assume the hit is radially in the middle of the cell
    // create the left and right hit local position (in cm again because of the transform matrix)
    double leftHitLocalPosition[3] = {-1 * smearedDistanceToWire * dd4hep::mm, 0, smearedZ * dd4hep::mm};
    double rightHitLocalPosition[3] = {smearedDistanceToWire * dd4hep::mm, 0, smearedZ * dd4hep::mm};
    // transform the left and right hit local position in global coordinate (still cm here)
    double leftHitGlobalPosition[3] = {0, 0, 0};
    double rightHitGlobalPosition[3] = {0, 0, 0};
    wireTransformMatrix.LocalToMaster(leftHitLocalPosition, leftHitGlobalPosition);
    wireTransformMatrix.LocalToMaster(rightHitLocalPosition, rightHitGlobalPosition);
    // std::cout << (leftHitGlobalPosition[2]==rightHitGlobalPosition[2]) <<  std::endl; // FIXME why are left and right
    // global z coordinates the same?
    debug() << "Global leftHit x " << leftHitGlobalPosition[0] << " --> Local leftHit x " << leftHitLocalPosition[0]
            << " in cm" << endmsg;
    debug() << "Global leftHit y " << leftHitGlobalPosition[1] << " --> Local leftHit y " << leftHitLocalPosition[1]
            << " in cm" << endmsg;
    debug() << "Global leftHit z " << leftHitGlobalPosition[2] << " --> Local leftHit z " << leftHitLocalPosition[2]
            << " in cm" << endmsg;
    debug() << "Global rightHit x " << rightHitGlobalPosition[0] << " --> Local rightHit x " << rightHitLocalPosition[0]
            << " in cm" << endmsg;
    debug() << "Global rightHit y " << rightHitGlobalPosition[1] << " --> Local rightHit y " << rightHitLocalPosition[1]
            << " in cm" << endmsg;
    debug() << "Global rightHit z " << rightHitGlobalPosition[2] << " --> Local rightHit z " << rightHitLocalPosition[2]
            << " in cm" << endmsg;
    // fill the output DriftChamberDigi (making sure we are back in mm)
    output_digi_hit.setCellID(cellID);
    edm4hep::Vector3d leftHitGlobalPositionVector =
        edm4hep::Vector3d(leftHitGlobalPosition[0] / dd4hep::mm, leftHitGlobalPosition[1] / dd4hep::mm,
                          leftHitGlobalPosition[2] / dd4hep::mm);
    edm4hep::Vector3d rightHitGlobalPositionVector =
        edm4hep::Vector3d(rightHitGlobalPosition[0] / dd4hep::mm, rightHitGlobalPosition[1] / dd4hep::mm,
                          rightHitGlobalPosition[2] / dd4hep::mm);
    output_digi_hit.setLeftPosition(leftHitGlobalPositionVector);
    output_digi_hit.setRightPosition(rightHitGlobalPositionVector);
    output_digi_hit.setTime(input_sim_hit.getTime()); // will apply smearing when we know more from R&D teams
    // output_digi_hit.setEDep(input_sim_hit.getEDep()); // will enable this when we know more from R&D teams

    // create the association between digitized and simulated hit
    auto digi_sim_association = digi_sim_associations->create();
    digi_sim_association.setDigi(output_digi_hit);
    digi_sim_association.setSim(input_sim_hit);

    // if required, populate debugging distributions
    if (m_debugMode) {
      // Fill the local coordinate digi hit
      auto output_digi_local_hit = output_digi_local_hits->create();
      output_digi_local_hit.setCellID(cellID);
      output_digi_local_hit.setDistanceToWire(smearedDistanceToWire);
      output_digi_local_hit.setZPositionAlongWire(smearedZ);
      // produce a vector instead of using directly smearedDistanceToWire or smearedZ for a more thorough testing
      dd4hep::rec::Vector3D leftHitLocalPositionVector(leftHitLocalPosition[0] / dd4hep::mm,
                                                       leftHitLocalPosition[1] / dd4hep::mm,
                                                       leftHitLocalPosition[2] / dd4hep::mm);
      dd4hep::rec::Vector3D rightHitLocalPositionVector(rightHitLocalPosition[0] / dd4hep::mm,
                                                        rightHitLocalPosition[1] / dd4hep::mm,
                                                        rightHitLocalPosition[2] / dd4hep::mm);
      leftHitSimHitDeltaDistToWire->push_back(leftHitLocalPositionVector.rho() - simHitLocalPositionVector.rho());
      leftHitSimHitDeltaLocalZ->push_back(leftHitLocalPositionVector.z() - simHitLocalPositionVector.z());
      rightHitSimHitDeltaDistToWire->push_back(rightHitLocalPositionVector.rho() - simHitLocalPositionVector.rho());
      rightHitSimHitDeltaLocalZ->push_back(rightHitLocalPositionVector.z() - simHitLocalPositionVector.z());
    }
  }
  debug() << "Output Digi Hit collection size: " << output_digi_hits->size() << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode DCHsimpleDigitizerExtendedEdm::finalize() { return StatusCode::SUCCESS; }
