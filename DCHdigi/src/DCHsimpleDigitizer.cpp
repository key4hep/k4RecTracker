#include "DCHsimpleDigitizer.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"

// ROOT
#include "Math/Cylindrical3D.h"

DECLARE_COMPONENT(DCHsimpleDigitizer)

DCHsimpleDigitizer::DCHsimpleDigitizer(const std::string& aName, ISvcLocator* aSvcLoc)
    : GaudiAlgorithm(aName, aSvcLoc), m_geoSvc("GeoSvc", "DCHsimpleDigitizer") {
  declareProperty("inputSimHits", m_input_sim_hits, "Input sim tracker hit collection name");
  declareProperty("outputDigiHits", m_output_digi_hits, "Output digitized tracker hit collection name");
}

DCHsimpleDigitizer::~DCHsimpleDigitizer() {}

StatusCode DCHsimpleDigitizer::initialize() {
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
  if (m_geoSvc->lcdd()->readouts().find(m_readoutName) == m_geoSvc->lcdd()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  // set the cellID decoder
  m_decoder = m_geoSvc->lcdd()->readout(m_readoutName).idSpec().decoder();
  // retrieve the volume manager
  m_volman = m_geoSvc->lcdd()->volumeManager();

  return StatusCode::SUCCESS;
}

StatusCode DCHsimpleDigitizer::execute() {
  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  debug() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Digitize the sim hits
  edm4hep::TrackerHitCollection* output_digi_hits = m_output_digi_hits.createAndPut();
  for (const auto& input_sim_hit : *input_sim_hits) {
    auto output_digi_hit = output_digi_hits->create();
    // smear the hit position: need to go in the wire local frame to smear in the direction aligned/perpendicular with the wire for z/distanceToWire, taking e.g. stereo angle into account
    // retrieve the cell detElement
    dd4hep::DDSegmentation::CellID cellID         = input_sim_hit.getCellID();
    auto                           cellDetElement = m_volman.lookupDetElement(cellID);
    // retrieve the wire (in DD4hep 1.23 there is no easy way to access the volume daughters we have to pass by detElements, in later versions volumes can be used)
    const std::string& wireDetElementName =
        Form("superLayer_%d_layer_%d_phi_%d_wire", m_decoder->get(cellID, "superLayer"),
             m_decoder->get(cellID, "layer"), m_decoder->get(cellID, "phi"));
    dd4hep::DetElement wireDetElement = cellDetElement.child(wireDetElementName);
    // get the transformation matrix used to place the wire
    const auto& wireTransformMatrix = wireDetElement.nominal().worldTransformation();
    // Retrieve global position in mm and apply unit transformation (translation matrix is tored in cm)
    double simHitGlobalPosition[3] = {input_sim_hit.getPosition().x * dd4hep::mm,
                                      input_sim_hit.getPosition().y * dd4hep::mm,
                                      input_sim_hit.getPosition().z * dd4hep::mm};
    double simHitLocalPosition[3]  = {0, 0, 0};
    // get the simHit coordinate in cm in the wire reference frame to be able to apply smearing of radius perpendicular to the wire
    wireTransformMatrix.MasterToLocal(simHitGlobalPosition, simHitLocalPosition);
    debug() << "Cell ID string: " << m_decoder->valueString(cellID) << endmsg;
    ;
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
    double smearedDistanceToWire = simHitLocalPositionVector.rho() + m_gauss_xy.shoot() * dd4hep::mm;
    // smear the z position (in local coordinate the z axis is aligned with the wire i.e. it take the stereo angle into account);
    double smearedZ = simHitLocalPositionVector.z() + m_gauss_z.shoot() * dd4hep::mm;
    // build the local position vector of the smeared hit using cylindrical coordinates. When we will have edm4hep::DCHit there will be probably no need
    ROOT::Math::Cylindrical3D digiHitLocalPositionVector(smearedDistanceToWire, smearedZ,
                                                         simHitLocalPositionVector.phi());
    debug() << "Local simHit distance to the wire: " << simHitLocalPositionVector.rho()
            << " Local digiHit distance to the wire: " << smearedDistanceToWire << " in cm" << endmsg;
    debug() << "Local simHit z: " << simHitLocalPositionVector.z()
            << " Local digiHit distance to the wire: " << smearedZ << " in cm" << endmsg;
    debug() << "Local simHit phi: " << simHitLocalPositionVector.phi()
            << " Local digiHit distance to the wire: " << digiHitLocalPositionVector.Phi() << endmsg;
    // go back to the global frame
    double digiHitLocalPosition[3]  = {digiHitLocalPositionVector.x(), digiHitLocalPositionVector.y(),
                                       digiHitLocalPositionVector.z()};
    double digiHitGlobalPosition[3] = {0, 0, 0};
    wireTransformMatrix.LocalToMaster(digiHitLocalPosition, digiHitGlobalPosition);
    // go back to mm
    edm4hep::Vector3d digiHitGlobalPositionVector(digiHitGlobalPosition[0] / dd4hep::mm,
                                                  digiHitGlobalPosition[1] / dd4hep::mm,
                                                  digiHitGlobalPosition[2] / dd4hep::mm);
    output_digi_hit.setPosition(digiHitGlobalPositionVector);
  }
  debug() << "Output Digi Hit collection size: " << output_digi_hits->size() << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode DCHsimpleDigitizer::finalize() { return StatusCode::SUCCESS; }
