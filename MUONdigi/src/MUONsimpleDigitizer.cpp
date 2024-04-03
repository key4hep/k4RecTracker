#include "MUONsimpleDigitizer.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"

// ROOT
//#include "Math/Cylindrical3D.h"

DECLARE_COMPONENT(MUONsimpleDigitizer)

MUONsimpleDigitizer::MUONsimpleDigitizer(const std::string& aName, ISvcLocator* aSvcLoc)
    : GaudiAlgorithm(aName, aSvcLoc), m_geoSvc("GeoSvc", "MUONsimpleDigitizer") {
  declareProperty("inputSimHits", m_input_sim_hits, "Input sim tracker hit collection name");
  declareProperty("outputDigiHits", m_output_digi_hits, "Output digitized tracker hit collection name");
}

MUONsimpleDigitizer::~MUONsimpleDigitizer() {}

StatusCode MUONsimpleDigitizer::initialize() {
  // Initialize random services
  if (service("RndmGenSvc", m_randSvc).isFailure()) {
    error() << "Couldn't get RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_gauss_x.initialize(m_randSvc, Rndm::Gauss(0., m_x_resolution)).isFailure()) {
    error() << "Couldn't initialize RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_gauss_y.initialize(m_randSvc, Rndm::Gauss(0., m_y_resolution)).isFailure()) {
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

StatusCode MUONsimpleDigitizer::execute() {
  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  debug() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Digitize the sim hits
  edm4hep::TrackerHit3DCollection* output_digi_hits = m_output_digi_hits.createAndPut();
  for (const auto& input_sim_hit : *input_sim_hits) {
    auto output_digi_hit = output_digi_hits->create();
    // smear the hit position:
    // retrieve the cell detElement
    dd4hep::DDSegmentation::CellID cellID         = input_sim_hit.getCellID();
    debug() << "Digitisation of " << m_readoutName << ", cellID: " << cellID << endmsg;
  //  auto                           cellDetElement = m_volman.lookupDetElement(cellID);
  
   /* const std::string& stripsDetElementName =
        Form("system_%d_layer_%d_phi_%d", m_decoder->get(cellID, "layer"), m_decoder->get(cellID, "phi"));
    dd4hep::DetElement stripsDetElement = cellDetElement.child(stripsDetElementName);
    const auto& stripsTransformMatrix = stripsDetElement.nominal().worldTransformation();
   */

    const auto& stripsTransformMatrix = m_volman.lookupVolumePlacement(cellID).matrix();

    // Retrieve global position in mm and apply unit transformation (translation matrix is tored in cm)
    double simHitGlobalPosition[3] = {input_sim_hit.getPosition().x * dd4hep::mm,
                                      input_sim_hit.getPosition().y * dd4hep::mm,
                                      input_sim_hit.getPosition().z * dd4hep::mm};
    dd4hep::rec::Vector3D simHitGlobalPositionVector(simHitGlobalPosition[0], simHitGlobalPosition[1],
                                                    simHitGlobalPosition[2]);                                  
    double simHitLocalPosition[3]  = {0, 0, 0};

    // get the simHit coordinate in cm in the strips reference frame 
    stripsTransformMatrix.MasterToLocal(simHitGlobalPosition, simHitLocalPosition);
    debug() << "Cell ID string: " << m_decoder->valueString(cellID) << endmsg;
    ;
    debug() << "Global simHit x " << simHitGlobalPosition[0] << " --> Local simHit x " << simHitLocalPosition[0]
            << " in cm" << endmsg;
    debug() << "Global simHit y " << simHitGlobalPosition[1] << " --> Local simHit y " << simHitLocalPosition[1]
            << " in cm" << endmsg;
    debug() << "Global simHit z " << simHitGlobalPosition[2] << " --> Local simHit z " << simHitLocalPosition[2]
            << " in cm" << endmsg;
    // build a vector to easily apply smearing of distance to the strips
    dd4hep::rec::Vector3D simHitLocalPositionVector(simHitLocalPosition[0], simHitLocalPosition[1],
                                                    simHitLocalPosition[2]);
     // smear xy position
    double smearedX = simHitLocalPositionVector.x() ;
    double smearedY = simHitLocalPositionVector.y() + m_gauss_x.shoot() * dd4hep::mm;    
    // smear the z position 
    double smearedZ = simHitLocalPositionVector.z() + m_gauss_y.shoot() * dd4hep::mm;

    double digiHitLocalPosition[3] = {smearedX, smearedY, smearedZ};

    // build the local position vector of the smeared hit. 
 //   dd4hep::rec::Vector3D  digiHitLocalPositionVector(smearedX, smearedY, smearedZ);
    debug() << "Local simHit x: " << simHitLocalPositionVector.x()
            << " Local digiHit x: " << smearedX << " in cm" << endmsg;
    debug() << "Local simHit y: " << simHitLocalPositionVector.y()
            << " Local digiHit y: " << smearedY << " in cm" << endmsg;
    debug() << "Local simHit z: " << simHitLocalPositionVector.z()
            << " Local digiHit z: " << smearedZ << " in cm" << endmsg;
    


    // go back to the global frame
 //   double digiHitLocalPosition[3]  = {digiHitLocalPositionVector.x(), digiHitLocalPositionVector.y(),
  //                                     digiHitLocalPositionVector.z()};
    double digiHitGlobalPosition[3] = {0, 0, 0};
    stripsTransformMatrix.LocalToMaster(digiHitLocalPosition, digiHitGlobalPosition);
    // go back to mm
    edm4hep::Vector3d digiHitGlobalPositionVector(digiHitGlobalPosition[0] / dd4hep::mm,
                                                  digiHitGlobalPosition[1] / dd4hep::mm,
                                                  digiHitGlobalPosition[2] / dd4hep::mm);
    output_digi_hit.setPosition(digiHitGlobalPositionVector);
    output_digi_hit.setCellID(cellID);
  }
  debug() << "Output Digi Hit collection size: " << output_digi_hits->size() << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode MUONsimpleDigitizer::finalize() { return StatusCode::SUCCESS; }
