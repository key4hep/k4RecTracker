#include "VTXdigitizer.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"

DECLARE_COMPONENT(VTXdigitizer)

VTXdigitizer::VTXdigitizer(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc), m_geoSvc("GeoSvc", "VTXdigitizer")  {
  declareProperty("inputSimHits", m_input_sim_hits, "Input sim vertex hit collection name");
  declareProperty("outputDigiHits", m_output_digi_hits, "Output digitized vertex hit collection name");
}

VTXdigitizer::~VTXdigitizer() {}

StatusCode VTXdigitizer::initialize() { 
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
  if (m_gauss_t.initialize(m_randSvc, Rndm::Gauss(0., m_t_resolution)).isFailure()) {
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

  return StatusCode::SUCCESS; }

StatusCode VTXdigitizer::execute() {
  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  verbose() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Digitize the sim hits
  edm4hep::TrackerHitCollection* output_digi_hits = m_output_digi_hits.createAndPut();
  for (const auto& input_sim_hit : *input_sim_hits) {
    auto output_digi_hit = output_digi_hits->create();

    warning() << "Test!!!" << endmsg;




    // smear the hit position: need to go in the local frame of the silicon sensor to smear in the direction along/perpendicular to the stave

    // retrieve the cell detElement
    dd4hep::DDSegmentation::CellID cellID         = input_sim_hit.getCellID();
    auto                           cellDetElement = m_volman.lookupDetElement(cellID);
    warning() << "cellID: " << cellID << endmsg;

    // retrieve the sensor          ###### (in DD4hep 1.23 there is no easy way to access the volume daughters we have to pass by detElements, in later versions volumes can be used)

    // // VTXIB
    // const std::string& sensorDetElementName =
    //     Form("VTXIB_layer%d_stave_sensors_module%d_sensorMotherVolume%d", m_decoder->get(cellID, "layer"),
    //           m_decoder->get(cellID, "module"), m_decoder->get(cellID, "sensor"));
    // warning() << "sensorDetElementName: " << sensorDetElementName << endmsg;

    // VTXD
    const std::string& sensorDetElementName =
        Form("VTXD_side%d_layer%d_petal%d_stave%d_sensors_module%d_sensor%d", m_decoder->get(cellID, "side"),
              m_decoder->get(cellID, "layer"), m_decoder->get(cellID, "petal"), m_decoder->get(cellID, "stave"), m_decoder->get(cellID, "module"), m_decoder->get(cellID, "sensor"));
    warning() << "sensorDetElementName: " << sensorDetElementName << endmsg;

    // // CLD Barrel
    // const std::string& sensorDetElementName =
    //     Form("VTXIB_layer%d_stave_sensors_module%d_sensorMotherVolume%d", m_decoder->get(cellID, "layer"),
    //           m_decoder->get(cellID, "module"), m_decoder->get(cellID, "sensor"));
    // warning() << "sensorDetElementName: " << sensorDetElementName << endmsg;    

    dd4hep::DetElement sensorDetElement = cellDetElement; // cellDetElement.child(sensorDetElementName);

    // get the transformation matrix used to place the wire
    const auto& sensorTransformMatrix = sensorDetElement.nominal().worldTransformation();

    // Retrieve global position in mm and apply unit transformation (translation matrix is tored in cm)
    double simHitGlobalPosition[3] = {input_sim_hit.getPosition().x * dd4hep::mm,
                                      input_sim_hit.getPosition().y * dd4hep::mm,
                                      input_sim_hit.getPosition().z * dd4hep::mm}; // Unit correct for vertex?
    double simHitLocalPosition[3]  = {0, 0, 0};

    // get the simHit coordinate in cm in the sensor reference frame to be able to apply smearing
    sensorTransformMatrix.MasterToLocal(simHitGlobalPosition, simHitLocalPosition);
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
    double smearedX = simHitLocalPositionVector.x() + m_gauss_x.shoot() * dd4hep::mm;
    double smearedY = simHitLocalPositionVector.y() + m_gauss_y.shoot() * dd4hep::mm;

    // go back to the global frame
    double digiHitLocalPosition[3]  = {smearedX, smearedY,
                                        simHitLocalPositionVector.z()};
    double digiHitGlobalPosition[3] = {0, 0, 0};
    sensorTransformMatrix.LocalToMaster(digiHitLocalPosition, digiHitGlobalPosition);
    
    // go back to mm
    edm4hep::Vector3d digiHitGlobalPositionVector(digiHitGlobalPosition[0] / dd4hep::mm,
                                                  digiHitGlobalPosition[1] / dd4hep::mm,
                                                  digiHitGlobalPosition[2] / dd4hep::mm);



    output_digi_hit.setEDep(input_sim_hit.getEDep());
    output_digi_hit.setPosition(digiHitGlobalPositionVector);
  }
  return StatusCode::SUCCESS;
}

StatusCode VTXdigitizer::finalize() { return StatusCode::SUCCESS; }
