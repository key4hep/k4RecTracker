#include "VTXdigitizer.h"

DECLARE_COMPONENT(VTXdigitizer)

VTXdigitizer::VTXdigitizer(const std::string& aName, ISvcLocator* aSvcLoc)
    : Gaudi::Algorithm(aName, aSvcLoc), m_geoSvc("GeoSvc", "VTXdigitizer") {
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
  if (m_gauss_time.initialize(m_randSvc, Rndm::Gauss(0., m_t_resolution)).isFailure()) {
    error() << "Couldn't initialize RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  // check if readout exists
  if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }

  // set the cellID decoder
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder(); // Can be used to access e.g. layer index: m_decoder->get(cellID, "layer"),

  // retrieve the volume manager
  m_volman = m_geoSvc->getDetector()->volumeManager();

  return StatusCode::SUCCESS;
}

StatusCode VTXdigitizer::execute(const EventContext&) const {
  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  verbose() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  unsigned nDismissedHits=0;

  // Digitize the sim hits
  edm4hep::TrackerHit3DCollection* output_digi_hits = m_output_digi_hits.createAndPut();
  for (const auto& input_sim_hit : *input_sim_hits) {
    auto output_digi_hit = output_digi_hits->create();

    // smear the hit position: need to go in the local frame of the silicon sensor to smear in the direction along/perpendicular to the stave

    // retrieve the cell detElement
    dd4hep::DDSegmentation::CellID cellID = input_sim_hit.getCellID();
    debug() << "Digitisation of " << m_readoutName << ", cellID: " << cellID << endmsg;

    // Get transformation matrix of sensor
    const auto& sensorTransformMatrix = m_volman.lookupVolumePlacement(cellID).matrix();

    // Retrieve global position in mm and apply unit transformation (translation matrix is stored in cm)
    double simHitGlobalPosition[3] = {input_sim_hit.getPosition().x * dd4hep::mm,
                                      input_sim_hit.getPosition().y * dd4hep::mm,
                                      input_sim_hit.getPosition().z * dd4hep::mm};
    dd4hep::rec::Vector3D simHitGlobalPositionVector(simHitGlobalPosition[0], simHitGlobalPosition[1],
                                                    simHitGlobalPosition[2]);
    double simHitLocalPosition[3]  = {0, 0, 0};

    if(m_forceHitsOntoSurface){
      dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();
      dd4hep::rec::SurfaceManager& surfMan = *theDetector.extension<dd4hep::rec::SurfaceManager>() ;
      dd4hep::DetElement det = m_geoSvc->getDetector()->detector(m_detectorName) ;
      _map = surfMan.map( det.name() ) ;

      if( ! _map )
        error() << " Could not find surface map for detector " << det.name() << " in SurfaceManager " << endmsg;

      dd4hep::rec::SurfaceMap::const_iterator sI = _map->find(cellID) ;
      if( sI == _map->end() )
        error() << " VTXdigitizer: no surface found for cellID " << m_decoder->valueString(cellID) << std::endl << endmsg;

      dd4hep::rec::Vector3D newPos ;
      const dd4hep::rec::ISurface* surf = sI->second ;    

      // Check if Hit is inside sensitive 
      if ( ! surf->insideBounds( simHitGlobalPositionVector ) ) {
        
        info() << "Hit at " << simHitGlobalPositionVector << " is not on surface " << *surf  
                << ". Distance: " << surf->distance(simHitGlobalPositionVector )
                << std::endl << endmsg;        

        if( m_forceHitsOntoSurface ){
          dd4hep::rec::Vector2D lv = surf->globalToLocal(simHitGlobalPositionVector) ;
          dd4hep::rec::Vector3D oldPosOnSurf = surf->localToGlobal( lv ) ; 
          
          info() << "Moved hit to " << oldPosOnSurf << ", distance " << (oldPosOnSurf-simHitGlobalPositionVector).r()
                  << std::endl << endmsg;
            
          simHitGlobalPositionVector = oldPosOnSurf ;
        } 
        else
          ++nDismissedHits;
      }
    }

    // get the simHit coordinate in cm in the sensor reference frame to be able to apply smearing
    sensorTransformMatrix.MasterToLocal(simHitGlobalPosition, simHitLocalPosition);
    debug() << "Cell ID string: " << m_decoder->valueString(cellID) << endmsg;
    ;
    debug() << "Global simHit x " << simHitGlobalPosition[0] << " [mm] --> Local simHit x " << simHitLocalPosition[0]
            << " [cm]" << endmsg;
    debug() << "Global simHit y " << simHitGlobalPosition[1] << " [mm] --> Local simHit y " << simHitLocalPosition[1]
            << " [cm]" << endmsg;
    debug() << "Global simHit z " << simHitGlobalPosition[2] << " [mm] --> Local simHit z " << simHitLocalPosition[2]
            << " [cm]" << endmsg;

    // build a vector to easily apply smearing of distance to the wire
    dd4hep::rec::Vector3D simHitLocalPositionVector(simHitLocalPosition[0], simHitLocalPosition[1],
                                                    simHitLocalPosition[2]);

    // Smear the hit in the local sensor coordinates
    double digiHitLocalPosition[3];
    if (m_readoutName == "VTXIBCollection" ||
        m_readoutName == "VTXOBCollection" ||
        m_readoutName == "VertexBarrelCollection" ||
        m_readoutName == "SiWrBCollection") {  // In barrel, the sensor box is along y-z
      digiHitLocalPosition[0] = simHitLocalPositionVector.x();
      digiHitLocalPosition[1] = simHitLocalPositionVector.y() + m_gauss_x.shoot() * dd4hep::mm;
      digiHitLocalPosition[2] = simHitLocalPositionVector.z() + m_gauss_y.shoot() * dd4hep::mm;
    } else if (m_readoutName == "VTXDCollection" ||
               m_readoutName == "VertexEndcapCollection" ||
               m_readoutName == "SiWrDCollection") {  // In the disks, the sensor box is already in x-y
      digiHitLocalPosition[0] = simHitLocalPositionVector.x() + m_gauss_x.shoot() * dd4hep::mm;
      digiHitLocalPosition[1] = simHitLocalPositionVector.y() + m_gauss_y.shoot() * dd4hep::mm;
      digiHitLocalPosition[2] = simHitLocalPositionVector.z();
    } else {
      error()
          << "VTX readout name (m_readoutName) unknown!"
          << endmsg;
      return StatusCode::FAILURE;
    }

    // go back to the global frame
    double digiHitGlobalPosition[3] = {0, 0, 0};
    sensorTransformMatrix.LocalToMaster(digiHitLocalPosition, digiHitGlobalPosition);

    // go back to mm
    edm4hep::Vector3d digiHitGlobalPositionVector(digiHitGlobalPosition[0] / dd4hep::mm,
                                                  digiHitGlobalPosition[1] / dd4hep::mm,
                                                  digiHitGlobalPosition[2] / dd4hep::mm);

    debug() << "Global digiHit x " << digiHitGlobalPositionVector[0] << " [mm] --> Local digiHit x "
            << digiHitLocalPosition[0] << " [cm]" << endmsg;
    debug() << "Global digiHit y " << digiHitGlobalPositionVector[1] << " [mm] --> Local digiHit y "
            << digiHitLocalPosition[1] << " [cm]" << endmsg;
    debug() << "Global digiHit z " << digiHitGlobalPositionVector[2] << " [mm] --> Local digiHit z "
            << digiHitLocalPosition[2] << " [cm]" << endmsg;
    debug() << "Moving to next hit... " << std::endl << endmsg;

    output_digi_hit.setEDep(input_sim_hit.getEDep());
    output_digi_hit.setPosition(digiHitGlobalPositionVector);

    // Apply time smearing
    output_digi_hit.setTime(input_sim_hit.getTime() + m_gauss_time.shoot());

    output_digi_hit.setCellID(cellID);
  }
  return StatusCode::SUCCESS;
}

StatusCode VTXdigitizer::finalize() { return StatusCode::SUCCESS; }
