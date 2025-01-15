#include "VTXdigitizerDetailed.h"

DECLARE_COMPONENT(VTXdigitizerDetailed)

VTXdigitizerDetailed::VTXdigitizerDetailed(const std::string& aName, ISvcLocator* aSvcLoc)
    : Gaudi::Algorithm(aName, aSvcLoc), m_geoSvc("GeoSvc", "VTXdigitizerDetailed") {
  declareProperty("inputSimHits", m_input_sim_hits, "Input sim vertex hit collection name");
  declareProperty("outputDigiHits", m_output_digi_hits, "Output digitized vertex hit collection name");
  declareProperty("outputSimDigiAssociation", m_output_sim_digi_link, "Output link between sim hits and digitized hits");
}

VTXdigitizerDetailed::~VTXdigitizerDetailed() {}

StatusCode VTXdigitizerDetailed::initialize() {
  // Initialize random services
  m_randSvc = service("RndmGenSvc", false);
  if (!m_randSvc) {
    error() << "Couldn't get RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  m_gauss_x_vec.resize(m_x_resolution.size());
  for (size_t i = 0; i < m_x_resolution.size(); ++i) {
    if (m_gauss_x_vec[i].initialize(m_randSvc, Rndm::Gauss(0., m_x_resolution[i])).isFailure()) {
      error() << "Couldn't initialize RndmGenSvc!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  m_gauss_y_vec.resize(m_y_resolution.size());
  for (size_t i = 0; i < m_y_resolution.size(); ++i) {  
    if (m_gauss_y_vec[i].initialize(m_randSvc, Rndm::Gauss(0., m_y_resolution[i])).isFailure()) {
      error() << "Couldn't initialize RndmGenSvc!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  m_gauss_t_vec.resize(m_t_resolution.size());
  for (size_t i = 0; i < m_t_resolution.size(); ++i) {  
    if (m_gauss_t_vec[i].initialize(m_randSvc, Rndm::Gauss(0., m_t_resolution[i])).isFailure()) {
      error() << "Couldn't initialize RndmGenSvc!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // check if readout exists
  if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }

  // set the cellID decoder
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder(); // Can be used to access e.g. layer index: m_decoder->get(cellID, "layer"),

  if (m_decoder->fieldDescription().find("layer") == std::string::npos){
    error() 
      << " Readout " << m_readoutName << " does not contain layer id!"
      << endmsg;
    return StatusCode::FAILURE;
  }
  
  // retrieve the volume manager
  m_volman = m_geoSvc->getDetector()->volumeManager();

  return StatusCode::SUCCESS;
}

StatusCode VTXdigitizerDetailed::execute(const EventContext&) const {
  // Get the input collection with Geant4 hits
    const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  verbose() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Digitize the sim hits
  edm4hep::TrackerHitPlaneCollection* output_digi_hits = m_output_digi_hits.createAndPut();
  edm4hep::TrackerHitSimTrackerHitLinkCollection* output_sim_digi_link_col = m_output_sim_digi_link.createAndPut();
  
  for (const auto& input_sim_hit : *input_sim_hits) {
    auto output_digi_hit = output_digi_hits->create();
    auto output_sim_digi_link = output_sim_digi_link_col->create();

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
        error() << " VTXdigitizerDetailed: no surface found for cellID " << m_decoder->valueString(cellID) << std::endl << endmsg;

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
      }
    }
    
    // get the simHit coordinate in cm in the sensor reference frame to be able to apply smearing
    sensorTransformMatrix.MasterToLocal(simHitGlobalPosition, simHitLocalPosition);
    debug() << "Cell ID string: " << m_decoder->valueString(cellID) << endmsg;
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
    int iLayer = m_decoder->get(cellID, "layer");      
    debug() << "readout: " << m_readoutName << ", layer id: " << iLayer << endmsg;
    if (m_readoutName == "VertexBarrelCollection" ||
        m_readoutName == "SiWrBCollection" ||
	m_readoutName == "InnerTrackerBarrelCollection" ||
	m_readoutName == "OuterTrackerBarrelCollection") {  // In barrel, the sensor box is along y-z
      digiHitLocalPosition[0] = simHitLocalPositionVector.x();
      digiHitLocalPosition[1] = simHitLocalPositionVector.y() + m_gauss_x_vec[iLayer].shoot() * dd4hep::mm;
      digiHitLocalPosition[2] = simHitLocalPositionVector.z() + m_gauss_y_vec[iLayer].shoot() * dd4hep::mm;
    } else if (m_readoutName == "VertexEndcapCollection" ||
               m_readoutName == "SiWrDCollection" ||
	       m_readoutName == "InnerTrackerEndcapCollection" ||
	       m_readoutName == "OuterTrackerEndcapCollection") {  // In the disks, the sensor box is already in x-y
      digiHitLocalPosition[0] = simHitLocalPositionVector.x() + m_gauss_x_vec[iLayer].shoot() * dd4hep::mm;
      digiHitLocalPosition[1] = simHitLocalPositionVector.y() + m_gauss_y_vec[iLayer].shoot() * dd4hep::mm;
      digiHitLocalPosition[2] = simHitLocalPositionVector.z();
    } else {
      error()
          << "VTX readout name (m_readoutName) unknown or xResolution/yResolution/tResolution not defining all detector layer resolutions!"
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
    output_digi_hit.setTime(input_sim_hit.getTime() + m_gauss_t_vec[iLayer].shoot());

    output_digi_hit.setCellID(cellID);

    // Set the link between sim and digi hit
    output_sim_digi_link.setFrom(output_digi_hit);
    output_sim_digi_link.setTo(input_sim_hit);
    
  }
  return StatusCode::SUCCESS;
}

StatusCode VTXdigitizerDetailed::finalize() { return StatusCode::SUCCESS; }

void VTXdigitizerDetailed::primary_ionization(const edm4hep::SimTrackerHit& hit, std::vector<ChargeDepositUnit>& ionizationPoints) const {
  /** Generate primary ionization along the track segment.
   *Divide the track into small sub-segments of 5 microns
   *Straight line approximation for trajectory of the incoming particle inside the active media
   */
  
  const float segmentLength = 0.005 * dd4hep::mm; //5microns in mm
  float pathLength = hit.getPathLength() * dd4hep::mm; // Path Length of the particle in the active material
  int numberOfSegments = int(pathLength / segmentLength); // Number of segments
  if (numberOfSegments <1) { numberOfSegments = 1; }

  float GeVperElectron = 3.61E-09 * dd4hep::GeV; // Mean ionization energy in silicon
  float charge;
  float eDep = hit.getEDep() * dd4hep::GeV; // Total energy deposited by the particle in this hit

  // Get the global position of the hit (defined by default in Geant4 as the mean between the entry and exit point in the active material)
  // and apply unit transformation (translation matrix is stored in cm
  double hitGlobalCentralPosition[3] = {hit.getPosition().x * dd4hep::mm,
                                        hit.getPosition().y * dd4hep::mm,
					hit.getPosition().z * dd4hep::mm}; 
  
  // Get the 3D momentum in global frame
  double hitGlobalMomentum[3] = {hit.getMomentum().x * dd4hep::GeV,
                                 hit.getMomentum().y * dd4hep::GeV,
				 hit.getMomentum().z * dd4hep::GeV};

  // Convert the global position and momentum into local frame
  // retrieve the cell detElement
  dd4hep::DDSegmentation::CellID cellID = hit.getCellID();
  
  // Get transformation matrix of sensor and get central position and momentum vector in local frame
  const auto& sensorTransformMatrix = m_volman.lookupVolumePlacement(cellID).matrix();
  double hitLocalCentralPosition[3] = {0, 0, 0};
  double hitLocalMomentum[3] = {0, 0, 0};
  
  sensorTransformMatrix.MasterToLocal(hitGlobalCentralPosition, hitLocalCentralPosition);
  sensorTransformMatrix.MasterToLocalVect(hitGlobalMomentum, hitLocalMomentum);

  // build vectors to simplify next calculation
  dd4hep::rec::Vector3D hitLocalCentralPositionVec(hitLocalCentralPosition[0],
						   hitLocalCentralPosition[1],
						   hitLocalCentralPosition[2]);
  dd4hep::rec::Vector3D hitLocalMomentumVec(hitLocalMomentum[0],
					    hitLocalMomentum[1],
					    hitLocalMomentum[2]);

  // Get particle direction in local frame normalized to path Length
  dd4hep::rec::Vector3D direction = (pathLength/hitLocalMomentumVec.r())*hitLocalMomentumVec; 

  // Get the entry point in the active material
  dd4hep::rec::Vector3D entryPoint = hitLocalCentralPositionVec - 0.5 * direction;

  // Maybe add a test that the point is inside the material or go to the closest material and do the same for exit point and then redefine the length
  
  // Calulate the charge deposited per segment
  ionizationPoints.resize(numberOfSegments);
  for (int i = 0; i != numberOfSegments; i++) {
    // Divide the segment into equal length sub-segments
    dd4hep::rec::Vector3D point = entryPoint + float((i + 0.5) / numberOfSegments) * direction;
    
    charge = eDep / GeVperElectron / float(numberOfSegments); // Here, do not consider energy fluctuation. To consider it, see the Geant4 Physics Reference Manual

    ChargeDepositUnit cdu(charge, point); // Define position,charge point
    ionizationPoints[i] = cdu; // save
    
  }// end loop over segments
} // End primary_ionization
  
