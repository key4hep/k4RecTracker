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
  // Initialise LocalNormalVectorDir with forcing the user to declare a value
  if ( !(m_LocalNormalVectorDir=="x" || m_LocalNormalVectorDir=="y" || m_LocalNormalVectorDir=="z") ) {
    error() << "LocalNormalVectorDir property should be declared as a string with the direction (x,y or z) of the normal vector to sensitive surface in the sensor local frame (may differ according to the geometry definition within k4geo). Add a - sign before the direction in case of indirect frame." << endmsg;
    return StatusCode::FAILURE;
  }
  
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

  // Get the sensor thickness per layer in mm
  //getSensorThickness();
  dd4hep::DetType type(m_geoSvc->getDetector()->detector(m_detectorName).typeFlag()); // Get detector Type
  if (type.is(dd4hep::DetType::BARREL)) // if this is a barrel detector
    getSensorThickness<dd4hep::rec::ZPlanarData>();
  else if (type.is(dd4hep::DetType::ENDCAP)) // If this is an Endcap Detector
    getSensorThickness<dd4hep::rec::ZDiskPetalsData>();
  else {
    error() << m_detectorName << " : Detector type should be BARREL or ENDCAP " << endmsg;
    return StatusCode::FAILURE;
  }

  // Initilize the diffusion parameters
  m_Dist50 = 0.050; // Define 50microns in mm
  
  return StatusCode::SUCCESS;
}

template<typename T> void VTXdigitizerDetailed::getSensorThickness() {
  /** Retrieve the sensor thickness per layer in millimeter
   */
  
  auto detElement = m_geoSvc->getDetector()->detector(m_detectorName);
  auto* theExtension = detElement.extension<T>();
  const std::vector<typename T::LayerLayout>& layers = theExtension->layers;
  m_sensorThickness.resize(layers.size());
  int layerCount = 0;
  for (auto const& layer : layers) {
    m_sensorThickness[layerCount] = (float)(layer.thicknessSensitive/dd4hep::mm);
    ++layerCount;
  }
  
} // End getSensorThickness

StatusCode VTXdigitizerDetailed::execute(const EventContext&) const {
  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  verbose() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Digitize the sim hits
  edm4hep::TrackerHitPlaneCollection* output_digi_hits = m_output_digi_hits.createAndPut();
  edm4hep::TrackerHitSimTrackerHitLinkCollection* output_sim_digi_link_col = m_output_sim_digi_link.createAndPut();
  
  for (const auto& input_sim_hit : *input_sim_hits) {
    
    // TEST new digitizer members
    
    std::vector<ChargeDepositUnit> ionizationPoints;
    std::vector<SignalPoint> collectionPoints;
    primary_ionization(input_sim_hit, ionizationPoints);
    drift(input_sim_hit, ionizationPoints, collectionPoints);
    
    // END TEST new digitizer members
      
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

  debug() << "Enter Primary_ionization, numberOfSegments=" << numberOfSegments << ", shift=" << pathLength << ", energy=" << eDep << "GeV" << endmsg;
  
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
  const dd4hep::DDSegmentation::CellID& cellID = hit.getCellID();

  // Get transformation matrix of sensor and get central position and momentum vector in local frame
  //const auto& sensorTransformMatrix = m_volman.lookupVolumePlacement(cellID).matrix();
  TGeoHMatrix sensorTransformMatrix = m_volman.lookupDetElement(cellID).nominal().worldTransformation();
  SetProperDirectFrame(sensorTransformMatrix); // Change coordinates to have z orthogonal to sensor with direct frame
  double hitLocalCentralPosition[3] = {0, 0, 0};
  double hitLocalMomentum[3] = {0, 0, 0};

  // get the simHit coordinate in cm in the sensor reference frame
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

void VTXdigitizerDetailed::drift(const edm4hep::SimTrackerHit& hit, const std::vector<ChargeDepositUnit>& ionizationPoints, std::vector<SignalPoint>& collectionPoints) const{
  /** Drift the charge segments to the sensor surface (collection plane)
   *  Include the effect of E-field and B-field
   *  Even for strips, consider a uniform E-field
   *  The sensor depth is explicitly considered to be along z
   *  For now, the drift is considered to be instantaneous. Later On, should consider time and time uncertainty
   */
  
  debug() << "Enter drift " << endmsg;

  collectionPoints.resize(ionizationPoints.size()); //set size

  dd4hep::rec::Vector3D driftDir = DriftDirection(hit); // Get the charge (electrons) drift direction
  if (driftDir.z() == 0.) {
    warning() << "Charge drift in z is 0" << endmsg;
    return;
  }

  float tanLorentzAngleX = driftDir.x(); // tan of Lorentz angle
  float tanLorentzAngleY = driftDir.y();
  float dir_z = driftDir.z(); // The z drift direction
  float cosLorentzAngleX = 1. / sqrt(1. + tanLorentzAngleX * tanLorentzAngleX); // Cosine
  float cosLorentzAngleY = 1. / sqrt(1. + tanLorentzAngleY * tanLorentzAngleY);

  // Get the sensor thickness
  const dd4hep::DDSegmentation::CellID& cellID = hit.getCellID();
  int layerID = m_decoder->get(cellID, "layer");
  double moduleThickness = m_sensorThickness[layerID];
  
  float SigmaX = 1.; // Charge spread
  float SigmaY = 1.;
  float DriftDistance; // Distance between charge generation and collection
  float DriftLength;   // Actual drift length
  float Sigma;

  for (unsigned int i = 0; i < ionizationPoints.size(); i++) {
    float SegX, SegY, SegZ; // Position
    SegX = ionizationPoints[i].x();
    SegY = ionizationPoints[i].y();
    SegZ = ionizationPoints[i].z();

    // Distance from the collection plane
    // Include explicitly the E drift direction (to positive z)
    DriftDistance = moduleThickness / 2. - (dir_z * SegZ);

    if (DriftDistance < 0.)
      DriftDistance = 0.;
    else if (DriftDistance > moduleThickness)
      DriftDistance = moduleThickness;

    // Assume full depletion
    float XDriftDueToMagField = DriftDistance * tanLorentzAngleX;
    float YDriftDueToMagField = DriftDistance * tanLorentzAngleY;

    // Shift cloud center
    float CloudCenterX = SegX + XDriftDueToMagField;
    float CloudCenterY = SegY + YDriftDueToMagField;

    // Calculate how long is the charge drift path
    DriftLength = sqrt(DriftDistance * DriftDistance
		       + XDriftDueToMagField * XDriftDueToMagField
		       + YDriftDueToMagField * YDriftDueToMagField);

    // What is the charge diffusion after this path
    // One considers a reference value for the diffusion on 50mum drift
    Sigma = sqrt(DriftLength / m_Dist50) * m_Sigma50;

    // Project the diffusion sigma on the collection plane
    SigmaX = Sigma / cosLorentzAngleX;
    SigmaY = Sigma / cosLorentzAngleY;

    SignalPoint sp(CloudCenterX, CloudCenterY, SigmaX, SigmaY, hit.getTime(), hit, ionizationPoints[i].charge());

    // Load the charge distribution parameters
    collectionPoints[i] = sp;
      
  } // End loop over ionization points
  
} // End drift

dd4hep::rec::Vector3D VTXdigitizerDetailed::DriftDirection(const edm4hep::SimTrackerHit& hit) const {
  /** Set the drif direction according to the B-field in local sensor frame
   *  Works for both barrel and endcap pixels
   *  Replace de sign convention to fit M.Swartz's formulaes
   *  Possibility to set the Lorentz angle through the option file (default=0.1)
   */

  // Retrieve the B-field at hit position in local frame
  // Get detector
  const dd4hep::Detector* detector = m_geoSvc->getDetector();

  // Get hit global position
  dd4hep::rec::Vector3D hitPosition(hit.getPosition().x * dd4hep::mm,
				    hit.getPosition().y * dd4hep::mm,
				    hit.getPosition().z * dd4hep::mm);

  // Get B field in global frame
  const auto& magneticField = detector->field();
  double BFieldGlobal[3];
  magneticField.magneticField(hitPosition,BFieldGlobal);
    
  // Get B field in local frame
  double BFieldLocal[3] = { 0, 0, 0 };
  const dd4hep::DDSegmentation::CellID& cellID = hit.getCellID();
  auto sensorTransformMatrix = m_volman.lookupDetElement(cellID).nominal().worldTransformation();
  SetProperDirectFrame(sensorTransformMatrix); // Change coordinates to have z orthogonal to sensor with direct frame
  sensorTransformMatrix.MasterToLocalVect(BFieldGlobal, BFieldLocal);

  dd4hep::rec::Vector3D BFieldLocalVec(BFieldLocal[0], BFieldLocal[1], BFieldLocal[2]);

  // Get the drift direction
  float dir_x = 0.;
  float dir_y = 0.;
  float dir_z = 0.;
  float scale = 0.;
  float alpha2 = m_tanLorentzAnglePerTesla * m_tanLorentzAnglePerTesla;

  // get the directions considering E-field along z and B field (Lorentz angle)
  // Consider here that electrons go to +z (E field in -z direction). Signs to be change if reverse, accordingly to CMS digitizer.
  dir_x = -m_tanLorentzAnglePerTesla * BFieldLocalVec.y() + alpha2 * BFieldLocalVec.z() * BFieldLocalVec.x();
  dir_y =  m_tanLorentzAnglePerTesla * BFieldLocalVec.x() + alpha2 * BFieldLocalVec.z() * BFieldLocalVec.y();
  dir_z =  1 + alpha2 * BFieldLocalVec.z() * BFieldLocalVec.z();
  scale = dir_z; // Normalize to z direction
  
  dd4hep::rec::Vector3D theDriftDirection(dir_x / scale, dir_y / scale, dir_z /scale);
  
  return theDriftDirection;
  
} // End DriftDirection

void VTXdigitizerDetailed::SetProperDirectFrame(TGeoHMatrix& sensorTransformMatrix) const {
  /** Change the sensorTransformMatrix to have a direct frame with z orthogonal to sensor surface
   */
  
  std::string LocalNormalVectorDir = m_LocalNormalVectorDir;
  bool IsDirect = true; // Is the origin frame direct ?
  if (LocalNormalVectorDir[0]=='-') {
    IsDirect = false;
    LocalNormalVectorDir = LocalNormalVectorDir[1];
  }  

  // If the orthogonal direction is along X or Y in local frame, rotate the frame to have Z orthogonal to sensors instead
  if (LocalNormalVectorDir=="x") {
    TGeoRotation rot("rot",90.,90.,0.);
    sensorTransformMatrix.Multiply(rot);
  }
  if (LocalNormalVectorDir=="y") {
    TGeoRotation rot("rot",0.,-90.,0.);
    sensorTransformMatrix.Multiply(rot);
  }

  // If the frame isn't direct, make it direct by reflecting the x axis. This is necessary to correctly calculte the drift in X-Y due to B-field
  if (!IsDirect) {
    sensorTransformMatrix.ReflectX(false);
  }
} // End SetProperDirectFrame
