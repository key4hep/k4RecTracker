#include "VTXdigitizerDetailed.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Readout.h"

DECLARE_COMPONENT(VTXdigitizerDetailed)

VTXdigitizerDetailed::VTXdigitizerDetailed(const std::string& aname, ISvcLocator* asvcLoc)
    : MultiTransformer(aname, asvcLoc, {KeyValue("inputSimHits", "")},
                       {KeyValue("outputDigiHits", ""), KeyValue("outputSimDigiAssociation", "")}) {}

StatusCode VTXdigitizerDetailed::initialize() {

  // Check that input and output collections names are declared by the user

  if (inputLocations("inputSimHits")[0].empty()) {
    error() << "You must specify the inputSimHits collection name!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (outputLocations("outputDigiHits")[0].empty()) {
    error() << "You must specify the outputDigiHits collection name!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (outputLocations("outputSimDigiAssociation")[0].empty()) {
    error() << "You must specify the outputSimDigiAssociation collection name!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Initialise the debugging histograms
  if (m_DebugHistos) {
    if (m_DebugFileName == "") {
      error() << "Please provide a name for the file containing the debug histograms with option DebugFileName."
              << endmsg;
      return StatusCode::FAILURE;
    }

    hErrorX = new TH1D("hErrorX", "Distance in X in local frame between the true hit and digitized one in mm", 10000,
                       -0.05, 0.05);
    hErrorX->SetDirectory(0);
    hErrorY = new TH1D("hErrorY", "Distance in Y in local frame between the true hit and digitized one in mm", 10000,
                       -0.05, 0.05);
    hErrorY->SetDirectory(0);
    hErrorZ = new TH1D("hErrorZ", "Distance in Z in local frame between the true hit and digitized one in mm", 10000,
                       -0.05, 0.05);
    hErrorZ->SetDirectory(0);
    hError = new TH1D("hError", "Distance between the true hit and digitized one in mm", 10000, 0., 0.05);
    hError->SetDirectory(0);

    // Before Threshold

    hChargeBeforeThreshold = new TH1D("hChargeBeforeThreshold", "Pixel Charge before threshold", 4000, 0., 20.);
    hChargeBeforeThreshold->SetXTitle("Pixel Charge Without Threshold (ke)");
    hChargeBeforeThreshold->SetYTitle("Pixels");
    hChargeBeforeThreshold->SetDirectory(0);

    hActivePixelCountBeforeThreshold =
        new TH1D("hActivePixelCountPerCluster", "Active Pixels Before Threshold per Cluster", 60, 0, 30);
    hActivePixelCountBeforeThreshold->SetXTitle("Active Pixels");
    hActivePixelCountBeforeThreshold->SetYTitle("Clusters");
    hActivePixelCountBeforeThreshold->SetDirectory(0);

    // After Threshold

    hChargeAboveThreshold = new TH1D("hChargeAboveThreshold", "Pixel Charge after threshold", 4000, 0., 20.);
    hChargeAboveThreshold->SetXTitle("Pixel Charge (ke)");
    hChargeAboveThreshold->SetYTitle("Pixels");
    hChargeAboveThreshold->SetDirectory(0);

    hActivePixelCountAfterThreshold =
        new TH1D("hActivePixelCountAfterThresholdPerCluster", "Active Pixels After Threshold per Cluster", 60, 0, 30);
    hActivePixelCountAfterThreshold->SetXTitle("Active Pixels");
    hActivePixelCountAfterThreshold->SetYTitle("Clusters");
    hActivePixelCountAfterThreshold->SetDirectory(0);

    // Charge per cluster of pixels

    hChargePerCluster = new TH1D("hChargePerCluster", "Charge per Cluster", 3000, 0., 30.);
    hChargePerCluster->SetXTitle("Charge in Cluster (ke)");
    hChargePerCluster->SetYTitle("Clusters");
    hChargePerCluster->SetDirectory(0);

    /// Drift due to magnetic field

    hXDriftDueToMagField = new TH1D("hXDriftDueToMagField", "X Drift due to magnetic field", 2000, -0.01, 0.01);
    hXDriftDueToMagField->SetXTitle("X Drift due to magnetic field (mm)");
    hXDriftDueToMagField->SetDirectory(0);

    hYDriftDueToMagField = new TH1D("hYDriftDueToMagField", "Y Drift due to magnetic field", 2000, -0.01, 0.01);
    hYDriftDueToMagField->SetXTitle("Y Drift due to magnetic field (mm)");
    hYDriftDueToMagField->SetDirectory(0);

    // Cluster and fired pixels per layer

    hClusterPerLayer = new TH1D("hClusterPerLayer", "Cluster per Layer", 10, 0, 10);
    hClusterPerLayer->SetDirectory(0);

    hActivePixelPerlayer = new TH1D("hActivePixelPerlayer", "Active Pixels per Layer", 10, 0, 10);
    hActivePixelPerlayer->SetXTitle("Layer");
    hActivePixelPerlayer->SetYTitle("Number of active pixels");
    hActivePixelPerlayer->SetDirectory(0);
  }

  // Initialize random services
  m_randSvc = service("RndmGenSvc", false);
  if (!m_randSvc) {
    error() << "Couldn't get RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Initialize Gaussian random generator for time smearing
  m_gauss_t_vec.resize(m_t_resolution.size());
  for (size_t i = 0; i < m_t_resolution.size(); ++i) {
    if (m_gauss_t_vec[i].initialize(m_randSvc, Rndm::Gauss(0., m_t_resolution[i])).isFailure()) {
      error() << "Couldn't initialize RndmGenSvc!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // Initialize Gaussian random generator for threshold smearing
  if (m_gauss_threshold.initialize(m_randSvc, Rndm::Gauss(m_Threshold, m_ThresholdSmearing)).isFailure()) {
    error() << "Couldn't initialize Gaussian generator for threshold smearing!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Initialize and check if GeoSvc is valid
  m_geoSvc = serviceLocator()->service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Failed to retrieve GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  info() << "GeoSvc successfully retrieved" << endmsg;

  // Get detector from geoSvc
  const auto detector = m_geoSvc->getDetector();

  // Set the cellID decoder

  const std::string metadataKey =
      podio::collMetadataParamName(inputLocations("inputSimHits")[0], edm4hep::labels::CellIDEncoding);
  std::string encodingString;
  try {
    auto optionalEncoding = k4FWCore::getParameter<std::string>(metadataKey, this);
    if (!optionalEncoding.has_value()) {
      error() << "CellIDEncoding metadata key '" << metadataKey << "' not found for collection '"
              << inputLocations("inputSimHits")[0] << "'. Cannot initialize decoder." << endmsg;
      return StatusCode::FAILURE;
    }
    encodingString = optionalEncoding.value();
  } catch (const std::exception& e) {
    error() << "An unexpected exception occurred while retrieving metadata: " << e.what() << endmsg;
    return StatusCode::FAILURE;
  }

  try {
    m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(encodingString);
  } catch (const std::runtime_error& e) {
    error() << "Failed to create BitFieldCoder for collection " << inputLocations("inputSimHits")[0]
            << " using encoding: " << encodingString << ": " << e.what() << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_decoder->fieldDescription().find("layer") == std::string::npos) {
    error() << " Collection " << inputLocations("inputSimHits")[0] << " does not contain layer id!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Retrieve the surface map from the SurfaceManager

  const auto surfMan = detector->extension<dd4hep::rec::SurfaceManager>();
  surfaceMap = surfMan->map(m_detectorName);

  if (!surfaceMap) {
    throw std::runtime_error(
        fmt::format("Could not find surface map for detector: {} in SurfaceManager", m_detectorName.value()));
  }

  // retrieve the volume manager
  m_volman = detector->volumeManager();

  // Get the sensor thickness and 2D size per layer in mm
  dd4hep::DetType type(detector->detector(m_detectorName).typeFlag()); // Get detector Type
  if (type.is(dd4hep::DetType::BARREL)) {                              // if this is a barrel detector
    getSensorThickness<dd4hep::rec::ZPlanarData>();
  } else if (type.is(dd4hep::DetType::ENDCAP)) { // If this is an Endcap Detector
    getSensorThickness<dd4hep::rec::ZDiskPetalsData>();
  } else {
    error() << m_detectorName << " : Detector type should be BARREL or ENDCAP " << endmsg;
    return StatusCode::FAILURE;
  }

  debug() << "Initializing VTXdigitizerDetailed with the following parameters:" << endmsg;
  debug() << "Detector Name: " << m_detectorName << endmsg;

  debug() << "Threshold: " << m_Threshold << endmsg;
  debug() << "Threshold Smearing : " << m_ThresholdSmearing << endmsg;

  return StatusCode::SUCCESS;
}

template <typename T>
void VTXdigitizerDetailed::getSensorThickness() {
  /** Retrieve the sensor thickness per layer in millimeter
   */

  auto detElement = m_geoSvc->getDetector()->detector(m_detectorName);
  auto* theExtension = detElement.extension<T>();
  const std::vector<typename T::LayerLayout>& layers = theExtension->layers;
  m_sensorThickness.resize(layers.size());
  int layerCount = 0;
  for (auto const& layer : layers) {
    m_sensorThickness[layerCount] = (float)(layer.thicknessSensitive / dd4hep::mm);
    ++layerCount;
  }
  debug() << "Sensor thickness per layer (mm): " << m_sensorThickness[0] << endmsg;
} // End getSensorThickness

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection>
VTXdigitizerDetailed::operator()(const edm4hep::SimTrackerHitCollection& inputSimHits) const {
  // Get the input collection with Geant4 hits
  verbose() << "Input Sim Hit collection size: " << inputSimHits.size() << endmsg;

  // Check if it is produced by Secondaries or Overlay
  int nSecondaryHits = 0;
  int nOverlayHits = 0;

  // Digitize the sim hits
  auto output_digi_hits = edm4hep::TrackerHitPlaneCollection();
  auto output_sim_digi_link_col = edm4hep::TrackerHitSimTrackerHitLinkCollection();

  for (const auto& input_sim_hit : inputSimHits) {

    if (input_sim_hit.getTime() > m_maxHitTime) {
      debug() << "Skipping SimHit with large time: " << input_sim_hit.getTime()
              << " ns (CellID: " << input_sim_hit.getCellID() << ")" << endmsg;
      continue;
    }

    // Get the normal vector direction in local frame during the first hit digitization
    if (m_LocalNormalVectorDir == "") {
      debug() << "First hit digitization, retrieve the local normal vector direction" << endmsg;
      GetNormalVectorLocal(input_sim_hit);
    } else {
      debug() << "Local normal vector direction already retrieved: " << m_LocalNormalVectorDir << endmsg;
    }

    // Check if the hit is a secondary or an overlay
    if (input_sim_hit.isProducedBySecondary()) {
      nSecondaryHits++;
    } else if (input_sim_hit.isOverlay()) {
      nOverlayHits++;
    }

    // Digitization process
    auto ionizationPoints = primary_ionization(input_sim_hit);

    auto collectionPoints = drift(input_sim_hit, ionizationPoints);

    auto hit_map = get_charge_per_pixel(input_sim_hit, collectionPoints);

    generate_output(input_sim_hit, output_digi_hits, output_sim_digi_link_col, hit_map);
  }

  debug() << "SimHits stats: " << nSecondaryHits << " produced by secondaries, " << nOverlayHits
          << " marked as overlay." << endmsg;

  // Cluster count per layer performed only when debug mode is enabled
  if (m_DebugHistos) {
    std::map<int, int> digisPerLayer;

    for (const auto& digi_hit : output_digi_hits) {
      uint64_t cellID = digi_hit.getCellID();
      int layer = m_decoder->get(cellID, "layer");
      // int side = m_decoder->get(cellID, "side");
      hClusterPerLayer->Fill(layer);
    }
  }
  info() << "Execution of VTXdigitizerDetailed completed." << endmsg;

  return std::make_tuple(std::move(output_digi_hits), std::move(output_sim_digi_link_col));
}

StatusCode VTXdigitizerDetailed::finalize() {

  if (m_DebugHistos) {
    Create_outputROOTfile_for_debugHistograms();
  }

  return StatusCode::SUCCESS;
}

void VTXdigitizerDetailed::GetNormalVectorLocal(const edm4hep::SimTrackerHit& input_sim_hit) const {
  /** Retrieve the normal vector direction in local frame (normal to sensor surface)
   *  This function is called during the first hit digitization to initialise the member m_LocalNormalVectorDir
   */

  const dd4hep::DDSegmentation::CellID& cellID = input_sim_hit.getCellID();
  std::uint64_t m_mask = (static_cast<std::uint64_t>(1) << 32) - 1;
  const std::uint64_t reduced_cellID = cellID & m_mask; // Mask to 32 bits only
  dd4hep::rec::SurfaceMap::const_iterator sI = surfaceMap->find(reduced_cellID);
  if (sI == surfaceMap->end()) {
    throw std::runtime_error(
        fmt::format("VTXdigitizerDetailed::operator(): no surface found for cellID : {}", reduced_cellID));
  }
  const dd4hep::rec::ISurface* surf = sI->second;
  // normal vector is by default in global frame
  dd4hep::rec::Vector3D normal_Global = surf->normal();
  TGeoHMatrix sensorTransformMatrix = m_volman.lookupDetElement(cellID).nominal().worldTransformation();
  double normal_Local_tmp[3] = {0, 0, 0};
  sensorTransformMatrix.MasterToLocalVect(normal_Global, normal_Local_tmp);
  dd4hep::rec::Vector3D normal_local(normal_Local_tmp[0], normal_Local_tmp[1], normal_Local_tmp[2]);

  // Determine which axis is the normal vector direction in local frame
  double epsilon = 1.0e-6; // A reasonable tolerance for floating-point comparisons

  if (std::abs(std::abs(normal_local.z()) - 1.0) < epsilon && std::abs(normal_local.x()) < epsilon &&
      std::abs(normal_local.y()) < epsilon) {
    m_LocalNormalVectorDir = "z";
  } else if (std::abs(std::abs(normal_local.y()) - 1.0) < epsilon && std::abs(normal_local.x()) < epsilon &&
             std::abs(normal_local.z()) < epsilon) {
    m_LocalNormalVectorDir = "y";
  } else if (std::abs(std::abs(normal_local.x()) - 1.0) < epsilon && std::abs(normal_local.y()) < epsilon &&
             std::abs(normal_local.z()) < epsilon) {
    m_LocalNormalVectorDir = "x";
  } else {
    error() << "Local normal vector is not aligned with local X, Y, or Z axis!" << endmsg;
  }

  debug() << "Local normal vector direction is: " << m_LocalNormalVectorDir << endmsg;

  return;
}

std::vector<VTXdigitizerDetailed::ChargeDepositUnit>
VTXdigitizerDetailed::primary_ionization(const edm4hep::SimTrackerHit& hit) const {
  /** Generate primary ionization along the track segment.
   *  Divide the track into small sub-segments of 1 microns
   *  Straight line approximation for trajectory of the incoming particle inside the active media
   *  The positions are given in mm
   */

  float pathLength = hit.getPathLength(); // Path Length of the particle in the active material in mm
  int numberOfSegments = int(pathLength / m_IonizationGranularity); // Number of segments
  if (numberOfSegments < 1) {
    numberOfSegments = 1;
  }

  float GeVperElectron = 3.61E-09; // Mean ionization energy in silicon [GeV]
  float charge;
  float eDep = hit.getEDep() * dd4hep::GeV; // Total energy deposited by the particle in this hit

  debug() << "Enter Primary_ionization, numberOfSegments=" << numberOfSegments << ", shift=" << pathLength
          << ", energy=" << eDep << "GeV" << endmsg;

  // Get the global position of the hit (defined by default in Geant4 as the mean between the entry and exit point in
  // the active material) and apply unit transformation (translation matrix is stored in cm
  double hitGlobalCentralPosition[3] = {hit.getPosition().x * dd4hep::mm, hit.getPosition().y * dd4hep::mm,
                                        hit.getPosition().z * dd4hep::mm};

  // Get the 3D momentum in global frames
  double hitGlobalMomentum[3] = {hit.getMomentum().x * dd4hep::GeV, hit.getMomentum().y * dd4hep::GeV,
                                 hit.getMomentum().z * dd4hep::GeV};

  // Convert the global position and momentum into local frame
  // retrieve the cell detElement
  const dd4hep::DDSegmentation::CellID& cellID = hit.getCellID();

  // Get transformation matrix of sensor and get central position and momentum vector in local frame
  // const auto& sensorTransformMatrix = m_volman.lookupVolumePlacement(cellID).matrix();
  TGeoHMatrix sensorTransformMatrix = m_volman.lookupDetElement(cellID).nominal().worldTransformation();
  double hitLocalCentralPosition[3] = {0, 0, 0};
  double hitLocalMomentum[3] = {0, 0, 0};

  // get the simHit coordinate in cm and momentum in the sensor reference frame
  sensorTransformMatrix.MasterToLocalVect(hitGlobalMomentum, hitLocalMomentum);
  SetLocalPos_In_ProperDirectFrame(hitLocalMomentum[0], hitLocalMomentum[1], hitLocalMomentum[2]);

  debug() << "Hit global momentum vector: (" << hitGlobalMomentum[0] << ", " << hitGlobalMomentum[1] << ", "
          << hitGlobalMomentum[2] << ")" << endmsg;
  debug() << "Hit local momentum vector: (" << hitLocalMomentum[0] << ", " << hitLocalMomentum[1] << ", "
          << hitLocalMomentum[2] << ")" << endmsg;

  sensorTransformMatrix.MasterToLocal(hitGlobalCentralPosition, hitLocalCentralPosition);
  SetLocalPos_In_ProperDirectFrame(hitLocalCentralPosition[0], hitLocalCentralPosition[1], hitLocalCentralPosition[2]);

  // build vectors to simplify next calculation and come back to mm
  dd4hep::rec::Vector3D hitLocalCentralPositionVec(hitLocalCentralPosition[0] / dd4hep::mm,
                                                   hitLocalCentralPosition[1] / dd4hep::mm,
                                                   hitLocalCentralPosition[2] / dd4hep::mm);
  dd4hep::rec::Vector3D hitLocalMomentumVec(hitLocalMomentum[0] / dd4hep::mm, hitLocalMomentum[1] / dd4hep::mm,
                                            hitLocalMomentum[2] / dd4hep::mm);

  // Get particle direction in local frame normalized to path Length in mm
  dd4hep::rec::Vector3D direction = (pathLength / hitLocalMomentumVec.r()) * hitLocalMomentumVec;

  // Get the entry point in the active material in mm
  dd4hep::rec::Vector3D entryPoint = hitLocalCentralPositionVec - 0.5 * direction;

  // Maybe add a test that the point is inside the material or go to the closest material and do the same for exit point
  // and then redefine the length

  // Calulate the charge deposited per segment
  std::vector<ChargeDepositUnit> ionizationPoints;
  ionizationPoints.resize(numberOfSegments);

  for (int i = 0; i != numberOfSegments; i++) {
    // Divide the segment into equal length sub-segments in mm
    dd4hep::rec::Vector3D point = entryPoint + float((i + 0.5) / numberOfSegments) * direction;

    charge = eDep / GeVperElectron / float(numberOfSegments); // Here, do not consider energy fluctuation. To consider
                                                              // it, see the Geant4 Physics Reference Manual

    ChargeDepositUnit cdu(charge, point); // Define position(mm),charge(Number of electrons) point
    ionizationPoints[i] = cdu;            // save

  } // end loop over segments

  return ionizationPoints;
} // End primary_ionization

std::vector<VTXdigitizerDetailed::SignalPoint>
VTXdigitizerDetailed::drift(const edm4hep::SimTrackerHit& hit,
                            const std::vector<ChargeDepositUnit>& ionizationPoints) const {
  /** Drift the charge segments to the sensor surface (collection plane)
   *  Include the effect of E-field and B-field
   *  Even for strips, consider a uniform E-field
   *  The sensor depth is explicitly considered to be along z
   *  For now, the drift is considered to be instantaneous. Later On, should consider time and time uncertainty
   *  The positions are in mm
   */

  debug() << "Enter drift " << endmsg;

  std::vector<SignalPoint> collectionPoints;
  collectionPoints.resize(ionizationPoints.size()); // set size

  dd4hep::rec::Vector3D driftDir = DriftDirection(hit); // Get the charge (electrons) drift direction
  if (driftDir.z() == 0.) {
    warning() << "Charge drift in z is 0" << endmsg;
    return collectionPoints;
  }

  float tanLorentzAngleX = driftDir.x(); // tan of Lorentz angle
  float tanLorentzAngleY = driftDir.y();
  float dir_z = driftDir.z();                                                   // The z drift direction
  float cosLorentzAngleX = 1. / sqrt(1. + tanLorentzAngleX * tanLorentzAngleX); // Cosine
  float cosLorentzAngleY = 1. / sqrt(1. + tanLorentzAngleY * tanLorentzAngleY);

  // Get the sensor thickness
  const dd4hep::DDSegmentation::CellID& cellID = hit.getCellID();
  int layerID = m_decoder->get(cellID, "layer");
  double moduleThickness = m_sensorThickness[layerID]; // In mm

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

    if (m_DebugHistos) {
      hXDriftDueToMagField->Fill(XDriftDueToMagField);
      hYDriftDueToMagField->Fill(YDriftDueToMagField);
    }

    // Shift cloud center
    float CloudCenterX = SegX + XDriftDueToMagField;
    float CloudCenterY = SegY + YDriftDueToMagField;

    // Calculate how long is the charge drift path
    DriftLength = sqrt(DriftDistance * DriftDistance + XDriftDueToMagField * XDriftDueToMagField +
                       YDriftDueToMagField * YDriftDueToMagField);

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

  return collectionPoints;

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

  // Get hit global position in mm (magnetic Field mapped in mm)
  dd4hep::rec::Vector3D hitPosition(hit.getPosition().x * dd4hep::mm, hit.getPosition().y * dd4hep::mm,
                                    hit.getPosition().z * dd4hep::mm);

  // Get B field in global frame
  const auto& magneticField = detector->field();
  double BFieldGlobal[3];
  magneticField.magneticField(hitPosition, BFieldGlobal);

  // Get B field in local frame
  double BFieldLocal[3] = {0, 0, 0};
  const dd4hep::DDSegmentation::CellID& cellID = hit.getCellID();
  auto sensorTransformMatrix = m_volman.lookupDetElement(cellID).nominal().worldTransformation();
  sensorTransformMatrix.MasterToLocalVect(BFieldGlobal, BFieldLocal);
  SetLocalPos_In_ProperDirectFrame(BFieldLocal[0], BFieldLocal[1], BFieldLocal[2]);

  // Retrieve local magnetic field in a vector and in Tesla
  dd4hep::rec::Vector3D BFieldLocalVec(BFieldLocal[0] / dd4hep::tesla, BFieldLocal[1] / dd4hep::tesla,
                                       BFieldLocal[2] / dd4hep::tesla);

  // Get the drift direction
  float dir_x = 0.;
  float dir_y = 0.;
  float dir_z = 0.;
  float scale = 0.;
  float alpha2 = m_tanLorentzAnglePerTesla * m_tanLorentzAnglePerTesla;

  // get the directions considering E-field along z and B field (Lorentz angle)
  // Consider here that electrons go to +z (E field in -z direction). Signs to be change if reverse, accordingly to CMS
  // digitizer.
  dir_x = -m_tanLorentzAnglePerTesla * BFieldLocalVec.y() + alpha2 * BFieldLocalVec.z() * BFieldLocalVec.x();
  dir_y = m_tanLorentzAnglePerTesla * BFieldLocalVec.x() + alpha2 * BFieldLocalVec.z() * BFieldLocalVec.y();
  dir_z = 1 + alpha2 * BFieldLocalVec.z() * BFieldLocalVec.z();
  scale = dir_z; // Normalize to z direction

  dd4hep::rec::Vector3D theDriftDirection(dir_x / scale, dir_y / scale, dir_z / scale);

  return theDriftDirection;

} // End DriftDirection

VTXdigitizerDetailed::hit_map_type
VTXdigitizerDetailed::get_charge_per_pixel(const edm4hep::SimTrackerHit& hit,
                                           const std::vector<SignalPoint>& collectionPoints) const {

  /** Get the map of recorded charges per pixel for a collection of drifted charges for a given hit
   */

  // Initialize the map
  hit_map_type hit_map;

  const dd4hep::DDSegmentation::CellID& cellID = hit.getCellID();
  int ilayer = m_decoder->get(cellID, "layer");

  const float PixSizeX = m_PixSizeX[ilayer]; // in mm
  const float PixSizeY = m_PixSizeY[ilayer]; // in mm

  std::uint64_t m_mask = (static_cast<std::uint64_t>(1) << 32) - 1;
  const std::uint64_t reduced_cellID = cellID & m_mask;

  const auto solid = m_volman.lookupDetElement(reduced_cellID).volume().solid();

  double dimX, dimY, dimZ; // Dimensions of the solid in cm
  if (std::string(solid.type()) == "TGeoBBox") {
    dd4hep::Box box(solid);
    dimX = 2 * box.x(); // Get the dimensions in cm
    dimY = 2 * box.y();
    dimZ = 2 * box.z();
    SetLocalPos_In_ProperDirectFrame(dimX, dimY, dimZ); // Set the local position in the proper direct frame
    dimX = fabs(dimX);
    dimY = fabs(dimY);
    dimZ = fabs(dimZ); // Ensure positive dimensions
  } else if (std::string(solid.type()) == "TGeoTrd2") {
    // We consider only the largest width and length (trapezoid considered as largest rectangle)
    dd4hep::Trapezoid trap(solid);
    dimX = 2 * std::max(trap.dX1(), trap.dX2()); // Get the dimensions in cm
    dimY = 2 * std::max(trap.dY1(), trap.dY2()); // Get the dimensions in cm
    dimZ = 2 * trap.dZ();
    SetLocalPos_In_ProperDirectFrame(dimX, dimY, dimZ); // Set the local position in the proper direct frame
    dimX = fabs(dimX);
    dimY = fabs(dimY);
    dimZ = fabs(dimZ); // Ensure positive dimensions
  } else {
    error() << "Solid type not recognized: " << solid.type() << endmsg;
  }

  float widthMin = -dimX / 2. / dd4hep::mm, widthMax = dimX / 2. / dd4hep::mm;
  float lengthMin = -dimY / 2. / dd4hep::mm, lengthMax = dimY / 2. / dd4hep::mm;

  debug() << "Sensor size: width [" << widthMin << ", " << widthMax << "] mm, length [" << lengthMin << ", "
          << lengthMax << "] mm" << endmsg;

  // Mapping the sensor's physical boundaries to pixel indices
  const int MinPixXSensor = static_cast<int>(std::floor((widthMin + 0.5f * PixSizeX) / PixSizeX));
  const int MaxPixXSensor = static_cast<int>(std::floor((widthMax + 0.5f * PixSizeX) / PixSizeX));
  const int MinPixYSensor = static_cast<int>(std::floor((lengthMin + 0.5f * PixSizeY) / PixSizeY));
  const int MaxPixYSensor = static_cast<int>(std::floor((lengthMax + 0.5f * PixSizeY) / PixSizeY));

  debug() << "Sensor bounds in Pixel Indices : X [" << MinPixXSensor << ", " << MaxPixXSensor << "], Y ["
          << MinPixYSensor << ", " << MaxPixYSensor << "]" << endmsg;

  // map to store the pixel integrals in the x and y directions
  std::map<int, float> x, y;

  // Assign signal per readout channel and store sorted by channel number (in modified local frame)
  // Iterate over collection points on the collection plane
  for (const auto& point : collectionPoints) {
    float CloudCenterX = point.x();   // Charge position in x in mm
    float CloudCenterY = point.y();   // Charge position in y
    float SigmaX = point.sigma_x();   // Charge spread in x in mm
    float SigmaY = point.sigma_y();   // Charge spread in y
    float Charge = point.amplitude(); // Charge amplitude in number of e-

    // Find the maximum cloud spread in 2D Plane, assume 3*sigma
    float CloudMinX = CloudCenterX - m_ClusterWidth * SigmaX;
    float CloudMaxX = CloudCenterX + m_ClusterWidth * SigmaX;
    float CloudMinY = CloudCenterY - m_ClusterWidth * SigmaY;
    float CloudMaxY = CloudCenterY + m_ClusterWidth * SigmaY;

    // Convert the maximum cloud spread into pixel IDs
    // Considers that the pixel of indice (0,0) is centered in position (0,0)
    int IpxCloudMinX = int(floor((CloudMinX + 0.5 * PixSizeX) / PixSizeX));
    int IpxCloudMaxX = int(floor((CloudMaxX + 0.5 * PixSizeX) / PixSizeX));
    int IpxCloudMinY = int(floor((CloudMinY + 0.5 * PixSizeY) / PixSizeY));
    int IpxCloudMaxY = int(floor((CloudMaxY + 0.5 * PixSizeY) / PixSizeY));

    x.clear(); // clear temporary integration arrays
    y.clear();

    // Integrate charge strips in x
    for (int ix = IpxCloudMinX; ix <= IpxCloudMaxX; ix++) {

      float LowerBound = (float(ix) - 0.5) * PixSizeX; // Lower bound of the strip
      float UpperBound = (float(ix) + 0.5) * PixSizeX; // Uper bound of the strip

      float TotalStripCharge =
          0.5 * (erf((UpperBound - CloudCenterX) / (sqrt(2) * SigmaX)) -
                 erf((LowerBound - CloudCenterX) /
                     (sqrt(2) * SigmaX))); // Charge proportion in the strip calculated from erf function

      x[ix] = TotalStripCharge;

    } // End Integrate charge strips in x

    // Integrate charge strips in y
    for (int iy = IpxCloudMinY; iy <= IpxCloudMaxY; iy++) {

      float LowerBound = (float(iy) - 0.5) * PixSizeY; // Lower bound of the strip
      float UpperBound = (float(iy) + 0.5) * PixSizeY; // Uper bound of the strip

      float TotalStripCharge =
          0.5 * (erf((UpperBound - CloudCenterY) / (sqrt(2) * SigmaY)) -
                 erf((LowerBound - CloudCenterY) /
                     (sqrt(2) * SigmaY))); // Charge proportion in the strip calculated from erf function

      y[iy] = TotalStripCharge;

    } // End Integrate charge strips in y

    // Get the 2D charge integrals by folding x and y strips
    for (int ix = IpxCloudMinX; ix <= IpxCloudMaxX; ix++) {
      for (int iy = IpxCloudMinY; iy <= IpxCloudMaxY; iy++) {

        // Check if the pixel is within the sensor bounds
        if (ix < MinPixXSensor || ix > MaxPixXSensor || iy < MinPixYSensor || iy > MaxPixYSensor) {
          debug() << "Pixel out of bounds (skipped) : " << ix << ":" << iy << endmsg;
          continue; // Skip pixels outside the sensor bounds
        }

        // Calculate the charge in the pixel
        // ChargeInPixel = Charge * x[ix] * y[iy]
        // where x[ix] and y[iy] are the charge fractions in the x and y directions
        float ChargeInPixel = Charge * x[ix] * y[iy];

        hit_map[ix][iy] += ChargeInPixel; // Add the charge to the pixel map
      } // end loop over y
    } // end loop over x
  } // End loop over charge collection

  return hit_map;
} // End get_charge_per_pixel

bool VTXdigitizerDetailed::Apply_Threshold(double& ChargeInE) const {
  /** Apply the threshold to the charge
   *  The threshold is defined in electrons
   *  The threshold is smeared with a Gaussian distribution if m_ThresholdSmearing > 0
   *  The threshold is applied to the charge
   *  If the charge is above the threshold, return true
   *  If the charge is below the threshold, return false
   */

  // Get the threshold in electrons
  double ThresholdInE = 0;

  if (m_ThresholdSmearing > 0) {
    ThresholdInE = m_gauss_threshold();
  } else {
    ThresholdInE = m_Threshold;
  }

  // In order to avoid negative threshold values, we set the minimum value to 0
  ThresholdInE = std::max(ThresholdInE, 0.0);

  debug() << "ChargeInE (weight): " << ChargeInE << ", ThresholdInE (after smearing): " << ThresholdInE << endmsg;

  return ChargeInE >= ThresholdInE;
}

void VTXdigitizerDetailed::generate_output(const edm4hep::SimTrackerHit& hit,
                                           edm4hep::TrackerHitPlaneCollection& output_digi_hits,
                                           edm4hep::TrackerHitSimTrackerHitLinkCollection& output_sim_digi_link_col,
                                           const hit_map_type& hit_map) const {
  /** Get a single point from the hit_map by getting the average position and time
   *  Return the position and uncertainty in global coordinates with recorded energy and time as a
   * TrackerhitplaneCollection For now, the position is taken from a weighted average, maybe try with a 2D Gaussian fit
   */

  // Get the pixel dimensions in mm along x and y in the (modified) local frame with z orthogonal

  const dd4hep::DDSegmentation::CellID& cellID = hit.getCellID();
  int ilayer = m_decoder->get(cellID, "layer");

  const float PixSizeX = m_PixSizeX[ilayer]; // in mm
  const float PixSizeY = m_PixSizeY[ilayer]; // in mm

  debug() << "SimHit CellID: " << cellID << endmsg;
  for (const auto& field : m_decoder->fields()) {
    const std::string& name = field.name();
    int value = m_decoder->get(cellID, name);
    debug() << "  " << name << " = " << value << endmsg;
  }
  debug() << "Hit position in mm : (" << hit.getPosition().x << ", " << hit.getPosition().y << ", "
          << hit.getPosition().z << ")" << endmsg;
  debug() << "is Secondary ? : " << hit.isProducedBySecondary() << endmsg;
  debug() << "is Overlay ? : " << hit.isOverlay() << endmsg;

  // Local position of the digitized hit
  double DigiLocalX = 0.;
  double DigiLocalY = 0.;

  using ChargeMap = std::unordered_map<int, double>;
  ChargeMap Qix, Qiy;

  double sumWeights = 0.; // Sum of all weights (sum of charges recorded in all pixels)
  int CountBeforeThreshold = 0;
  int CountAfterThreshold = 0;
  int nPixels = 0;

  debug() << "------ Pixel map per cluster (hit_map) ------" << endmsg;
  // Loop over hit map to apply threshold and accumulate weights
  // loop to load the weights per x and y layers
  for (auto const& ix : hit_map) {
    for (auto const& iy : ix.second) {

      double weight = iy.second;
      int x = ix.first;
      int y = iy.first;
      debug() << "Pixel (" << x << ", " << y << ") : Charge = " << weight << endmsg;
      nPixels++;
      if (m_DebugHistos) {
        ++CountBeforeThreshold;
        hChargeBeforeThreshold->Fill(weight / 1e3); // Fill the histogram with the charge before threshold in (ke)
      }

      if (weight > 0 && Apply_Threshold(weight)) // Check if the charge is above the threshold
      {

        if (m_DebugHistos) {
          ++CountAfterThreshold;
          hChargeAboveThreshold->Fill(weight / 1e3); // Fill the histogram with the charge above threshold in (ke)
        }

        Qix[ix.first] += weight;
        Qiy[iy.first] += weight;
        sumWeights += weight;
      }
    }
  } // end loop over y
  debug() << "--------------------------------------" << endmsg;
  debug() << "Total number of pixels with signal : " << nPixels << endmsg;
  if (m_DebugHistos) {
    hActivePixelCountBeforeThreshold->Fill(CountBeforeThreshold);
  }

  if (sumWeights <= 0.) {
    debug() << "sumWeights <= 0 Skip output." << endmsg;
    return;
  }

  if (m_DebugHistos) {
    hActivePixelCountAfterThreshold->Fill(CountAfterThreshold);
    hChargePerCluster->Fill(sumWeights / 1e3);
    hActivePixelPerlayer->Fill(ilayer, CountAfterThreshold);
  }

  double sumWeightsSqX = 0.; // Sum of the square of weights along x (used later for position resolution)
  double sumWeightsSqY = 0.;

  double pixLocalX, pixLocalY; // Initiate pixels position
  // Loop to get the weighted average of the X position
  for (auto const& ix : Qix) {

    pixLocalX =
        ix.first *
        PixSizeX; // Pixel local X position // Considers that the pixel of indice (0,0) is centered in position (0,0)
    DigiLocalX += ix.second * pixLocalX;
    sumWeightsSqX += ix.second * ix.second; // Load the square sum of the weights

  } // loop over x

  DigiLocalX = DigiLocalX / sumWeights;

  // Loop to get the weighted average of the Y position
  for (auto const& iy : Qiy) {

    pixLocalY = iy.first * PixSizeY; // Pixel local Y position
    DigiLocalY += iy.second * pixLocalY;
    sumWeightsSqY += iy.second * iy.second; // Load the square sum of the weights

  } // loop over y

  DigiLocalY = DigiLocalY / sumWeights;

  // Get the position resolution in Local frame
  double resU = (PixSizeX / sqrt(12)) * sqrt(sumWeightsSqX) /
                sumWeights; // Position uncertainty along X. Consider Uniform distribution of charges in a pixel
  double resV = (PixSizeY / sqrt(12)) * sqrt(sumWeightsSqY) / sumWeights;

  double dirLocalX[3] = {1., 0., 0.}; // X direction in local frame
  double dirLocalY[3] = {0., 1., 0.}; // Y direction in local frame
  double dirGlobalU[3];               // Global direction U (corresponding to local X)
  double dirGlobalV[3];               // Global direction V (cooresponding to local Y)
  double DigiLocalPos[3] = {DigiLocalX * dd4hep::mm, DigiLocalY * dd4hep::mm,
                            0.}; // Local Position in cm, necessary for transformation to global with transform matrix
  double DigiGlobalPos[3];       // Global position in cm

  debug() << "DigiLocalX: " << DigiLocalPos[0] << ", DigiLocalY: " << DigiLocalPos[1] << endmsg;

  // Get the transformation matrix between local and global frames
  auto sensorTransformMatrix = m_volman.lookupDetElement(cellID).nominal().worldTransformation();

  // Transfrom the local positions and directions in ProperDirectFrame (z normal to sensor) to DD4Hep Standard Frame
  SetLocalPos_In_StandardDD4hepFrame(dirLocalX[0], dirLocalX[1], dirLocalX[2]);
  SetLocalPos_In_StandardDD4hepFrame(dirLocalY[0], dirLocalY[1], dirLocalY[2]);
  SetLocalPos_In_StandardDD4hepFrame(DigiLocalPos[0], DigiLocalPos[1], DigiLocalPos[2]);

  // Transform the local positions and directions to global ones
  sensorTransformMatrix.LocalToMasterVect(dirLocalX, dirGlobalU);
  sensorTransformMatrix.LocalToMasterVect(dirLocalY, dirGlobalV);
  sensorTransformMatrix.LocalToMaster(DigiLocalPos, DigiGlobalPos);

  // Go back to mm
  edm4hep::Vector3d DigiGlobalPosVec(DigiGlobalPos[0] / dd4hep::mm, DigiGlobalPos[1] / dd4hep::mm,
                                     DigiGlobalPos[2] / dd4hep::mm);

  // Get global directions of the sensor for the HitPlane
  dd4hep::rec::Vector3D dirGlobalUVec(dirGlobalU[0], dirGlobalU[1], dirGlobalU[2]);
  dd4hep::rec::Vector3D dirGlobalVVec(dirGlobalV[0], dirGlobalV[1], dirGlobalV[2]);

  float u_direction[2];
  u_direction[0] = dirGlobalUVec.theta();
  u_direction[1] = dirGlobalUVec.phi();

  float v_direction[2];
  v_direction[0] = dirGlobalVVec.theta();
  v_direction[1] = dirGlobalVVec.phi();

  // Create and fill the output collection of digitized hits
  auto output_digi_hit = output_digi_hits.create();
  auto output_sim_digi_link = output_sim_digi_link_col.create();

  float GeVperElectron = 3.61E-09;                                    // Mean ionization energy in silicon [GeV]
  output_digi_hit.setEDep(sumWeights * GeVperElectron / dd4hep::GeV); // Energy in GeV from number of charges
  // setEDepError ?

  output_digi_hit.setPosition(DigiGlobalPosVec);

  // Apply time smearing
  output_digi_hit.setTime(hit.getTime() + m_gauss_t_vec[ilayer].shoot());

  output_digi_hit.setCellID(cellID);

  // Set sensor directions
  output_digi_hit.setU(u_direction);
  output_digi_hit.setV(v_direction);

  // Set the position uncertainty
  output_digi_hit.setDu(resU);
  output_digi_hit.setDv(resV);

  // Set the link between sim and digi hit
  output_sim_digi_link.setFrom(output_digi_hit);
  output_sim_digi_link.setTo(hit);

  debug() << "Output Digi hit : " << endmsg;
  debug() << "EDep : " << output_digi_hit.getEDep() << endmsg;
  debug() << "Position (mm) : (" << output_digi_hit.getPosition().x << ", " << output_digi_hit.getPosition().y << ", "
          << output_digi_hit.getPosition().z << ")" << endmsg;
  debug() << "Time (ns) : " << output_digi_hit.getTime() << endmsg;
  debug() << "CellID : " << output_digi_hit.getCellID() << endmsg;
  debug() << "Vector U Theta : " << u_direction[0] << endmsg;
  debug() << "Vector U Phi : " << u_direction[1] << endmsg;
  debug() << "Vector V Theta : " << v_direction[0] << endmsg;
  debug() << "Vector V Phi : " << v_direction[1] << endmsg;
  debug() << "Du (mm) : " << output_digi_hit.getDu() << endmsg;
  debug() << "Dv (mm) : " << output_digi_hit.getDv() << endmsg;
  // Fill the debug histograms if needed
  if (m_DebugHistos) {

    // Get the global position of the hit (defined by default in Geant4 as the mean between the entry and exit point in
    // the active material) and apply unit transformation (translation matrix is stored in cm)
    double hitGlobalCentralPosition[3] = {hit.getPosition().x * dd4hep::mm, hit.getPosition().y * dd4hep::mm,
                                          hit.getPosition().z * dd4hep::mm};

    double hitLocalCentralPosition[3] = {0, 0, 0};
    // get the simHit coordinate in cm in the sensor reference frame
    sensorTransformMatrix.MasterToLocal(hitGlobalCentralPosition, hitLocalCentralPosition);
    SetLocalPos_In_ProperDirectFrame(hitLocalCentralPosition[0], hitLocalCentralPosition[1],
                                     hitLocalCentralPosition[2]);

    double DistX = DigiLocalX - hitLocalCentralPosition[0] / dd4hep::mm;
    double DistY = DigiLocalY - hitLocalCentralPosition[1] / dd4hep::mm;
    double DistZ = 0. - hitLocalCentralPosition[2] / dd4hep::mm;
    hErrorX->Fill(DistX);
    hErrorY->Fill(DistY);
    hErrorZ->Fill(DistZ);
    hError->Fill(sqrt(DistX * DistX + DistY * DistY + DistZ * DistZ));

  } // End Debug Histos

} // End get_hist_signal_point

void VTXdigitizerDetailed::SetLocalPos_In_ProperDirectFrame(double& x, double& y, double& z) const {
  /** Change coordinates of a point to have a direct frame with z orthogonal to sensor surface
   */

  std::string LocalNormalVectorDir = m_LocalNormalVectorDir;

  double x_tmp, y_tmp, z_tmp;

  // If the orthogonal direction is along X or Y in local frame, rotate the frame to have Z orthogonal to sensors
  // instead
  if (LocalNormalVectorDir == 'x') {
    // X->Z / Y->X / Z->Y
    z_tmp = x;
    x_tmp = y;
    y_tmp = z;
  } else if (LocalNormalVectorDir == 'y') {
    // X->X / Y->Z / Z->-Y
    z_tmp = -y;
    x_tmp = x;
    y_tmp = z;
  } else {
    // X->X / Y->Y / Z->Z
    z_tmp = z;
    x_tmp = x;
    y_tmp = y;
  }

  x = x_tmp;
  y = y_tmp;
  z = z_tmp;

} // SetLocalPos_In_ProperDirectFrame

void VTXdigitizerDetailed::SetLocalPos_In_StandardDD4hepFrame(double& x, double& y, double& z) const {
  /** Change coordinates of a point in Proper Direct Frame with z orthogonal to sensor surface
   *  back to the standard DD4hep local frame
   */

  std::string LocalNormalVectorDir = m_LocalNormalVectorDir;

  double x_tmp = x, y_tmp = y, z_tmp = z;

  // If the orthogonal direction is along X or Y in local frame, rotate the frame to have the good axie orthogonal to
  // sensors
  if (LocalNormalVectorDir == 'x') {
    // Z->X / X->Y / Y->Z
    x_tmp = z;
    y_tmp = x;
    z_tmp = y;
  } else if (LocalNormalVectorDir == 'y') {
    // X->X / Z->Y / Y->-Z
    x_tmp = x;
    y_tmp = -z;
    z_tmp = y;
  } else {
    x_tmp = x;
    y_tmp = y;
    z_tmp = z;
  }

  x = x_tmp;
  y = y_tmp;
  z = z_tmp;
} // SetLocalPos_In_StandardDD4hepFrame

void VTXdigitizerDetailed::Create_outputROOTfile_for_debugHistograms() const {
  /** This is an internal function to save the debug histograms in the corresponding rootfile
   */

  // save current ROOT directory
  TDirectory* currentDir = gDirectory;

  // save the debug histograms in a file
  // file is saved and closed when going out of scope
  auto filename = m_DebugFileName.value().c_str();
  std::unique_ptr<TFile> ofile{TFile::Open(filename, "recreate")};
  if (!ofile || ofile->IsZombie()) {
    error() << "Error: Could not open file " << filename << std::endl;
    return;
  }
  ofile->cd();
  hErrorX->Write();
  hErrorY->Write();
  hErrorZ->Write();
  hError->Write();
  hXDriftDueToMagField->Write();
  hYDriftDueToMagField->Write();

  hChargeBeforeThreshold->Write();
  hActivePixelCountBeforeThreshold->Write();
  hChargeAboveThreshold->Write();
  hActivePixelCountAfterThreshold->Write();
  hChargePerCluster->Write();
  hClusterPerLayer->Write();
  hActivePixelPerlayer->Write();

  // Restore previous ROOT directory
  if (currentDir && (not currentDir->IsDestructed()))
    currentDir->cd();
  return;

} // End Create_outputROOTfile_for_debugHistograms
