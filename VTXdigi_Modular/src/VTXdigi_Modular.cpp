// VTXdigi_Modular/src/VTXdigi_Modular.cpp
#include "VTXdigi_Modular.h"

DECLARE_COMPONENT(VTXdigi_Modular)

VTXdigi_Modular::VTXdigi_Modular(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValues("SimTrackHitCollectionName", {"UNDEFINED_SimTrackHitCollectionName"}),
                        KeyValues("HeaderName", {"UNDEFINED_HeaderName"}),},
                       {KeyValues("TrackerHitCollectionName", {"UNDEFINED_TrackerHitCollectionName"}),
                        KeyValues("SimTrkHitRelationsCollection", {"UNDEFINED_SimTrkHitRelationsCollection"})}) {
  info() << "Constructed successfully" << endmsg;
}

StatusCode VTXdigi_Modular::initialize() {
  info() << "INITIALIZING VTXdigi_Modular..." << endmsg;

  info() << "OutputLevel set to " << msgSvc()->outputLevel(name()) << endmsg;
  // TODO: implement if-clause for debug/verbose messages in hot loops (to avoid constructing the message string when the message won't be printed) -> this will improve performance significantly

  // if (!VTXdigi_tools::ToolTest())
  //   return StatusCode::FAILURE;

  InitServicesAndGeometry();

  InitLayersAndSensors();

  if (m_debugHistograms)
    InitHistograms();

  // This needs to come in after the properties, geometry and services have all been initialized 
  verbose() << "Initializing charge collection method: " << m_chargeCollectionMethod.value() << endmsg;
  m_chargeCollector = VTXdigi_tools::CreateChargeCollector(*this, m_chargeCollectionMethod);

  info() << " - Initialized successfully." << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode VTXdigi_Modular::finalize() {
  info() << "FINALIZING VTXdigi_Modular..." << endmsg;

  PrintCountersSummary();

  debug() << " - finalized successfully." << endmsg;
  return StatusCode::SUCCESS;
}


/* ---- Event loop ---- */

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> VTXdigi_Modular::operator()
  (const edm4hep::SimTrackerHitCollection& simTrackerHits, const edm4hep::EventHeaderCollection& headers) const {
  if (!CheckEventSetup(simTrackerHits, headers)) {
    return std::make_tuple(edm4hep::TrackerHitPlaneCollection(), edm4hep::TrackerHitSimTrackerHitLinkCollection());
  }

  /* TODO: Implement fast digitization, where simHit pos is simply smeared by a Gaussian.
  * this makes a lot of the init etc unnecessary, so maybe keep it in a separate class? We cannot simply add a ChargeCollector implementation to do this, because ChargeCollectors act on pixels, but fast digitization acts on the position itself (ofc we could emulate this in clusters, but thats incredibly inefficient and stinks) */

  std::unordered_map<dd4hep::DDSegmentation::CellID, std::vector<VTXdigi_tools::SimHitWrapper>> sensorSimHits; // map from sensor (shortened cellID) to hits
  for (const auto& simTrackerHit : simTrackerHits) {
    if (CheckSimhitLayer(simTrackerHit)){
      const dd4hep::DDSegmentation::CellID cellID = VTXdigi_tools::GetCellID_short(simTrackerHit);
      sensorSimHits[cellID].emplace_back(simTrackerHit, m_cellIdDecoder); // simTrackerHits are copied here. Pointers to these are passed around (eg. in hit/pixel/cluster objects).
    }
  }
  debug() << " - Found simTrackerHits on " << sensorSimHits.size() << " individual sensors. Digitising sensor-wise..." << endmsg;

  auto digiHits = edm4hep::TrackerHitPlaneCollection();
  auto digiHitLinks = edm4hep::TrackerHitSimTrackerHitLinkCollection();
  std::vector<VTXdigi_tools::SimHitWrapper> hitsSensor;
  
  /* loop over sensors */
  for (const auto& [cellID, simHits] : sensorSimHits) {
    debug() << "   - Processing sensor with cellID " << cellID << " (layer " << simHits.back().layer() << "). Has " << simHits.size() << " simHits." << endmsg;
  
    const TGeoHMatrix trafoMatrix = VTXdigi_tools::ComputeSensorTrafoMatrix(cellID, m_volumeManager, m_sensorNormalRotation); // transformation from global detector to local sensor frame
    
    /* Fill hit map with charge from all simHits on this sensor, using the selected charge collection method */
    VTXdigi_tools::HitMap hitMap(m_pixelCount);

    for (const VTXdigi_tools::SimHitWrapper& simHit : simHits) {
      const dd4hep::rec::Vector3D pos_global = VTXdigi_tools::ConvertVector(simHit.hitPtr()->getPosition());
      simHit.SetTruthPos(VTXdigi_tools::GlobalToLocal(pos_global, trafoMatrix)); // do this only now to not compute the trafo matrix twice. TruthPos might be shifted by the charge collection algorithm later.
      debug() << "     - Processing simHit, charge dep. " << simHit.charge() << " e at (" << simHit.truthPos().x() << ", " << simHit.truthPos().y() << ", " << simHit.truthPos().z() << ") local, (" << pos_global.x() << ", " << pos_global.y() << ", " << pos_global.z() << ") global" << endmsg;

      m_chargeCollector->FillHit(simHit, hitMap, trafoMatrix); // uses the selected charge collection method

      if (m_debugHistograms.value())
        FillHistograms_perSimHit(simHit);
    }

    
    if (m_smearing_charge.value() > 0.f) 
      hitMap.ApplyChargeSmearing(m_rndm_charge);
    if (m_threshold.value() > 0.f)
      hitMap.ApplyThreshold(m_threshold.value());

    std::vector<VTXdigi_tools::Cluster> clusters = Clusterize(hitMap);

    CreateDigiHits(digiHits, digiHitLinks, cellID, trafoMatrix, clusters);
  } /* loop over sensors */
  
  debug() << " - Finished digitization. Created " << digiHits.size() << " digiHits from " << simTrackerHits.size() << " simTrackerHits." << endmsg;
  return std::make_tuple(std::move(digiHits), std::move(digiHitLinks));
} // operator()

/* ---- Initialization & finalization functions ---- */

void VTXdigi_Modular::InitServicesAndGeometry() {
  /* Sets the members: 
   *  - m_uidService
   *  - m_geoService
   *  - m_cellIdDecoder
   *  - m_detector
   *  - m_surfaceMap
   *  - m_volumeManager
   *  - m_subDetector
   */
  if (m_threshold.value() < 0.f)
    throw GaudiException("Threshold " + std::to_string(m_threshold.value()) + " e- is negative.", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
  if (m_smearing_charge.value() < 0.f)
    throw GaudiException("Charge smearing sigma " + std::to_string(m_smearing_charge.value()) + " e- is negative.", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  if (m_threshold.value() <= 5 * m_smearing_charge.value())
    warning() << "Threshold " << m_threshold.value() << " e- is less than 5 times the charge smearing sigma " << m_smearing_charge.value() << " e-. This digitiser does only apply smearing to pixels that collect any charge from simHits, so it cannot simulate random firing pixels. (doing this by drawing a noise for every pixel in the detector for every event would be INCREDIBLY slow. A work-around to simulate random pixels firing might be implemented in the future)." << endmsg;
  if (m_smearing_time.value() < 0.f)
    throw GaudiException("Time smearing sigma " + std::to_string(m_smearing_time.value()) + " ns is negative.", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  if (m_positionUncertainty.value().empty())
    info() << "Cluster position uncertainty not set. Doing simple charge-weighted uncertainty estimation." << endmsg;
  else if (m_positionUncertainty.value().size() == 2)
    info() << "Cluster position uncertainty set to (" << m_positionUncertainty.value().at(0) << " mm, " << m_positionUncertainty.value().at(1) << " mm) in u and v direction." << endmsg;
  else if (m_positionUncertainty.value().size() == 6)
    info() << "Cluster position uncertainty set to (" << m_positionUncertainty.value().at(0) << ", " << m_positionUncertainty.value().at(1) << ", " << m_positionUncertainty.value().at(2) << ") mm and (" << m_positionUncertainty.value().at(3) << ", " << m_positionUncertainty.value().at(4) << ", " << m_positionUncertainty.value().at(5) << ") mm for cluster lengths of (1, 2, 3), in u and v direction, respectively." << endmsg;
  else
    throw GaudiException("Property ClusterPositionUncertainty must be either empty (for charge-weighted estimation), have exactly 2 values (for fixed uncertainty in u and v), or six values (for cluster-length based estimation).", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  m_randomService = service("RndmGenSvc", false);
  if (!m_randomService) 
    throw GaudiException("Unable to get RndmGenSvc.", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  if (m_rndm_charge.initialize(m_randomService, Rndm::Gauss(0., m_smearing_charge.value())).isFailure())
    throw GaudiException("Unable to initialize random number generator for charge smearing.", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
  if (m_rndm_time.initialize(m_randomService, Rndm::Gauss(0., m_smearing_time.value())).isFailure())
    throw GaudiException("Unable to initialize random number generator for time smearing.", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
  
  if (m_subDetName.value() == m_undefinedString)
    throw GaudiException("Property SubDetectorName is not set!", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  m_geoService = serviceLocator()->service(m_geoServiceName);
  if (!m_geoService)
    throw GaudiException("Unable to retrieve the GeoSvc", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
  
  /* TODO: get cellID without using geoService, a la Jessy (if this is advantageous) */
  std::string cellIDstr = m_geoService->constantAsString(m_encodingStringVariable.value());
  m_cellIdDecoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(cellIDstr);
  if (!m_cellIdDecoder)
    throw GaudiException("Unable to retrieve the cellID decoder", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
  
  m_detector = m_geoService->getDetector();
  if (!m_detector)
    throw GaudiException("Unable to retrieve the DD4hep detector from GeoSvc", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
  
  const dd4hep::rec::SurfaceManager* surfaceManager = m_detector->extension<dd4hep::rec::SurfaceManager>();
  if (!surfaceManager)
    throw GaudiException("Unable to retrieve the SurfaceManager from the DD4hep detector", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  m_surfaceMap = surfaceManager->map(m_subDetName.value());
  if (!m_surfaceMap)
    throw GaudiException("Unable to retrieve the simSurface map for subdetector " + m_subDetName.value(), "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  m_volumeManager = m_detector->volumeManager();
  if (!m_volumeManager.isValid())
    throw GaudiException("Unable to retrieve the VolumeManager from the DD4hep detector", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  
  { /* DD4hep has a transformation from the global detector coordinates to each sensors local system. The definition of the local system might change.
    * We have to make sure that a sensor always sits in the u-v plane in the local system, eg. we might have to swap the axes of the local system. 
    * We accomplish this with a transformation matrix. This is valid for all sensors in the subdetector */
    
    /* Set the sensor local transformation matrix, st. the sensors normal vector is parallel to w-axis (by simply looking at the first sensor in the map) */
    auto surfaceMapIter = m_surfaceMap->begin();
    if (surfaceMapIter == m_surfaceMap->end())
      throw GaudiException("Surface map for subdetector " + m_subDetName.value() + " is empty.", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
    const unsigned long cellID = surfaceMapIter->first;
    const dd4hep::rec::ISurface* surface = surfaceMapIter->second;
  
    TGeoHMatrix sensorTrafoMatrix = m_volumeManager.lookupDetElement(cellID).nominal().worldTransformation();
    double tempVec[3];
    
    const double epsilon = 1.0e-6; // reasonable for comparing to 1
    
    /* first: rotate sensor normal onto W-axis*/
    sensorTrafoMatrix.MasterToLocalVect(surface->normal().unit(), tempVec);
    dd4hep::rec::Vector3D n_local(tempVec[0], tempVec[1], tempVec[2]);

    if (std::abs(n_local.x() - 1.0) < epsilon) {
      debug() << "   - Local sensor normal vector is (1,0,0). Defining rotation matrix to rotate it to (0,0,1)." << endmsg;
      m_sensorNormalRotation = TGeoRotation("rot",90.,90.,0.);
    } 
    else if (std::abs(n_local.y() - 1.0) < epsilon) {
      debug() << "   - Local sensor normal vector is (0,1,0). Defining rotation matrix to rotate it to (0,0,1)." << endmsg;
      m_sensorNormalRotation = TGeoRotation("rot",0.,-90.,0.);
    } 
    else if (std::abs(n_local.z() - 1.0) < epsilon) {
      debug() << "   - Local sensor normal vector is already parallel to (0,0,1)." << endmsg;
      // no rotation needed
    } 
    else {
      error() << "Sensor local normal vector is not aligned with any local axis!" << endmsg;
    }

    /* TODO: this only makes sure that the sensor normal vector is perpendicular to the surface.
    * It does NOT make sure that the local x and y are not swapped and have the correct polarity.
    * I tried coming up with the correct transformation matrix based on the axis rotations, but tbh I am quite stumped by how the rotation acts on the vectors, and how I can fix it. */
    TGeoHMatrix M = VTXdigi_tools::ComputeSensorTrafoMatrix(cellID, m_volumeManager, m_sensorNormalRotation);
    M.MasterToLocalVect(surface->u().unit(), tempVec);
    if (std::abs(tempVec[0]-1.) > epsilon || std::abs(tempVec[1]) > epsilon || std::abs(tempVec[2]) > epsilon)
      warning() << "After rotation, local sensor U direction (" << tempVec[0] << "," << tempVec[1] << "," << tempVec[2] << ") is not parallel to (1,0,0) axis!" << endmsg;
    M.MasterToLocalVect(surface->v().unit(), tempVec);
    if (std::abs(tempVec[0]) > epsilon || std::abs(tempVec[1]-1.) > epsilon || std::abs(tempVec[2]) > epsilon)
      warning() << "After rotation, local sensor V direction (" << tempVec[0] << "," << tempVec[1] << "," << tempVec[2] << ") is not parallel to (0,1,0) axis!" << endmsg;
    M.MasterToLocalVect(surface->normal().unit(), tempVec);
    if (std::abs(tempVec[0]) > epsilon || std::abs(tempVec[1]) > epsilon || std::abs(tempVec[2]-1.) > epsilon)
      warning() << "After rotation, local sensor normal direction (" << tempVec[0] << "," << tempVec[1] << "," << tempVec[2] << ") is not parallel to (0,0,1) axis!" << endmsg;
  }

  /* IDEA / Allegro: 
   *   - subDet is child of detector, eg. Vertex 
   *   - subDetChild is child of subDet, eg. VertexBarrel
   *   - subDetChildChild are layers 
   * If this is not the case, layers might be direct children of subDet: */

  debug() << " - Retrieving subdetector " << m_subDetName.value() << " from DD4hep detector..." << endmsg;
  verbose() << "   - The detector has the following subDetectors: " << endmsg;
  for (const auto& [subDetName, subDet] : m_detector->detectors()) {
    verbose() << "   - " << subDetName << endmsg;
  }

  if (m_subDetChildName.value() != m_undefinedString) {
    /* IDEA/Allegro setup */
    const dd4hep::DetElement subDet = m_detector->detector(m_subDetName.value());
    m_subDetector = subDet.child(m_subDetChildName.value());
  }
  else {
    /* layers are direct children of subDet */
    m_subDetector = m_detector->detector(m_subDetName.value());
  }
  if (!m_subDetector)
    throw GaudiException("Unable to retrieve the subdetector DetElement " + m_subDetName.value(), "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  debug() << " - Retrieved all necessary services and dd4hep detector elements." << endmsg;
}

void VTXdigi_Modular::InitLayersAndSensors() {
  /* Get the detector & sensor geometry information from the DD4hep detector
  * 
  * Sets the members: 
  *  - m_pixelPitch 
  *  - m_sensorThickness
  *  - m_layers
  *  - */
  verbose() << " - Retrieving layers, subdetector geometry, and sensor size and pitch..." << endmsg;

  /* Find layers that we want to digitize */
  verbose() << "   - Determining relevant layers... " << endmsg;
  std::vector<int> availableLayers;
  for (const auto& [layerName, layerObj] : m_subDetector.children()) {
    dd4hep::VolumeID layerVolumeID = layerObj.volumeID();
    const int layer = m_cellIdDecoder->get(layerVolumeID, "layer");
    availableLayers.push_back(layer);
  }
  if (availableLayers.empty())
    throw GaudiException("No layers found in subdetector " + m_subDetName.value(), "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);

  if (m_layers.value().empty()) {
    /* digitize all layers in m_subDetector */
    m_layers.value() = availableLayers;
  } else {
    /* If layers-to-digitize is specified: check that all requested layers exist */

    for (const auto layer : m_layers.value()) {
      if (std::find(availableLayers.begin(), availableLayers.end(), layer) == availableLayers.end()) {
        throw GaudiException("Requested layer " + std::to_string(layer) + " to digitize, but this layer does not exist in subdetector " + m_subDetName.value(), "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
      }
    }
    debug() << "   - Digitizing only specific layers, as defined in Gaudi property by the user. All requested layers were found." << endmsg;
  }
  info() << "   - Digitizing " << m_layers.value().size() << " layers: " << m_layers.value() << " in subdetector " << m_subDetName.value() << "." << endmsg;

  /* Get pixel pitch from the segmentation (from the readout that matches our simHitCollection) 
  *  I don't know of a better way to do this (that also works for the IDEA detector model) */
  std::string simHitCollectionName;
  if (this->getProperty("SimTrackHitCollectionName", simHitCollectionName).isFailure())
    throw GaudiException("Could not retrieve SimTrackHitCollectionName property while checking geometry consistency.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);

  dd4hep::Detector::HandleMap readoutHandleMap = m_detector->readouts();
  int readoutCount = 0;
  std::string matchedReadoutKey;
  /* loop over readouts, see if one matches our simHitCollection */
  for (const auto& [readoutKey, readoutHandle] : readoutHandleMap) {
    if (simHitCollectionName.find(readoutKey) != std::string::npos) {
      ++readoutCount;
      if (readoutCount != 1)
        continue;
      matchedReadoutKey = readoutKey;

      const dd4hep::Segmentation& segmentation = m_detector->readout(readoutKey).segmentation();
      if (!segmentation.isValid())
        throw GaudiException("Segmentation for readout " + readoutKey + " is not valid while checking geometry consistency.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
      const auto cellDimensions = segmentation.cellDimensions(0); // this assumes all cells have the same dimensions (ie. only one sensor type in this readout)
      m_pixelPitch.first = cellDimensions.at(0) * 10; // convert cm to mm
      m_pixelPitch.second = cellDimensions.at(1) * 10;
      /* TODO: get pixel count from segmentation*/
    }
    else {
      verbose() << "     - Readout \"" << readoutKey << "\" does NOT MATCH SimTrackHitCollectionName \"" << simHitCollectionName << "\". skipping." << endmsg;
    }
  } // end loop over readouts
  if (readoutCount == 0)
    throw GaudiException("Could not find any readout matching SimTrackHitCollectionName " + simHitCollectionName + " in detector while checking geometry consistency.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
  else if (readoutCount > 1)
    warning() << "Found multiple (" << readoutCount << ") readouts matching SimTrackHitCollectionName \"" << simHitCollectionName << "\" in detector while checking geometry consistency. Used the first one found. Enable verbose messages for more info." << endmsg;


  /* TODO: Check that local sensor coordinates (u,v,w) are correctly defined wrt. to the global coordinates (already done in VTXdigi_Allpix2 master branch, but VERY clunky). This requires deep understanding of coordinate systems. I think there is a easy way to do this. I have not figured it out yet */


  /* Loop over all sensors in the subDetector and check that they have the same dimensions
  *  This takes << 1s for the IDEA VTX */
  int moduleNumber = 0, sensorNumber = 0;
  bool membersDefined = false;
  for (const auto& [layerKey, layerObj] : m_subDetector.children()) {
    dd4hep::VolumeID layerVolumeID = layerObj.volumeID();
    int layer = m_cellIdDecoder->get(layerVolumeID, "layer");
    /* If list of layers is given: skip layers that are not on the list*/
    if (!m_layers.value().empty()) {
      if (std::find(m_layers.value().begin(), m_layers.value().end(), layer) == m_layers.value().end()) {
        verbose() << "   - Skipping layer " << layerKey << " (layer " << layer << ", volumeID " << layerVolumeID << " with " << layerObj.children().size() << " modules) as it is not in the LayersToDigitize list." << endmsg;
        continue;
      }
    }
    verbose() << "   - Found layer \"" << layerKey << "\" (layer " << layer << ", volumeID " << layerVolumeID << ", " << layerObj.children().size() << " modules) for subDetector \"" << m_subDetName.value() << "\"." << endmsg;

    for (const auto& [moduleKey, moduleObj] : layerObj.children()) {
      ++moduleNumber;
      for (const auto& [sensorKey, sensorObj] : moduleObj.children()) {
        ++sensorNumber;
        dd4hep::VolumeID sensorVolumeID = sensorObj.volumeID();

        const auto surfaceIt = m_surfaceMap->find(sensorObj.volumeID());
        if (surfaceIt == m_surfaceMap->end()) {
          throw GaudiException("Could not find surface for sensor " + sensorKey + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layer) + " of subDetector " + m_subDetName.value() + " while checking geometry consistency.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
        }
        dd4hep::rec::ISurface* surface = surfaceIt->second;
        if (!surface) {
          throw GaudiException("Surface pointer for sensor " + sensorKey + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layer) + " of subDetector " + m_subDetName.value() + " is null while checking geometry consistency.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
        }

        const float sensorLength_u = surface->length_along_u() * 10; // convert to mm
        const float sensorLength_v = surface->length_along_v() * 10; 
        const float sensorThickness = (surface->innerThickness() + surface->outerThickness()) * 10; //

        if (!membersDefined) {
          /* Set members based on the first sensor we find */
          m_sensorThickness = sensorThickness;
          m_sensorLength.first = sensorLength_u; 
          m_sensorLength.second = sensorLength_v;
          float pixelCountU = m_sensorLength.first / m_pixelPitch.first;
          float pixelCountV = m_sensorLength.second / m_pixelPitch.second;
          if (abs(pixelCountU - std::round(pixelCountU)) > 0.001 || abs(pixelCountV - std::round(pixelCountV)) > 0.001)
            throw GaudiException("Sensor side length (" + std::to_string(m_sensorLength.first) + " x " + std::to_string(m_sensorLength.second) + ") mm and pixel pitch (" + std::to_string(m_pixelPitch.first) + " x " + std::to_string(m_pixelPitch.second) + ") mm result in a non-integer pixel count (" + std::to_string(pixelCountU) + " x " + std::to_string(pixelCountV) + ") in subDetector " + m_subDetName.value() + ".", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
          m_pixelCount.first = std::round(pixelCountU);
          m_pixelCount.second = std::round(pixelCountV);
          membersDefined = true;
          verbose() << "     - Found sensor: " << sensorKey << ", volumeID: " << sensorVolumeID << " (sensor " << sensorNumber << " in layer " << layer << "). Setting expected dimensions to (" << m_sensorLength.first << " x " << m_sensorLength.second << " x " << m_sensorThickness << ") mm3." << endmsg;
        } // if !membersDefined
        else {
          /* For all other sensors: check for consistency with first sensor */
          if (std::abs(sensorLength_u - m_sensorLength.first) > 0.001 || std::abs(sensorLength_v - m_sensorLength.second) > 0.001 || std::abs(sensorThickness - m_sensorThickness) > 0.001)
            throw GaudiException("Sensor dimension mismatch found in sensor " + sensorKey + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layer) + " of subDetector " + m_subDetName.value() + ": expected dimensions of (" + std::to_string(m_sensorLength.first) + " x " + std::to_string(m_sensorLength.second) + " x " + std::to_string(m_sensorThickness) + ") mm3, but found (" + std::to_string(sensorLength_u) + " x " + std::to_string(sensorLength_v) + " x " + std::to_string(sensorThickness) + ") mm3. This algorithm expects exactly one type of sensor per subDetector. Use different instances of the algorithm if different layers consist of different sensors.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
        } // if membersDefined
      } // loop over sensors
    } // loop over modules
  } // loop over layers
  
  info() << " - Retrieved sensor parameters: area (" << m_sensorLength.first << " x " << m_sensorLength.second << ") mm, thickness " << m_sensorThickness << " mm, pixel pitch (" << m_pixelPitch.first << " x " << m_pixelPitch.second << ") mm, pixel count (" << m_pixelCount.first << " x " << m_pixelCount.second << "). All " << sensorNumber << " sensors in the relevant layers share these parameters." << endmsg;

  if (abs(m_depletedRegionDepthCenter) > m_sensorThickness/2.f)
    warning() << "Depleted region depth center " << m_depletedRegionDepthCenter << " mm is outside the sensor (which is  " << m_sensorThickness << " mm thick). Note that DepletedRegionDepthCenter=0 lies in the middle of the sensor, ie. at half it's thickness." << endmsg;
} // InitLayersAndSensors()

void VTXdigi_Modular::InitHistograms() {
  /* Define axes globally to make adjusting them easier
  * TODO: Make some of these adjustable via Gaudi Parameters? Might not be necessary.*/
  Gaudi::Accumulators::Axis<float> axis_z{200, -200, 200};
  Gaudi::Accumulators::Axis<float> axis_cosTheta{100, 0, 1};
  Gaudi::Accumulators::Axis<float> axis_theta{4*180, 0, 180};
  Gaudi::Accumulators::Axis<float> axis_phi{4*180, -180, 180};
  
  Gaudi::Accumulators::Axis<float> axis_moduleID{2000, -0.5f, 1999.5f};
  Gaudi::Accumulators::Axis<float> axis_clusterSize{60, 0.5f, 60.5f};
  Gaudi::Accumulators::Axis<float> axis_E{1000, 0, m_sensorThickness*2000.f};
  Gaudi::Accumulators::Axis<float> axis_charge{1000, 0, m_sensorThickness*500000.f};
  Gaudi::Accumulators::Axis<float> axis_momentum_keV{10000, 0.f, 1000.f};
  Gaudi::Accumulators::Axis<float> axis_momentum_MeV{10000, 0.f, 1000.f};
  Gaudi::Accumulators::Axis<float> axis_momentum_GeV{10000, 0.f, 1000.f};
  
  Gaudi::Accumulators::Axis<float> axis_residual{1600, -200.f, 200.f};
  Gaudi::Accumulators::Axis<float> axis_residual_abs{800, 0.f, 200.f};
  Gaudi::Accumulators::Axis<float> axis_residual_pixels{2020, -50.5f, 50.5f}; // 20 in-pix bins
  Gaudi::Accumulators::Axis<float> axis_pdg{1401, -700.5f, 700.5f};

  Gaudi::Accumulators::Axis<float> axis_pixels_u{
    static_cast<unsigned int>(m_pixelCount.first),
    -0.5f,
    static_cast<float>(m_pixelCount.first+0.5)};
  Gaudi::Accumulators::Axis<float> axis_pixels_v{
    static_cast<unsigned int>(m_pixelCount.second),
    -0.5f,
    static_cast<float>(m_pixelCount.second+0.5)};

  Gaudi::Accumulators::Axis<float> axis_pathLength{500, 0.f, m_sensorThickness*1000.f*10.f};
  Gaudi::Accumulators::Axis<float> axis_pathTravel{1000, -m_sensorThickness*1000.f*10.f, m_sensorThickness*1000.f*10.f};
    
  /* Fill histograms per layer */
  for (int layer : m_layers.value()) {
    if (layer == 0) {
      axis_z = Gaudi::Accumulators::Axis<float>{100, -96.5, 96.5}; // Want to cover layer 0 of IDEA vertex det perfectly to avoid binning-edge-effects, so we use a the correct length of 185 mm
    }

    std::array< std::unique_ptr< Gaudi::Accumulators::StaticHistogram< 1, Gaudi::Accumulators::atomicity::full, float > >, hist1dArrayLen > hist1d;

    hist1d.at(hist1d_simHitE).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_energyDep",
        "SimHit deposited energy - Layer " + std::to_string(layer) + ";Energy [keV];Entries",
        axis_E
      }
    );
    hist1d.at(hist1d_simHitCharge).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_chargeDep",
        "SimHit deposited charge - Layer " + std::to_string(layer) + ";Charge [e-];Entries",
        axis_charge
      }
    );
    hist1d.at(hist1d_simHitMomentum_keV).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_momentum_keV",
        "SimHit particle momentum at hit position - Layer " + std::to_string(layer) + ";Momentum [keV/c];Entries",
        axis_momentum_keV
      }
    );
    hist1d.at(hist1d_simHitMomentum_MeV).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_momentum_MeV",
        "SimHit particle momentum at hit position - Layer " + std::to_string(layer) + ";Momentum [MeV/c];Entries",
        axis_momentum_MeV
      }
    );
    hist1d.at(hist1d_simHitMomentum_GeV).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_momentum_GeV",
        "SimHit particle momentum at hit position - Layer " + std::to_string(layer) + ";Momentum [GeV/c];Entries",
        axis_momentum_GeV
      }
    );

    hist1d.at(hist1d_simHit_z).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_global_z",
        "Global z-position of simHits - Layer " + std::to_string(layer) + ";z [mm];Entries",
        axis_z
      }
    );
    hist1d.at(hist1d_simHit_z_createdInGenerator).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_global_z_createdInGenerator",
        "Global z-position of simHits from particles created in the generator - Layer " + std::to_string(layer) + ";z [mm];Entries",
        axis_z
      }
    );
    hist1d.at(hist1d_simHit_z_createdInSimulation).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_global_z_createdInSimulation",
        "Global z-position of simHits from particles created in simulation - Layer " + std::to_string(layer) + ";z [mm];Entries",
        axis_z
      }
    );
    
    hist1d.at(hist1d_simHitPDG).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_PDG",
        "SimHit particle PDG code - Layer " + std::to_string(layer) + ";PDG;Entries",
        axis_pdg
      }
    );


    hist1d.at(hist1d_digiHitCharge).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_chargeDep",
        "DigiHit charge - Layer " + std::to_string(layer) + ";Charge [e-];Entries",
        axis_charge
      }
    );
    hist1d.at(hist1d_digiHitsPerSimHit).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_digiHits_per_simHit",
        "Raw number of digiHits created per simHit - Layer " + std::to_string(layer) + ";DigiHits per SimHit;Entries",
        axis_clusterSize // not technically correct, but works
      }
    );

    hist1d.at(hist1d_clusterSize).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize",
        "Number of pixels per cluster (after clustering algorithm) - Layer " + std::to_string(layer) + ";Cluster size [pixels];Entries",
        axis_clusterSize
      }
    );
    hist1d.at(hist1d_clusterSize_u).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_u",
        "Cluster length in u direction (after clustering algorithm) - Layer " + std::to_string(layer) + ";Cluster size [pixels];Entries",
        axis_clusterSize
      }
    );
    hist1d.at(hist1d_clusterSize_v).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_v",
        "Cluster length in v direction (after clustering algorithm) - Layer " + std::to_string(layer) + ";Cluster size [pixels];Entries",
        axis_clusterSize
      }
    );
    hist1d.at(hist1d_clusterSize_createdInGenerator).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_createdInGenerator",
        "Number of pixels per cluster (for simHits created in generator) - Layer " + std::to_string(layer) + ";Cluster size [pixels];Entries",
        axis_clusterSize
      }
    );
    hist1d.at(hist1d_clusterSize_createdInSimulation).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_createdInSimulation",
        "Number of pixels per cluster (for simHits created in simulation) - Layer " + std::to_string(layer) + ";Cluster size [pixels];Entries",
        axis_clusterSize
      }
    );

    hist1d.at(hist1d_residual_u).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_u",
        "Residual (u_digiHit - u_simHit) in local u direction - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_u_singlePixelCluster).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_u_singlePixelCluster",
        "Residual (u_digiHit - u_simHit) in local u direction for clusters with exactly one pixel - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_u_multiPixelCluster).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_u_multiPixelCluster",
        "Residual (u_digiHit - u_simHit) in local u direction for clusters with more than one pixel - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_u_createdInGenerator).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_u_createdInGenerator",
        "Residual (u_digiHit - u_simHit) in local u direction, from particles created in the generator - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_u_createdInSimulation).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_u_createdInSimulation",
        "Residual (u_digiHit - u_simHit) in local u direction, from particles created in simulation - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_v).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_v",
        "Residual (v_digiHit - v_simHit) in local v direction - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_v_singlePixelCluster).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_v_singlePixelCluster",
        "Residual (v_digiHit - v_simHit) in local v direction for clusters with exactly one pixel - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_v_multiPixelCluster).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_v_multiPixelCluster",
        "Residual (v_digiHit - v_simHit) in local v direction for clusters with more than one pixel - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_v_createdInGenerator).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_v_createdInGenerator",
        "Residual (v_digiHit - v_simHit) in local v direction, from particles created in the generator - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_v_createdInSimulation).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_v_createdInSimulation",
        "Residual (v_digiHit - v_simHit) in local v direction, from particles created in simulation - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_w).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_w",
        "Residual (w_digiHit - w_simHit) in local w direction (ie. sensor normal direction) - Layer " + std::to_string(layer) + ";Residual w [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_r).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_r",
        "Residual magnitude (|r_digiHit - r_simHit|) - Layer " + std::to_string(layer) + ";Residual r [um];Entries",
        axis_residual_abs
      }
    );

    hist1d.at(hist1d_clusterPosUncertainty_u).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterPosUncertainty_u",
        "DigiHit position uncertainty in u direction - Layer " + std::to_string(layer) + ";Position uncertainty u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_clusterPosUncertainty_v).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterPosUncertainty_v",
        "DigiHit position uncertainty in v direction - Layer " + std::to_string(layer) + ";Position uncertainty v [um];Entries",
        axis_residual
      }
    );

    hist1d.at(hist1d_pixelDistToClusterCenter_u).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_pixelDistanceToClusterCenter_u",
        "Distance of pixel to center of cluster in u direction - Layer " + std::to_string(layer) + ";Distance to cluster center u [pix];Entries",
        axis_residual_pixels
      }
    );
    hist1d.at(hist1d_pixelDistToClusterCenter_v).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_pixelDistanceToClusterCenter_v",
        "Distance of pixel to center of cluster in v direction - Layer " + std::to_string(layer) + ";Distance to cluster center v [pix];Entries",
        axis_residual_pixels
      }
    );

    hist1d.at(hist1d_pathTravel_u).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_u",
        "Path travel length inside the sensor active volume in u direction - Layer " + std::to_string(layer) + ";Path length in u [um];Entries",
        axis_pathTravel
      }
    );
    hist1d.at(hist1d_pathTravel_v).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_v",
        "Path travel length inside the sensor active volume in v direction - Layer " + std::to_string(layer) + ";Path length in v [um];Entries",
        axis_pathTravel
      }
    );
    hist1d.at(hist1d_pathTravel_r).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_r",
        "Path travel length inside the sensor active volume - Layer " + std::to_string(layer) + ";Path length [um];Entries",
        axis_pathLength
      }
    ); 

    m_hist1d.emplace(layer, std::move(hist1d));

    /* -- 1d-profile-hist -- */

    std::array< std::unique_ptr< Gaudi::Accumulators::StaticProfileHistogram<1,Gaudi::Accumulators::atomicity::full,float>>, histProfile1dArrayLen> histProfile1d;

    histProfile1d.at(histProfile1d_digiHitCharge_vs_global_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_chargeDep_vs_global_z",
        "DigiHit charge - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Charge [e-]",
        axis_z
      }
    );

    histProfile1d.at(histProfile1d_clusterSize_vs_global_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_vs_global_z",
        "Cluster size - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Pixels per cluster",
        axis_z
      }
    );
    histProfile1d.at(histProfile1d_clusterSize_u_vs_global_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_u_vs_global_z",
        "Cluster length along u - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Cluster length in u [pix]",
        axis_z
      }
    );
    histProfile1d.at(histProfile1d_clusterSize_v_vs_global_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_v_vs_global_z",
        "Cluster length along v - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Cluster length in v [pix]",
        axis_z
      }
    );

    histProfile1d.at(histProfile1d_residual_u_vs_global_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_u_vs_global_z",
        "Residual u (u_digiHit - u_simHit) vs. global z position - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Residual u [um]",
        axis_z
      }
    );
    histProfile1d.at(histProfile1d_residual_v_vs_global_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_v_vs_global_z",
        "Residual v (v_digiHit - v_simHit) vs. global z position - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Residual v [um]",
        axis_z
      }
    );
    histProfile1d.at(histProfile1d_residual_r_vs_global_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_r_vs_global_z",
        "Residual r (|r_digiHit - r_simHit|) vs. global z position - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Residual r [um]",
        axis_z
      }
    );

    histProfile1d.at(histProfile1d_pathTravel_u_vs_global_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_u_vs_global_z",
        "Path travel length inside the sensor active volume in u direction vs. global z position - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Path length in u [um]",
        axis_z
      }
    );
    histProfile1d.at(histProfile1d_pathTravel_v_vs_global_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_v_vs_global_z",
        "Path travel length inside the sensor active volume in v direction vs. global z position - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Path length in v [um]",
        axis_z
      }
    );
    histProfile1d.at(histProfile1d_pathTravel_r_vs_global_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_r_vs_global_z",
        "Path travel length inside the sensor active volume vs. global z position - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Path length [um]",
        axis_z
      }
    );

    m_histProfile1d.emplace(layer, std::move(histProfile1d));

    /* -- 2d-hist -- */

    std::array< std::unique_ptr< Gaudi::Accumulators::StaticHistogram< 2, Gaudi::Accumulators::atomicity::full, float > >, hist2dArrayLen>  hist2d;

    hist2d.at(hist2d_digiHitCharge_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_chargeDep_vs_global_z_2D",
        "DigiHit charge - Layer " + std::to_string(layer) + ";Global z position [mm];Charge [e-]",
        axis_z,
        axis_charge
      }
    );

    hist2d.at(hist2d_hitMap_simHits).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_hitMap",
        "SimHit hit map (pixel in which this simHit position lies) - Layer " + std::to_string(layer) + ";Pixel u;Pixel v;Entries",
        axis_pixels_u,
        axis_pixels_v
      }
    );
    hist2d.at(hist2d_hitMap_simHits_createdInGenerator).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_hitMap_createdInGenerator",
        "SimHit hit map (pixel in which this simHit position lies) for particles created in the generator - Layer " + std::to_string(layer) + ";Pixel u;Pixel v;Entries",
        axis_pixels_u,
        axis_pixels_v
      }
    );
    hist2d.at(hist2d_hitMap_simHits_createdInSimulation).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_hitMap_createdInSimulation",
        "SimHit hit map (pixel in which this simHit position lies) for particles created in the simulation - Layer " + std::to_string(layer) + ";Pixel u;Pixel v;Entries",
        axis_pixels_u,
        axis_pixels_v
      }
    );
    hist2d.at(hist2d_hitMap_pixelHits).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_hitMap",
        "Hit map (pixels that collect charge over the threshold) - Layer " + std::to_string(layer) + ";Pixel u;Pixel v;Entries",
        axis_pixels_u,
        axis_pixels_v
      }
    );
    hist2d.at(hist2d_clusterSize_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_vs_global_z_2D",
        "Cluster size vs. global hit z - Layer " + std::to_string(layer) + ";Global z position [mm];Pixels per cluster;Entries",
        axis_z,
        axis_clusterSize
      }
    );
    hist2d.at(hist2d_clusterSize_u_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_u_vs_global_z_2D",
        "Cluster length along u vs. global hit z - Layer " + std::to_string(layer) + ";Global z position [mm];Cluster length in u [pix];Entries",
        axis_z,
        axis_clusterSize
      }
    );
    hist2d.at(hist2d_clusterSize_v_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_v_vs_global_z_2D",
        "Cluster length along v vs. global hit z - Layer " + std::to_string(layer) + ";Global z position [mm];Cluster length in v [pix];Entries",
        axis_z,
        axis_clusterSize
      }
    );
    hist2d.at(hist2d_clusterSize_vs_global_z_createdInGenerator).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_vs_global_z_createdInGenerator_2D",
        "Cluster size vs. global hit z (simHit particle created in generator) - Layer " + std::to_string(layer) + ";Global z position [mm];Pixels per cluster;Entries",
        axis_z,
        axis_clusterSize
      }
    );
    hist2d.at(hist2d_clusterSize_vs_global_z_createdInSim).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_vs_global_z_createdInSimulation_2D",
        "Cluster size vs. global hit z (simHit particle created in simulation) - Layer " + std::to_string(layer) + ";Global z position [mm];Pixels per cluster;Entries",
        axis_z,
        axis_clusterSize
      }
    );

    hist2d.at(hist2d_residual_u_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_u_vs_global_z_2D",
        "Residual u (u_digiHit - u_simHit) vs. global z position - Layer " + std::to_string(layer) + ";Global z position [mm];Residual u [um];Entries",
        axis_z,
        axis_residual
      }
    );
    hist2d.at(hist2d_residual_v_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_v_vs_global_z_2D",
        "Residual v (v_digiHit - v_simHit) vs. global z position - Layer "+ std::to_string(layer) + ";Global z position [mm];Residual v [um];Entries",
        axis_z,
        axis_residual
      }
    );
    hist2d.at(hist2d_residual_r_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_r_vs_global_z_2D",
        "Residual r (|r_digiHit - r_simHit|) vs. global z position - Layer " + std::to_string(layer) + ";Global z position [mm];Residual r [um];Entries",
        axis_z,
        axis_residual_abs
      }
    );

    hist2d.at(hist2d_residual_vs_clusterPosUncertainty_u).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_vs_clusterPosUncertainty_u_2D",
        "Residual r (|r_digiHit - r_simHit|) vs. cluster position uncertainty in u direction - Layer " + std::to_string(layer) + ";Cluster position uncertainty u [um];Residual u [um];Entries",
        axis_residual_abs,
        axis_residual_abs
      }
    );
    hist2d.at(hist2d_residual_vs_clusterPosUncertainty_v).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_vs_clusterPosUncertainty_v_2D",
        "Residual r (|r_digiHit - r_simHit|) vs. cluster position uncertainty in v direction - Layer " + std::to_string(layer) + ";Cluster position uncertainty v [um];Residual v [um];Entries",
        axis_residual_abs,
        axis_residual_abs
      }
    );

    hist2d.at(hist2d_pathTravel_u_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_u_vs_global_z_2D",
        "Path travel length inside the sensor active volume in u direction vs. global z position - Layer " + std::to_string(layer) + ";Global z position [mm];Path length in u [um];Entries",
        axis_z,
        axis_pathTravel
      }
    );
    hist2d.at(hist2d_pathTravel_u_vs_global_z_createdInGenerator).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_u_vs_global_z_createdInGenerator_2D",
        "Path travel length inside the sensor active volume in u direction vs. global z position, for particles created in the generator - Layer " + std::to_string(layer) + ";Global z position [mm];Path length in u [um];Entries",
        axis_z,
        axis_pathTravel
      }
    );
    hist2d.at(hist2d_pathTravel_u_vs_global_z_createdInSimulation).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_u_vs_global_z_createdInSimulation_2D",
        "Path travel length inside the sensor active volume in u direction vs. global z position, for particles created in the simulation - Layer " + std::to_string(layer) + ";Global z position [mm];Path length in u [um];Entries",
        axis_z,
        axis_pathTravel
      }
    );
    hist2d.at(hist2d_pathTravel_v_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_v_vs_global_z_2D",
        "Path travel length inside the sensor active volume in v direction vs. global z position - Layer " + std::to_string(layer) + ";Global z position [mm];Path length in v [um];Entries",
        axis_z,
        axis_pathTravel
      }
    );
    hist2d.at(hist2d_pathTravel_v_vs_global_z_createdInGenerator).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_v_vs_global_z_createdInGenerator_2D",
        "Path travel length inside the sensor active volume in v direction vs. global z position, for particles created in the generator - Layer " + std::to_string(layer) + ";Global z position [mm];Path length in v [um];Entries",
        axis_z,
        axis_pathTravel
      }
    );
    hist2d.at(hist2d_pathTravel_v_vs_global_z_createdInSimulation).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_v_vs_global_z_createdInSimulation_2D",
        "Path travel length inside the sensor active volume in v direction vs. global z position, for particles created in the simulation - Layer " + std::to_string(layer) + ";Global z position [mm];Path length in v [um];Entries",
        axis_z,
        axis_pathTravel
      }
    );
    hist2d.at(hist2d_pathTravel_w_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_w_vs_global_z_2D",
        "Path travel length inside the sensor active volume in w direction vs. global z position - Layer " + std::to_string(layer) + ";Global z position [mm];Path length in w [um];Entries",
        axis_z,
        axis_pathTravel
      }
    );
    hist2d.at(hist2d_pathTravel_w_vs_global_z_createdInGenerator).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_w_vs_global_z_createdInGenerator_2D",
        "Path travel length inside the sensor active volume in w direction vs. global z position, for particles created in the generator - Layer " + std::to_string(layer) + ";Global z position [mm];Path length in w [um];Entries",
        axis_z,
        axis_pathTravel
      }
    );
    hist2d.at(hist2d_pathTravel_w_vs_global_z_createdInSimulation).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_w_vs_global_z_createdInSimulation_2D",
        "Path travel length inside the sensor active volume in w direction vs. global z position, for particles created in the simulation - Layer " + std::to_string(layer) + ";Global z position [mm];Path length in w [um];Entries",
        axis_z,
        axis_pathTravel
      }
    );
    hist2d.at(hist2d_pathTravel_r_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_r_vs_global_z_2D",
        "Path travel length inside the sensor active volume vs. global z position - Layer " + std::to_string(layer) + ";Global z position [mm];Path length [um];Entries",
        axis_z,
        axis_pathLength
      }
    );

    hist2d.at(hist2d_simHit_xy).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_xy_2D",
        "SimHit x vs y position (global) - Layer " + std::to_string(layer) + ";X position [mm];Y position [mm];Entries",
        axis_z,
        axis_z
      }
    );
    hist2d.at(hist2d_simHit_xz).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_xz_2D",
        "SimHit x vs z position (global) - Layer " + std::to_string(layer) + ";X position [mm];Z position [mm];Entries",
        axis_z,
        axis_z
      }
    );
    hist2d.at(hist2d_simHit_yz).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_yz_2D",
        "SimHit y vs z position (global) - Layer " + std::to_string(layer) + ";Y position [mm];Z position [mm];Entries",
        axis_z,
        axis_z
      }
    );

    m_hist2d.emplace(layer, std::move(hist2d));

  } /* loop over layers */

  /* global (eg covering all layers) histograms */

  m_hist1dglobal.at(hist1dglobal_pathTravel_r).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Global/simHit_pathLength",
        "Path travel length in sensor active volume, as computed by VTXdigi_tools reconstruction;Path length [um];Entries",
        axis_pathLength
      }
    );
  m_hist1dglobal.at(hist1dglobal_pathTravel_r_Geant4).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
      "Global/simHit_pathLength_Geant4",
      "Path length in sensor active volume, as given by Geant4;Path length [um];Entries",
      axis_pathLength
    }
  );
  m_hist1dglobal.at(hist1dglobal_pathTravel_r_ratio).reset(
    new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
      "Global/simHit_pathLength_ratio",
      "Path length in sensor active volume divided by path length given by Geant4;Path length / path length Geant4;Entries",
      {500, 0.f, 2.f}
    }
  );

}

/* ---- Eventloop functions ---- */

bool VTXdigi_Modular::CheckEventSetup(const edm4hep::SimTrackerHitCollection& simTrackerHits, const edm4hep::EventHeaderCollection& headers) const {
  if (m_counter_eventsRead.value() % m_infoPrintInterval.value() == 0)
  info() << "PROCESSING event [run " << headers.at(0).getRunNumber() << ", event " << headers.at(0).getEventNumber() << ", found " << simTrackerHits.size() << " simHits]. " << m_counter_eventsRead.value() << " events so far." << endmsg;
  /* events are not necessarily numbered sequentially... */
  ++m_counter_eventsRead;

  if (simTrackerHits.size()==0) {
    debug() << " - No SimTrackerHits in collection, returning empty output collections" << endmsg;
    ++m_counter_eventsRejected_noSimHits;
    return false;
  }

  ++m_counter_eventsAccepted;
  return true;
}

bool VTXdigi_Modular::CheckSimhitLayer(const edm4hep::SimTrackerHit& simHit) const {
  const int layer = m_cellIdDecoder->get(simHit.getCellID(), "layer");
  if (m_layers.value().size()>0) { 
    if (std::find(m_layers.value().begin(), m_layers.value().end(),  layer) == m_layers.value().end()) {
      // verbose() << " - Dismissing simHit in layer " << layer << ". (not in the list of layers to digitize)" << endmsg;
      ++m_counter_simHitsRejected_LayerNotToBeDigitized;
      return false;
    }
  }
  else {
    // verbose() << " - All layers are to be digitized, as property \"LayersToDigitize\" is not set ." << endmsg;
  }

  ++m_counter_simHitsAccepted;
  return true;
}

std::vector<VTXdigi_tools::Cluster> VTXdigi_Modular::Clusterize(const VTXdigi_tools::HitMap& hitMap) const {
  if (!hitMap.Hits().size()) {
    debug() << "     - No pixels with charge found on this sensor." << endmsg;
    return {};
  }

  if (m_clusterize.value()) {
    debug() << "     - Clusterizing " << hitMap.Hits().size() << " hits with a total charge of " << hitMap.GetTotalCharge() << " e." << endmsg;
    return hitMap.ComputeClusters();
  }
  else {
    return hitMap.ComputeClusters_singePixels();
  }
}

void VTXdigi_Modular::CreateDigiHits(edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitLinks, const dd4hep::DDSegmentation::CellID& cellID, const TGeoHMatrix& trafoMatrix, const std::vector<VTXdigi_tools::Cluster>& clusters) const {

  // Called for each cluster on a sensor, so cellID and direction vectors are same among these clusters
  const dd4hep::rec::Vector3D direction_u_3d = VTXdigi_tools::LocalToGlobal(dd4hep::rec::Vector3D(1, 0, 0), trafoMatrix);
  const edm4hep::Vector2f direction_u = edm4hep::Vector2f(direction_u_3d.theta(), direction_u_3d.phi()); // edm4hep expects direction in to given as (theta, phi)

  const dd4hep::rec::Vector3D direction_v_3d = VTXdigi_tools::LocalToGlobal(dd4hep::rec::Vector3D(0, 1, 0), trafoMatrix);
  const edm4hep::Vector2f direction_v = edm4hep::Vector2f(direction_v_3d.theta(), direction_v_3d.phi()); // edm4hep expects direction in to given as (theta, phi)

  for (auto& cluster : clusters) {
    if (cluster.simHits.empty()) {
      error() << "Cluster with no contributing simHits found." << endmsg;
    }
    if (cluster.charge <= 0) {
      error() << "Cluster with non-positive charge found." << endmsg;
    }

    edm4hep::MutableTrackerHitPlane digiHit = digiHits.create();
    ++m_counter_digiHitsCreated;
    /* TODO: do we want cellID with segmentation? (instead of using the short cellID that only encodes the sensor and not the pixel) */
    digiHit.setCellID(cellID);
    digiHit.setEDep(cluster.charge / VTXdigi_tools::kChargePerkeV);
    

    // collect position
    const std::pair<float, float> clusterPos_index = cluster.ComputePos();
    const dd4hep::rec::Vector3D clusterPos_local = VTXdigi_tools::ComputePosFromPixIndex_local(clusterPos_index, m_sensorLength, m_pixelPitch, m_depletedRegionDepthCenter);
    const dd4hep::rec::Vector3D clusterPos_global = VTXdigi_tools::LocalToGlobal(clusterPos_local, trafoMatrix);
    debug() << "     - Found cluster with " << cluster.pixels.size() << " pixels, charge " << cluster.charge << ", center at (" << clusterPos_index.first << ", " << clusterPos_index.second << "). Has " << cluster.simHits.size() << " contributing simHits." << endmsg;
    digiHit.setPosition(VTXdigi_tools::ConvertVector(clusterPos_global));
    
    // pos uncertainty
    digiHit.setU(direction_u);
    digiHit.setV(direction_v);
    if (m_positionUncertainty.value().empty()) {
      std::pair<float, float> clusterPos_unc = cluster.ComputePosUncertainty_ChargeWeighted(clusterPos_index);
      digiHit.setDu(clusterPos_unc.first * m_pixelPitch.first);
      digiHit.setDv(clusterPos_unc.second * m_pixelPitch.second);
    }
    else if (m_positionUncertainty.value().size() == 2) {
      digiHit.setDu(m_positionUncertainty.value().at(0));
      digiHit.setDv(m_positionUncertainty.value().at(1));
    }
    else if (m_positionUncertainty.value().size() == 6) {
      int clusterSize_u = cluster.GetSize(0);
      if (clusterSize_u == 1)
        digiHit.setDu(m_positionUncertainty.value().at(0));
      else if (clusterSize_u == 2)
        digiHit.setDu(m_positionUncertainty.value().at(1));
      else
        digiHit.setDu(m_positionUncertainty.value().at(2));

      int clusterSize_v = cluster.GetSize(1);
      if (clusterSize_v == 1)
        digiHit.setDv(m_positionUncertainty.value().at(3));
      else if (clusterSize_v == 2)
        digiHit.setDv(m_positionUncertainty.value().at(4));
      else
        digiHit.setDv(m_positionUncertainty.value().at(5));
    }
    debug() << "         - Set digiHit position uncertainty to (" << digiHit.getDu() << ", " << digiHit.getDv() << ") mm." << endmsg;

    // collect timestamp & create links
    float timeStamp = std::numeric_limits<float>::max(); // TODO: use earliest timestamp among simHits
    for (const auto& simHit : cluster.simHits) {
      auto link = digiHitLinks.create();
      link.setFrom(digiHit);
      link.setTo(*(simHit->hitPtr()));

      timeStamp = std::min(timeStamp, simHit->hitPtr()->getTime());
    }
    if (m_smearing_time > 0.f)
      timeStamp += m_rndm_time;
    digiHit.setTime(timeStamp); 
    

    /* TODO: estimate uncertainty from pitch/sqrt(12), or from cluster size in u/v*/

    if (m_debugHistograms.value()) {
      FillHistograms_perDigiHit(cluster.simHits, digiHit, trafoMatrix, cluster.pixels.size(), cluster.GetSize(0), cluster.GetSize(1));

      for (const VTXdigi_tools::Pixel* pix : cluster.pixels)
        FillHistograms_perPixel(cellID, *pix, clusterPos_index);
    }
  } /* loop over clusters */
}

/* ---- Histogramming functions ---- */

void VTXdigi_Modular::FillHistograms_perSimHit(const VTXdigi_tools::SimHitWrapper& simHit) const {
  /* executed once for each simHit, no cuts applied before */

  const int layer = simHit.layer();
  const dd4hep::rec::Vector3D simHitMomentum = VTXdigi_tools::ConvertVector(simHit.hitPtr()->getMomentum()); // in GeV
  const dd4hep::rec::Vector3D simHitPos_global = VTXdigi_tools::ConvertVector(simHit.hitPtr()->getPosition());
  const TGeoHMatrix trafoMatrix = VTXdigi_tools::ComputeSensorTrafoMatrix(simHit.cellID(), m_volumeManager, m_sensorNormalRotation);
  const dd4hep::rec::Vector3D simHitPos_local = VTXdigi_tools::GlobalToLocal(simHitPos_global, trafoMatrix);
  const std::pair<int, int> i_uv = VTXdigi_tools::ComputePixelIndices(simHitPos_local, m_pixelPitch, m_pixelCount);

  ++(*m_hist1d.at(layer).at(hist1d_simHitE))[simHit.hitPtr()->getEDep() * (dd4hep::GeV / dd4hep::keV)]; 
  ++(*m_hist1d.at(layer).at(hist1d_simHitCharge))[simHit.charge()]; 
  ++(*m_hist1d.at(layer).at(hist1d_simHitMomentum_keV))[simHitMomentum.r() * 1.E6]; 
  ++(*m_hist1d.at(layer).at(hist1d_simHitMomentum_MeV))[simHitMomentum.r() * 1.E3]; 
  ++(*m_hist1d.at(layer).at(hist1d_simHitMomentum_GeV))[simHitMomentum.r()]; 

  ++(*m_hist1d.at(layer).at(hist1d_simHit_z))[simHitPos_global.z()];

  if ( simHit.isCreatedInGenerator() ) {
    ++(*m_hist1d.at(layer).at(hist1d_simHit_z_createdInGenerator))[simHitPos_global.z()];
    ++(*m_hist2d.at(layer).at(hist2d_hitMap_simHits_createdInGenerator))[{i_uv.first, i_uv.second}];
  }
  else {
    ++(*m_hist1d.at(layer).at(hist1d_simHit_z_createdInSimulation))[simHitPos_global.z()];
    ++(*m_hist2d.at(layer).at(hist2d_hitMap_simHits_createdInSimulation))[{i_uv.first, i_uv.second}];
  }

  ++(*m_hist1d.at(layer).at(hist1d_simHitPDG))[simHit.hitPtr()->getParticle().getPDG()];
  
  ++(*m_hist2d.at(layer).at(hist2d_hitMap_simHits))[{i_uv.first, i_uv.second}];
  
  ++(*m_hist2d.at(layer).at(hist2d_simHit_xy))[{simHitPos_global.x(), simHitPos_global.y()}];
  ++(*m_hist2d.at(layer).at(hist2d_simHit_xz))[{simHitPos_global.x(), simHitPos_global.z()}];
  ++(*m_hist2d.at(layer).at(hist2d_simHit_yz))[{simHitPos_global.y(), simHitPos_global.z()}];
}

void VTXdigi_Modular::FillHistograms_perPixel(const dd4hep::DDSegmentation::CellID& cellID, const VTXdigi_tools::Pixel& pix, const std::pair<float, float> clusterPos_local) const {
  /* executed once for each pixel */

  const int layer = VTXdigi_tools::GetLayer(cellID, m_cellIdDecoder);
  const std::pair<int, int> i_uv = pix.index;

  ++(*m_hist2d.at(layer).at(hist2d_hitMap_pixelHits))[{i_uv.first, i_uv.second}];

  if (clusterPos_local.first != 0 && clusterPos_local.second != 0) { 
    /* cluster at (0,0) means that hits weren't clustered */
    const float distToClusterCenter_u = static_cast<float>(pix.index.first) - clusterPos_local.first; // in pix
    const float distToClusterCenter_v = static_cast<float>(pix.index.second) - clusterPos_local.second;

    ++(*m_hist1d.at(layer).at(hist1d_pixelDistToClusterCenter_u))[distToClusterCenter_u]; // convert to um
    ++(*m_hist1d.at(layer).at(hist1d_pixelDistToClusterCenter_v))[distToClusterCenter_v];
  }
}

void VTXdigi_Modular::FillHistograms_perDigiHit(const std::unordered_set<const VTXdigi_tools::SimHitWrapper*>& simHits, const edm4hep::TrackerHitPlane& digiHit, const TGeoHMatrix& trafoMatrix, const int clusterSize, const int clusterSize_u, const int clusterSize_v) const {
  /* executed once for each digiHit */
  const int layer = VTXdigi_tools::GetLayer(digiHit.getCellID(), m_cellIdDecoder);
  const dd4hep::rec::Vector3D pos_global = VTXdigi_tools::ConvertVector(digiHit.getPosition());
  const dd4hep::rec::Vector3D pos_local = VTXdigi_tools::GlobalToLocal(pos_global, trafoMatrix);
  
  ++(*m_hist1d.at(layer).at(hist1d_digiHitCharge))[digiHit.getEDep() * VTXdigi_tools::kChargePerkeV]; 

  ++(*m_hist1d.at(layer).at(hist1d_clusterSize))[clusterSize];
  ++(*m_hist1d.at(layer).at(hist1d_clusterSize_u))[clusterSize_u];
  ++(*m_hist1d.at(layer).at(hist1d_clusterSize_v))[clusterSize_v];
  
  (*m_histProfile1d.at(layer).at(histProfile1d_clusterSize_vs_global_z))[pos_global.z()] += clusterSize;
  (*m_histProfile1d.at(layer).at(histProfile1d_clusterSize_u_vs_global_z))[pos_global.z()] += clusterSize_u;
  (*m_histProfile1d.at(layer).at(histProfile1d_clusterSize_v_vs_global_z))[pos_global.z()] += clusterSize_v;
  
  ++(*m_hist2d.at(layer).at(hist2d_clusterSize_vs_global_z))[{pos_global.z(), clusterSize}];
  ++(*m_hist2d.at(layer).at(hist2d_clusterSize_u_vs_global_z))[{pos_global.z(), clusterSize_u}];
  ++(*m_hist2d.at(layer).at(hist2d_clusterSize_v_vs_global_z))[{pos_global.z(), clusterSize_v}];

  ++(*m_hist1d.at(layer).at(hist1d_clusterPosUncertainty_u))[ digiHit.getDu() * 1000.f ]; // convert to um
  ++(*m_hist1d.at(layer).at(hist1d_clusterPosUncertainty_v))[ digiHit.getDv() * 1000.f ];

  /* per simhit that is connected to this digiHit (cluster) */
  for (const auto& simHit : simHits) {
    const dd4hep::rec::Vector3D simHitPos_local = simHit->truthPos();
    const dd4hep::rec::Vector3D simHitPos_global = VTXdigi_tools::LocalToGlobal(simHitPos_local, trafoMatrix);
    // const dd4hep::rec::Vector3D simHitPos_global = VTXdigi_tools::ConvertVector(simHit->hitPtr()->getPosition());
    // const dd4hep::rec::Vector3D simHitPos_local = VTXdigi_tools::GlobalToLocal(simHitPos_global, trafoMatrix); 

    const dd4hep::rec::Vector3D residual_local = pos_local - simHitPos_local; // residual = observed - predicted

    const float hit_z = simHitPos_global.z();
    
    (*m_histProfile1d.at(layer).at(histProfile1d_residual_u_vs_global_z))[hit_z] += residual_local.x()*1000.f; // convert from mm to um
    (*m_histProfile1d.at(layer).at(histProfile1d_residual_v_vs_global_z))[hit_z] += residual_local.y()*1000.f;
    (*m_histProfile1d.at(layer).at(histProfile1d_residual_r_vs_global_z))[hit_z] += residual_local.r()*1000.f;

    ++(*m_hist2d.at(layer).at(hist2d_residual_u_vs_global_z))[{hit_z, residual_local.x()*1000.f}];
    ++(*m_hist2d.at(layer).at(hist2d_residual_v_vs_global_z))[{hit_z, residual_local.y()*1000.f}];
    ++(*m_hist2d.at(layer).at(hist2d_residual_r_vs_global_z))[{hit_z, residual_local.r()*1000.f}];

    ++(*m_hist1d.at(layer).at(hist1d_residual_u))[ residual_local.x()*1000.f ];
    ++(*m_hist1d.at(layer).at(hist1d_residual_v))[ residual_local.y()*1000.f ];
    ++(*m_hist1d.at(layer).at(hist1d_residual_w))[ residual_local.z()*1000.f ];
    ++(*m_hist1d.at(layer).at(hist1d_residual_r))[ residual_local.r()*1000.f ];
    
    ++(*m_hist2d.at(layer).at(hist2d_residual_vs_clusterPosUncertainty_u))[{std::abs(digiHit.getDu()) * 1000.f, std::abs(residual_local.x())*1000.f}]; // convert to um
    ++(*m_hist2d.at(layer).at(hist2d_residual_vs_clusterPosUncertainty_v))[{std::abs(digiHit.getDv()) * 1000.f, std::abs(residual_local.y())*1000.f}];

    if (clusterSize == 1) {
      ++(*m_hist1d.at(layer).at(hist1d_residual_u_singlePixelCluster))[ residual_local.x()*1000.f ];
      ++(*m_hist1d.at(layer).at(hist1d_residual_v_singlePixelCluster))[ residual_local.y()*1000.f ];
    }
    else if (clusterSize > 1) {
      ++(*m_hist1d.at(layer).at(hist1d_residual_u_multiPixelCluster))[ residual_local.x()*1000.f ];
      ++(*m_hist1d.at(layer).at(hist1d_residual_v_multiPixelCluster))[ residual_local.y()*1000.f ];
    }

    if ( simHit->isCreatedInGenerator() ) {
      ++(*m_hist1d.at(layer).at(hist1d_clusterSize_createdInGenerator))[clusterSize];
      ++(*m_hist2d.at(layer).at(hist2d_clusterSize_vs_global_z_createdInGenerator))[{hit_z, clusterSize}];

      ++(*m_hist1d.at(layer).at(hist1d_residual_u_createdInGenerator))[ residual_local.x()*1000.f ];
      ++(*m_hist1d.at(layer).at(hist1d_residual_v_createdInGenerator))[ residual_local.y()*1000.f ];
    }
    else {
      ++(*m_hist1d.at(layer).at(hist1d_clusterSize_createdInSimulation))[clusterSize];
      ++(*m_hist2d.at(layer).at(hist2d_clusterSize_vs_global_z_createdInSim))[{hit_z, clusterSize}];

      ++(*m_hist1d.at(layer).at(hist1d_residual_u_createdInSimulation))[ residual_local.x()*1000.f ];
      ++(*m_hist1d.at(layer).at(hist1d_residual_v_createdInSimulation))[ residual_local.y()*1000.f ];
    }
  }
}

void VTXdigi_Modular::FillHistograms_fromChargeCollector_perSimHit(const int layer, const dd4hep::rec::Vector3D& pathTravel, const float pathLength_Geant4, const dd4hep::rec::Vector3D& truthPos_local, const TGeoHMatrix& trafoMatrix, const bool createdInGenerator) const {
  if (!m_debugHistograms.value()) return;

  const float pathLength = pathTravel.r(); // in mm
  const float factor_um_per_mm = 1000.f;
  const dd4hep::rec::Vector3D truthPos_global = VTXdigi_tools::LocalToGlobal(truthPos_local, trafoMatrix);
  
  ++(*m_hist1dglobal.at(hist1dglobal_pathTravel_r))[pathLength*factor_um_per_mm]; // convert from mm to um
  ++(*m_hist1dglobal.at(hist1dglobal_pathTravel_r_Geant4))[pathLength_Geant4*factor_um_per_mm];
  if (pathLength_Geant4 != 0.f) {
    ++(*m_hist1dglobal.at(hist1dglobal_pathTravel_r_ratio))[pathLength / pathLength_Geant4];
  }
    
  ++(*m_hist1d.at(layer).at(hist1d_pathTravel_u))[pathTravel.x()*factor_um_per_mm];
  ++(*m_hist1d.at(layer).at(hist1d_pathTravel_v))[pathTravel.y()*factor_um_per_mm];
  ++(*m_hist1d.at(layer).at(hist1d_pathTravel_r))[pathLength*factor_um_per_mm];
    
  (*m_histProfile1d.at(layer).at(histProfile1d_pathTravel_u_vs_global_z))[truthPos_global.z()] += pathTravel.x()*factor_um_per_mm;
  (*m_histProfile1d.at(layer).at(histProfile1d_pathTravel_v_vs_global_z))[truthPos_global.z()] += pathTravel.y()*factor_um_per_mm;
  (*m_histProfile1d.at(layer).at(histProfile1d_pathTravel_r_vs_global_z))[truthPos_global.z()] += pathLength*factor_um_per_mm;

  ++(*m_hist2d.at(layer).at(hist2d_pathTravel_u_vs_global_z))[{truthPos_global.z(), pathTravel.x()*factor_um_per_mm}];
  ++(*m_hist2d.at(layer).at(hist2d_pathTravel_v_vs_global_z))[{truthPos_global.z(), pathTravel.y()*factor_um_per_mm}];
  ++(*m_hist2d.at(layer).at(hist2d_pathTravel_w_vs_global_z))[{truthPos_global.z(), pathTravel.z()*factor_um_per_mm}];
  ++(*m_hist2d.at(layer).at(hist2d_pathTravel_r_vs_global_z))[{truthPos_global.z(), pathLength*factor_um_per_mm}];

  if (createdInGenerator) {
    ++(*m_hist2d.at(layer).at(hist2d_pathTravel_u_vs_global_z_createdInGenerator))[{truthPos_global.z(), pathTravel.x()*factor_um_per_mm}];
    ++(*m_hist2d.at(layer).at(hist2d_pathTravel_v_vs_global_z_createdInGenerator))[{truthPos_global.z(), pathTravel.y()*factor_um_per_mm}];
    ++(*m_hist2d.at(layer).at(hist2d_pathTravel_w_vs_global_z_createdInGenerator))[{truthPos_global.z(), pathTravel.z()*factor_um_per_mm}];
  }
  else {
    ++(*m_hist2d.at(layer).at(hist2d_pathTravel_u_vs_global_z_createdInSimulation))[{truthPos_global.z(), pathTravel.x()*factor_um_per_mm}];
    ++(*m_hist2d.at(layer).at(hist2d_pathTravel_v_vs_global_z_createdInSimulation))[{truthPos_global.z(), pathTravel.y()*factor_um_per_mm}];
    ++(*m_hist2d.at(layer).at(hist2d_pathTravel_w_vs_global_z_createdInSimulation))[{truthPos_global.z(), pathTravel.z()*factor_um_per_mm}];
  }
}

void VTXdigi_Modular::PrintCountersSummary() const {
  const int colWidths[] = {65, 10};  
  info() << " Counters summary: " << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Events read"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_eventsRead.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Events rejected (no simHits)"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_eventsRejected_noSimHits.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "Events accepted"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_eventsAccepted.value() << " |" << endmsg;

  info() << " | " << std::setw(colWidths[0]) << std::left << "SimTrackerHits read"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_simHitsRead.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "SimTrackerHits rejected (layer ignored)"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_simHitsRejected_LayerNotToBeDigitized.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "SimTrackerHits accepted"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_simHitsAccepted.value() << " |" << endmsg;

  info()	<<	" | "	<<	std::setw(colWidths[0])	<<	std::left	<<	"Digi hits created"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_digiHitsCreated.value() << " |" << endmsg;
}

