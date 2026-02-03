#include "VTXdigi_Modular.h"

DECLARE_COMPONENT(VTXdigi_Modular)

/* Notes ~ Jona 2025-09
 * - all lengths in mm (edm4hep already does this)
 *   BUT: dd4hep uses cm internally, so convert when passing values to/from dd4hep via dd4hep::mm = 0.1
 *        mm -> cm: a [cm] = dd4hep::mm * a [mm]
 *        cm -> mm: a [mm] = 1/dd4hep::mm * a [cm]

 * - Vectors can be given in 
 *        a) dd4hep::rec::Vector3D <- fully featured vector, overloads operators *+- etc
 *        b) edm4hep::Vector3d <- natively used by edm4hep (where simHit, digiHit are from)
 *      -> generally use dd4hep::rec::Vector3D, convert via ConvertVector() where edm4hep::Vector3d is needed
 * - Indices named i_ ... refer to pixels. Indices named j_ ... refer to in-pixel bins (for charge deposition)
 * - Reference frames: 
 *        - global detector frame, use (x,y,z)
 *             - z along beamline
 *        - local sensor frame: (u,v,w)
 *             - u,v span sensor plane, (for ARCADIA in barrel: v along z)
 *             - w normal to sensor plane
 * - Energies in keV, but deposited energy is always converted to the electron charge equivalent [e-], 3.65 eV per eh-pair
 * - Charges are given as either 
 *        - "raw charge" (as given by Geant4/Allpix2, before thresholding and noise) or 
 *        - "measured charge" (after thresholding and noise)
 */

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

  InitServicesAndGeometry();
  InitLayersAndSensors();

  /* This needs to come in after the properties, geometry and services have all been initialized */
  m_chargeCollector = VTXdigi_tools::CreateChargeCollector(*this, m_chargeCollectionMethod);

  VTXdigi_tools::PixelChargeMatrix pixelMatrix(0,0); // just to test that the class compiles

  info() << " - Initialized successfully." << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode VTXdigi_Modular::finalize() {
  info() << "FINALIZING VTXdigi_Modular..." << endmsg;
  debug() << " - finalized successfully." << endmsg;
  return StatusCode::SUCCESS;
} 


/* -- Event loop -- */

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection> VTXdigi_Modular::operator()
  (const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const {
  debug() << "STARTING event with " << simHits.size() << " simHits." << endmsg;
  
  /* TODO: */

  /* output collections */
  auto digiHits = edm4hep::TrackerHitPlaneCollection();
  auto digiHitsLinks = edm4hep::TrackerHitSimTrackerHitLinkCollection();
  
  m_chargeCollector->Collect();

  return std::make_tuple(std::move(digiHits), std::move(digiHitsLinks));
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
  m_uidService = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidService)
    throw GaudiException("Unable to get UniqueIDGenSvc", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  m_geoService = serviceLocator()->service(m_geoServiceName);
  if (!m_geoService)
    throw GaudiException("Unable to retrieve the GeoSvc", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
  
  /* TODO: get cellID without using geoService, a la Jessy (if this is advantageous) ~ Jona */
  std::string cellIDstr = m_geoService->constantAsString(m_encodingStringVariable.value());
  m_cellIdDecoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(cellIDstr);
  if (!m_cellIdDecoder)
    throw GaudiException("Unable to retrieve the cellID decoder", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
  
  m_detector = m_geoService->getDetector();
  if (!m_detector)
    throw GaudiException("Unable to retrieve the DD4hep detector from GeoSvc", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
  
  const dd4hep::rec::SurfaceManager* simSurfaceManager = m_detector->extension<dd4hep::rec::SurfaceManager>();
  if (!simSurfaceManager)
    throw GaudiException("Unable to retrieve the SurfaceManager from the DD4hep detector", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
  
  m_surfaceMap = simSurfaceManager->map(m_subDetName.value());
  if (!m_surfaceMap)
    throw GaudiException("Unable to retrieve the simSurface map for subdetector " + m_subDetName.value(), "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  m_volumeManager = m_detector->volumeManager();
  if (!m_volumeManager.isValid())
    throw GaudiException("Unable to retrieve the VolumeManager from the DD4hep detector", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  if (m_subDetName.value() == m_undefinedString)
    throw GaudiException("Property SubDetectorName is not set!", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  /* IDEA / Allegro: 
   *   - subDet is child of detector, eg. Vertex 
   *   - subDetChild is child of subDet, eg. VertexBarrel
   *   - subDetChildChild are layers 
   * If this is not the case, layers might be direct children of subDet: */

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

  debug() << " - Successfully retrieved all necessary services and dd4hep detector elements." << endmsg;
}

void VTXdigi_Modular::InitLayersAndSensors() {
  /* Get the detector & sensor geometry information from the DD4hep detector
  * 
  * Sets the members: 
  *  - m_pixelPitch 
  *  - m_sensorThickness
  *  - m_layerCount
  *  - */
  verbose() << " - Retrieving layers, subdetector geometry, and sensor size and pitch..." << endmsg;

  /* Find layers that we want to digitize */
  verbose() << "   - Determining relevant layers... " << endmsg;
  m_layerCount = 0;
  std::vector<int> availableLayers;
  for (const auto& [layerName, layer] : m_subDetector.children()) {
    dd4hep::VolumeID layerVolumeID = layer.volumeID();
    const int layerNumber = m_cellIdDecoder->get(layerVolumeID, "layer");
    availableLayers.push_back(layerNumber);
  }
  if (availableLayers.empty())
    throw GaudiException("No layers found in subdetector " + m_subDetName.value(), "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);

  if (m_layersToDigitize.value().empty()) {
    /* digitize all layers in m_subDetector */
    m_layerCount = availableLayers.size();
    m_layersToDigitize = availableLayers;
  } else {
    /* If layers-to-digitize is specified: check that all requested layers exist */
    m_layerCount = m_layersToDigitize.value().size();

    for (const auto layerNumber : m_layersToDigitize.value()) {
      if (std::find(availableLayers.begin(), availableLayers.end(), layerNumber) == availableLayers.end()) {
        throw GaudiException("Requested layer " + std::to_string(layerNumber) + " to digitize, but this layer does not exist in subdetector " + m_subDetName.value(), "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
      }
    }
    debug() << " - Digitizing only specific layers, as defined in Gaudi property by the user. All requested layers were found." << endmsg;
  }
  info() << " - Digitizing " << m_layerCount << " layers: " << m_layersToDigitize.value() << " in subdetector " << m_subDetName.value() << "." << endmsg;

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
      m_pixelPitch.at(0) = cellDimensions.at(0) * 10; // convert cm to mm
      m_pixelPitch.at(1) = cellDimensions.at(1) * 10;
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
  for (const auto& [layerKey, layer] : m_subDetector.children()) {
    dd4hep::VolumeID layerVolumeID = layer.volumeID();
    int layerNumber = m_cellIdDecoder->get(layerVolumeID, "layer");
    /* If list of layers-to-be-digitised is given: skip layers that are not on the list*/
    if (!m_layersToDigitize.value().empty()) {
      if (std::find(m_layersToDigitize.value().begin(), m_layersToDigitize.value().end(), layerNumber) == m_layersToDigitize.value().end()) {
        verbose() << " - Skipping layer " << layerKey << " (layerNumber " << layerNumber << ", volumeID " << layerVolumeID << " with " << layer.children().size() << " modules) as it is not in the LayersToDigitize list." << endmsg;
        continue;
      }
    }
    verbose() << " - Found layer \"" << layerKey << "\" (layerNumber " << layerNumber << ", volumeID " << layerVolumeID << ", " << layer.children().size() << " modules) for subDetector \"" << m_subDetName.value() << "\"." << endmsg;

    for (const auto& [moduleKey, module] : layer.children()) {
      ++moduleNumber;
      for (const auto& [sensorKey, sensor] : module.children()) {
        ++sensorNumber;
        dd4hep::VolumeID sensorVolumeID = sensor.volumeID();

        const auto surfaceIt = m_surfaceMap->find(sensor.volumeID());
        if (surfaceIt == m_surfaceMap->end()) {
          throw GaudiException("Could not find surface for sensor " + sensorKey + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layerNumber) + " of subDetector " + m_subDetName.value() + " while checking geometry consistency.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
        }
        dd4hep::rec::ISurface* surface = surfaceIt->second;
        if (!surface) {
          throw GaudiException("Surface pointer for sensor " + sensorKey + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layerNumber) + " of subDetector " + m_subDetName.value() + " is null while checking geometry consistency.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
        }

        const float sensorLength_u = surface->length_along_u() * 10; // convert to mm
        const float sensorLength_v = surface->length_along_v() * 10; 
        const float sensorThickness = (surface->innerThickness() + surface->outerThickness()) * 10;

        if (!membersDefined) {
          /* Set members based on the first sensor we find */
          m_sensorThickness = sensorThickness;
          m_sensorLength.at(0) = sensorLength_u; 
          m_sensorLength.at(1) = sensorLength_v;
          float pixelCountU = m_sensorLength.at(0) / m_pixelPitch.at(0);
          float pixelCountV = m_sensorLength.at(1) / m_pixelPitch.at(1);
          if (abs(pixelCountU - std::round(pixelCountU)) > 0.001 || abs(pixelCountV - std::round(pixelCountV)) > 0.001)
            throw GaudiException("Sensor side length (" + std::to_string(m_sensorLength.at(0)) + " x " + std::to_string(m_sensorLength.at(1)) + ") mm and pixel pitch (" + std::to_string(m_pixelPitch.at(0)) + " x " + std::to_string(m_pixelPitch.at(1)) + ") mm result in a non-integer pixel count (" + std::to_string(pixelCountU) + " x " + std::to_string(pixelCountV) + ") in subDetector " + m_subDetName.value() + ".", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
          m_pixelCount.at(0) = std::round(pixelCountU);
          m_pixelCount.at(1) = std::round(pixelCountV);
          membersDefined = true;
          verbose() << "     - Found sensor: " << sensorKey << ", volumeID: " << sensorVolumeID << " (sensor " << sensorNumber << " in layer " << layerNumber << "). Setting expected dimensions to (" << m_sensorLength.at(0) << " x " << m_sensorLength.at(1) << " x " << m_sensorThickness << ") mm3." << endmsg;
        } // if !membersDefined
        else {
          /* For all other sensors: check for consistency with first sensor */
          if (std::abs(sensorLength_u - m_sensorLength.at(0)) > 0.001 || std::abs(sensorLength_v - m_sensorLength.at(1)) > 0.001 || std::abs(sensorThickness - m_sensorThickness) > 0.001)
            throw GaudiException("Sensor dimension mismatch found in sensor " + sensorKey + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layerNumber) + " of subDetector " + m_subDetName.value() + ": expected dimensions of (" + std::to_string(m_sensorLength.at(0)) + " x " + std::to_string(m_sensorLength.at(1)) + " x " + std::to_string(m_sensorThickness) + ") mm3, but found (" + std::to_string(sensorLength_u) + " x " + std::to_string(sensorLength_v) + " x " + std::to_string(sensorThickness) + ") mm3. This algorithm expects exactly one type of sensor per subDetector. Use different instances of the algorithm if different layers consist of different sensors.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
        } // if membersDefined
      } // loop over sensors
    } // loop over modules
  } // loop over layers
  
  info() << " - Retrieved sensor parameters: area (" << m_sensorLength.at(0) << " x " << m_sensorLength.at(1) << ") mm, thickness " << m_sensorThickness << " mm, pixel pitch (" << m_pixelPitch.at(0) << " x " << m_pixelPitch.at(1) << ") mm, pixel count (" << m_pixelCount.at(0) << " x " << m_pixelCount.at(1) << "). All " << sensorNumber << " sensors in the relevant layers share these parameters." << endmsg;
} // InitLayersAndSensors()