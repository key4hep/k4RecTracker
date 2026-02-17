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
 * 
 * - I am sorry for the unholy mix of US and British spelling. I tried to keep the code consistent in US spelling. ~ Jona
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

  if (!VTXdigi_tools::ToolTest())
    return StatusCode::FAILURE;

  InitServicesAndGeometry();

  InitLayersAndSensors();

  if (m_debugHistograms)
    InitHistograms();

  /* This needs to come in after the properties, geometry and services have all been initialized */
  verbose() << "Initializing charge collection method: " << m_chargeCollectionMethod.value() << endmsg;
  m_chargeCollector = VTXdigi_tools::CreateChargeCollector(*this, m_chargeCollectionMethod);

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
  if (!CheckEventSetup(simHits, headers))
    return std::make_tuple(edm4hep::TrackerHitPlaneCollection(), edm4hep::TrackerHitSimTrackerHitLinkCollection());

  std::unordered_map<dd4hep::DDSegmentation::CellID, std::vector<VTXdigi_tools::Hit>> sensorHits; // map from sensor (shortened cellID) to hits
  for (const auto& simHit : simHits) {
    if (CheckSimhitLayer(simHit)){
      const dd4hep::DDSegmentation::CellID cellID = VTXdigi_tools::GetCellID_short(simHit);
      sensorHits[cellID].emplace_back(simHit, m_surfaceMap, m_cellIdDecoder);
      /* This copies each simHit into the hit container. 
      Because simHits are const references (rvalues), we cannot store references (pointers) to them directly in the Hit objects.
      I can't think of a significantly better way to handle this where we can delete select hits without invalidating references. 
      Hit could contain a index i that refers to the i-th simHit (instead of the simHit itself), but that would make the code too complicated (for now).
      I hope linking digiHits to simHits still works after copying the simHits. If not, we need to rethink the approach, maybe use simHitIndex in Hit.
      This is very similar to what is done in DCHdigi, so it's probably fine. ~ Jona 2026-02 */
      if (m_debugHistograms.value())
        FillHistograms_perSimHit(sensorHits[cellID].back());
    }
  }
  debug() << " - Found hits on " << sensorHits.size() << " individual sensors. Digitising sensor-wise..." << endmsg;

  auto digiHits = edm4hep::TrackerHitPlaneCollection();
  auto digiHitLinks = edm4hep::TrackerHitSimTrackerHitLinkCollection();
  std::vector<VTXdigi_tools::Hit> hitsSensor;

  /* loop over sensors */
  for (const auto& [cellID, hits] : sensorHits) {
    debug() << "   - Processing sensor with cellID " << cellID << " (layer " << hits.back().layer() << "). Has " << hits.size() << " hits." << endmsg;
  
    TGeoHMatrix trafoMatrix = VTXdigi_tools::ComputeSensorTrafoMatrix(cellID, m_volumeManager, m_sensorNormalRotation); // transformation from global to local sensor coordinates
    VTXdigi_tools::HitMap hitMap(m_pixelCount); // (integer) matrix to hold charge collected in this sensor, to be digitized at the end

    /* loop over hits on this sensor, deposit charges into hitMap */
    for (const VTXdigi_tools::Hit& hit : hits) {
      m_chargeCollector->FillHit(hit, hitMap, trafoMatrix); // use the chosen charge collection method
    }

    /* Create digiHits from the hitMap */
    if (!m_clusterize.value()) {
      for (const VTXdigi_tools::PixelHit& pix : hitMap.GetPixelsWithCharge()) {
        const dd4hep::rec::Vector3D pixPos_local = VTXdigi_tools::ComputePixelPos_local(pix.index, m_sensorLength, m_pixelPitch, m_depletedRegionDepthCenter);
        const dd4hep::rec::Vector3D pixPos_global = VTXdigi_tools::LocalToGlobal(pixPos_local, trafoMatrix);

        if (m_debugHistograms.value())
          FillHistograms_perPixel(hits.back().layer(), pix, {0,0}); // cluster center{0,0} means no clustering

        /* Match pixel hits to simHits. 
        * TODO: So far, we simply match each pixel hit to the closest simHit. At high occupancy (ttbar) this will lead to wrong associations (for 2 simHits that are very close), might deteriorate truth-matching of tracks or similar.
        * (instead, maybe save a list of all source-simHits for each pixel? This is the correct way to do it, but quite a bit more expensive. This would also mean we cannot simply remove each used hit from the hits vector) */
        const VTXdigi_tools::Hit* closestHit = &hits.at(0);
        float closestDist = std::numeric_limits<float>::max();
        bool firsthit = false;
        for (const VTXdigi_tools::Hit& hit : hits) {
          if (!firsthit) {
            firsthit = true;
            continue; // skip first hit, already set as closestHit
          }
          const dd4hep::rec::Vector3D hitPos_local = VTXdigi_tools::GlobalToLocal(VTXdigi_tools::ConvertVector(hit.simHit().getPosition()), trafoMatrix);
          const float dist = (hitPos_local - pixPos_local).r();
          if (dist < closestDist){
            closestHit = &hit;
            closestDist = dist;
          }
        }

        VTXdigi_tools::CreateDigiHit(closestHit->simHit(), digiHits, digiHitLinks, pixPos_global, pix.charge);

        FillHistograms_perDigiHit(*closestHit, pixPos_local, pix, trafoMatrix);
      } /* loop over firing pixels */
    } /* if not clusterising */
    else { 
      std::vector<VTXdigi_tools::PixelHit> pixHits = hitMap.GetPixelsWithCharge();

      /* For now: do direct-neighbor clustering.
      * TODO: implement diagonal neighbors? Could be interesting to test, at least. */

      std::vector<VTXdigi_tools::PixelHit> clusterPixs;
      clusterPixs.reserve(8); // avoid too many reallocations. 8 should include around 80-90% of clusters, i guess. 

      while (pixHits.size() > 0) {
        clusterPixs.clear();

        /* start new cluster with first pixel in list */
        clusterPixs.push_back(pixHits.back());
        pixHits.pop_back();

        /* loop over remaining pixels, add direct neighbors to cluster 
        Loop over an iterator, back to front*/


        for (auto pixIter = pixHits.begin(); pixIter != pixHits.end();) {
          /* TODO: I think there is a elegant (and more efficient, less loopy) recursive way to do this where each added pixel checks for neighbors in the set of remaining pixel hits. ~ Jona 2026-02 */
          if (IsDirectNeighbor(*pixIter, clusterPixs)) {
            clusterPixs.push_back(*pixIter);
            std::iter_swap(pixIter, pixHits.end() - 1);
            pixHits.pop_back();
          } else {
            ++pixIter; // only increment if no element was removed, otherwise we would skip the element that was just swapped into the current position
          }
        } /* loop over pixel hits (that are not clusterized yet) */

        
        const std::array<float, 2> clusterCenter = ComputeClusterPos_Weighted(clusterPixs);
        const dd4hep::rec::Vector3D clusterPos_local = VTXdigi_tools::ComputePos_local(clusterCenter, m_sensorLength, m_pixelPitch, m_depletedRegionDepthCenter);
        const dd4hep::rec::Vector3D clusterPos_global = VTXdigi_tools::LocalToGlobal(clusterPos_local, trafoMatrix);

        int clusterCharge = 0;
        for (const auto& pix : clusterPixs) {
          clusterCharge += pix.charge;
        }

        debug() << "     - Found cluster with " << clusterPixs.size() << " pixels, charge " << clusterCharge << ", center at (" << clusterCenter[0] << ", " << clusterCenter[1] << "), local pos (" << clusterPos_local.x() << ", " << clusterPos_local.y() << ", " << clusterPos_local.z() << "). " << endmsg;
        
        /* TODO: Match cluster to simHit. For now, simply match to closest simHit as for single pixels. This is probably less bad than for single pixels, because clusters are bigger and thus more likely to be close to their source simHit than single pixels are. But still not perfect. 
        TODO: implement better matching, maybe by saving a list of all source simHits for each pixel and then combining those lists for all pixels in the cluster? */

        /* TODO: create digiHit per cluster. */
        
        if (m_debugHistograms.value()) {
          // FillHistograms_perDigiHit(...);
          for (const VTXdigi_tools::PixelHit& pix : clusterPixs)
            FillHistograms_perPixel(hits.back().layer(), pix, clusterCenter); // per-pixel histograms
        } /* debug histograms */

      } /* loop over clusters */
    } /* if clusterizing */
  } /* loop over sensors */

  /* TODO: Impletement FAST charge collection method (simply smear simHit position)
  * I would simply create a function that returns digiHits and digiHitLinks directly from simHits, to be returned here.
  * this makes a lot of the init etc unnecessary, so maybe have a separate class for that? ~ Jona 2026-02 */

  debug() << " - Finished digitization. Created " << digiHits.size() << " digiHits." << endmsg;
  return std::make_tuple(std::move(digiHits), std::move(digiHitLinks));
} // operator()



bool IsDirectNeighbor(const VTXdigi_tools::PixelHit& pix, std::vector<VTXdigi_tools::PixelHit>& clusterPixs) {
  for (const auto& clusterPix : clusterPixs) {
    const int di_u = std::abs(pix.index[0] - clusterPix.index[0]);
    const int di_v = std::abs(pix.index[1] - clusterPix.index[1]);

    if ((di_u == 1 && di_v == 0) || (di_u == 0 && di_v == 1)) {
      return true;
    }
  }
  return false;
}

std::array<float, 2> ComputeClusterPos_Weighted(const std::vector<VTXdigi_tools::PixelHit>& clusterPixs) {
  std::array<float, 2> pos{0.f, 0.f};
  int totalCharge = 0;
  for (const auto& pix : clusterPixs) {
    pos[0] += pix.index[0] * pix.charge;
    pos[1] += pix.index[1] * pix.charge;
    totalCharge += pix.charge;
  }
  if (totalCharge <= 0)
    throw std::runtime_error("ComputeClusterPos_Weighted: total charge of cluster is 0 or negative, cannot compute weighted position");
  pos[0] /= totalCharge;
  pos[1] /= totalCharge;
  return pos;
}





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
  if (m_subDetName.value() == m_undefinedString)
    throw GaudiException("Property SubDetectorName is not set!", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

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
    * It does NOT make sure that the local x and y are not swapped and have the correct polarity
    * I do NOT have the energy to think about Euler angles now. ~Jona, 2026-02 
    * This at least prints a warning if it isn't correct. */
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
          m_sensorLength.at(0) = sensorLength_u; 
          m_sensorLength.at(1) = sensorLength_v;
          float pixelCountU = m_sensorLength.at(0) / m_pixelPitch.at(0);
          float pixelCountV = m_sensorLength.at(1) / m_pixelPitch.at(1);
          if (abs(pixelCountU - std::round(pixelCountU)) > 0.001 || abs(pixelCountV - std::round(pixelCountV)) > 0.001)
            throw GaudiException("Sensor side length (" + std::to_string(m_sensorLength.at(0)) + " x " + std::to_string(m_sensorLength.at(1)) + ") mm and pixel pitch (" + std::to_string(m_pixelPitch.at(0)) + " x " + std::to_string(m_pixelPitch.at(1)) + ") mm result in a non-integer pixel count (" + std::to_string(pixelCountU) + " x " + std::to_string(pixelCountV) + ") in subDetector " + m_subDetName.value() + ".", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
          m_pixelCount.at(0) = std::round(pixelCountU);
          m_pixelCount.at(1) = std::round(pixelCountV);
          membersDefined = true;
          verbose() << "     - Found sensor: " << sensorKey << ", volumeID: " << sensorVolumeID << " (sensor " << sensorNumber << " in layer " << layer << "). Setting expected dimensions to (" << m_sensorLength.at(0) << " x " << m_sensorLength.at(1) << " x " << m_sensorThickness << ") mm3." << endmsg;
        } // if !membersDefined
        else {
          /* For all other sensors: check for consistency with first sensor */
          if (std::abs(sensorLength_u - m_sensorLength.at(0)) > 0.001 || std::abs(sensorLength_v - m_sensorLength.at(1)) > 0.001 || std::abs(sensorThickness - m_sensorThickness) > 0.001)
            throw GaudiException("Sensor dimension mismatch found in sensor " + sensorKey + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layer) + " of subDetector " + m_subDetName.value() + ": expected dimensions of (" + std::to_string(m_sensorLength.at(0)) + " x " + std::to_string(m_sensorLength.at(1)) + " x " + std::to_string(m_sensorThickness) + ") mm3, but found (" + std::to_string(sensorLength_u) + " x " + std::to_string(sensorLength_v) + " x " + std::to_string(sensorThickness) + ") mm3. This algorithm expects exactly one type of sensor per subDetector. Use different instances of the algorithm if different layers consist of different sensors.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
        } // if membersDefined
      } // loop over sensors
    } // loop over modules
  } // loop over layers
  
  info() << " - Retrieved sensor parameters: area (" << m_sensorLength.at(0) << " x " << m_sensorLength.at(1) << ") mm, thickness " << m_sensorThickness << " mm, pixel pitch (" << m_pixelPitch.at(0) << " x " << m_pixelPitch.at(1) << ") mm, pixel count (" << m_pixelCount.at(0) << " x " << m_pixelCount.at(1) << "). All " << sensorNumber << " sensors in the relevant layers share these parameters." << endmsg;
} // InitLayersAndSensors()

void VTXdigi_Modular::InitHistograms() {
  /* Define axes globally to make adjusting them easier
  * TODO: Make some of these adjustable via Gaudi Parameters? Might not be necessary.*/
  Gaudi::Accumulators::Axis<float> axis_z{200, -200, 200};
  Gaudi::Accumulators::Axis<float> axis_cosTheta{100, 0, 1};
  Gaudi::Accumulators::Axis<float> axis_theta{4*180, 0, 180};
  Gaudi::Accumulators::Axis<float> axis_phi{4*180, -180, 180};
  
  Gaudi::Accumulators::Axis<float> axis_moduleID{2000, -0.5f, 1999.5f};
  Gaudi::Accumulators::Axis<float> axis_clusterSize{30, -0.5f, 29.5f};
  Gaudi::Accumulators::Axis<float> axis_E{1000, 0, m_sensorThickness*2000.f};
  Gaudi::Accumulators::Axis<float> axis_charge{1000, 0, m_sensorThickness*500000.f};
  Gaudi::Accumulators::Axis<float> axis_momentum_keV{10000, 0.f, 1000.f};
  Gaudi::Accumulators::Axis<float> axis_momentum_MeV{10000, 0.f, 1000.f};
  Gaudi::Accumulators::Axis<float> axis_momentum_GeV{10000, 0.f, 1000.f};
  
  Gaudi::Accumulators::Axis<float> axis_residual{2000, -1000.f, 1000.f};
  Gaudi::Accumulators::Axis<float> axis_residual_pixels{2020, -50.5f, 50.5f}; // 20 in-pix bins
  Gaudi::Accumulators::Axis<float> axis_pdg{1401, -700.5f, 700.5f};

  Gaudi::Accumulators::Axis<float> axis_pixels_u{
    static_cast<unsigned int>(m_pixelCount.at(0)),
    -0.5f,
    static_cast<float>(m_pixelCount.at(0)+0.5)};
  Gaudi::Accumulators::Axis<float> axis_pixels_v{
    static_cast<unsigned int>(m_pixelCount.at(1)),
    -0.5f,
    static_cast<float>(m_pixelCount.at(1)+0.5)};
    
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
        "Layer" + std::to_string(layer) + "/digiHits_per_simHit",
        "Raw number of digiHits created per simHit - Layer " + std::to_string(layer) + ";DigiHits per SimHit;Entries",
        axis_clusterSize // not technically correct, but works
      }
    );

    hist1d.at(hist1d_clusterSize).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/clusterSize",
        "Number of pixels per cluster (after clustering algorithm) - Layer " + std::to_string(layer) + ";Cluster size [pixels];Entries",
        axis_clusterSize
      }
    );
    hist1d.at(hist1d_clusterSize_createdInGenerator).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/clusterSize_createdInGenerator",
        "Number of pixels per cluster (for simHits created in generator) - Layer " + std::to_string(layer) + ";Cluster size [pixels];Entries",
        axis_clusterSize
      }
    );
    hist1d.at(hist1d_clusterSize_createdInSimulation).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/clusterSize_createdInSimulation",
        "Number of pixels per cluster (for simHits created in simulation) - Layer " + std::to_string(layer) + ";Cluster size [pixels];Entries",
        axis_clusterSize
      }
    );

    hist1d.at(hist1d_residualU).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_u",
        "Residual (u_digiHit - u_simHit) in local u direction - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residualV).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_v",
        "Residual (v_digiHit - v_simHit) in local v direction - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residualW).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_w",
        "Residual (w_digiHit - w_simHit) in local w direction (ie. sensor normal direction) - Layer " + std::to_string(layer) + ";Residual w [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residualR).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residual_r",
        "Residual magnitude (|r_digiHit - r_simHit|) - Layer " + std::to_string(layer) + ";Residual r [um];Entries",
        axis_residual
      }
    );

    hist1d.at(hist1d_pixelDistToClusterCenterU).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/pixelHit_distanceToClusterCenter_u",
        "Distance of pixel to center of cluster in u direction - Layer " + std::to_string(layer) + ";Distance to cluster center u [pix];Entries",
        axis_residual_pixels
      }
    );
    hist1d.at(hist1d_pixelDistToClusterCenterV).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/pixelHit_distanceToClusterCenter_v",
        "Distance of pixel to center of cluster in v direction - Layer " + std::to_string(layer) + ";Distance to cluster center v [pix];Entries",
        axis_residual_pixels
      }
    );

    m_hist1d.emplace(layer, std::move(hist1d));

    std::array< std::unique_ptr< Gaudi::Accumulators::StaticHistogram< 2, Gaudi::Accumulators::atomicity::full, float > >, hist2dArrayLen>  hist2d;

    hist2d.at(hist2d_hitMap_simHits).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_hitMap",
        "SimHit hit map (pixel in which this simHit position lies) - Layer " + std::to_string(layer) + ";Pixel u;Pixel v;Entries",
        axis_pixels_u,
        axis_pixels_v
      }
    );
    hist2d.at(hist2d_hitMap_pixelHits).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/hitMap",
        "Hit map (pixels that collect charge over the threshold) - Layer " + std::to_string(layer) + ";Pixel u;Pixel v;Entries",
        axis_pixels_u,
        axis_pixels_v
      }
    );
    hist2d.at(hist2d_clusterSize_vs_hit_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/clusterSize_vs_hit_z_2D",
        "Cluster size vs. global hit z - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Pixels per cluster;Entries",
        axis_z,
        axis_clusterSize
      }
    );
    hist2d.at(hist2d_clusterSize_vs_hit_z_createdInGenerator).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/clusterSize_vs_hit_z_createdInGenerator_2D",
        "Cluster size vs. global hit z (simHit particle created in generator) - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Pixels per cluster;Entries",
        axis_z,
        axis_clusterSize
      }
    );
    hist2d.at(hist2d_clusterSize_vs_hit_z_createdInSim).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/clusterSize_vs_hit_z_createdInSimulation_2D",
        "Cluster size vs. global hit z (simHit particle created in simulation) - Layer " + std::to_string(layer) + ";SimHit global z position [mm];Pixels per cluster;Entries",
        axis_z,
        axis_clusterSize
      }
    );

    m_hist2d.emplace(layer, std::move(hist2d));

  } /* loop over layers */
}

/* ---- Eventloop functions ---- */

bool VTXdigi_Modular::CheckEventSetup(const edm4hep::SimTrackerHitCollection& simHits, const edm4hep::EventHeaderCollection& headers) const {
  if (m_counter_eventsRead.value() % m_infoPrintInterval.value() == 0)
  info() << "PROCESSING event [run " << headers.at(0).getRunNumber() << ", event " << headers.at(0).getEventNumber() << ", found " << simHits.size() << " simHits]. " << m_counter_eventsRead.value() << " events so far." << endmsg;
  /* events are not necessarily numbered sequentially... */
  ++m_counter_eventsRead;

  if (simHits.size()==0) {
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
      verbose() << " - Dismissing simHit in layer " << layer << ". (not in the list of layers to digitize)" << endmsg;
      ++m_counter_simHitsRejected_LayerNotToBeDigitized;
      return false;
    }
  }
  else {
    verbose() << " - All layers are to be digitized, as property \"LayersToDigitize\" is not set ." << endmsg;
  }

  ++m_counter_simHitsAccepted;
  return true;
}

void VTXdigi_Modular::FillHistograms_perSimHit(const VTXdigi_tools::Hit& hit) const {
  /* executed once for each simHit */

  const int layer = hit.layer();
  const dd4hep::rec::Vector3D simHitMomentum = VTXdigi_tools::ConvertVector(hit.simHit().getMomentum()); // in GeV
  const dd4hep::rec::Vector3D simHitPos_global = VTXdigi_tools::ConvertVector(hit.simHit().getPosition());
  const TGeoHMatrix trafoMatrix = VTXdigi_tools::ComputeSensorTrafoMatrix(hit.cellID(), m_volumeManager, m_sensorNormalRotation);
  const dd4hep::rec::Vector3D simHitPos_local = VTXdigi_tools::GlobalToLocal(simHitPos_global, trafoMatrix);
  const std::array<int, 2> pix_i = VTXdigi_tools::ComputePixelIndices(simHitPos_local, m_pixelPitch, m_pixelCount);

  ++(*m_hist1d.at(layer).at(hist1d_simHitE))[hit.simHit().getEDep() * (dd4hep::GeV / dd4hep::keV)]; 
  ++(*m_hist1d.at(layer).at(hist1d_simHitCharge))[hit.charge()]; 
  ++(*m_hist1d.at(layer).at(hist1d_simHitMomentum_keV))[simHitMomentum.r() * 1.E6]; 
  ++(*m_hist1d.at(layer).at(hist1d_simHitMomentum_MeV))[simHitMomentum.r() * 1.E3]; 
  ++(*m_hist1d.at(layer).at(hist1d_simHitMomentum_GeV))[simHitMomentum.r()]; 
  ++(*m_hist1d.at(layer).at(hist1d_simHitPDG))[hit.simHit().getParticle().getPDG()];

  ++(*m_hist2d.at(layer).at(hist2d_hitMap_simHits))[{pix_i.at(0), pix_i.at(1)}];
}

void VTXdigi_Modular::FillHistograms_perPixel(const int layer, const VTXdigi_tools::PixelHit& pix, const std::array<float, 2> clusterPos_local) const {
  /* executed once for each pixel */

  ++(*m_hist2d.at(layer).at(hist2d_hitMap_pixelHits))[{pix.index.at(0), pix.index.at(1)}];

  if (clusterPos_local[0] != 0 && clusterPos_local[1] != 0) { 
    /* cluster at (0,0) means that hits weren't clustered */
    const float distToClusterCenter_u = static_cast<float>(pix.index.at(0)) - clusterPos_local[0]; // in pix
    const float distToClusterCenter_v = static_cast<float>(pix.index.at(1)) - clusterPos_local[1];

    ++(*m_hist1d.at(layer).at(hist1d_pixelDistToClusterCenterU))[distToClusterCenter_u]; // convert to um
    ++(*m_hist1d.at(layer).at(hist1d_pixelDistToClusterCenterV))[distToClusterCenter_v];
  }
}

void VTXdigi_Modular::FillHistograms_perDigiHit(const VTXdigi_tools::Hit& hit, const dd4hep::rec::Vector3D& pos_local, const VTXdigi_tools::PixelHit& pix, const TGeoHMatrix& trafoMatrix) const {
  /* executed once for each digiHit */

  const int layer = hit.layer();
  const dd4hep::rec::Vector3D simHitPos_global = VTXdigi_tools::ConvertVector(hit.simHit().getPosition());
  const dd4hep::rec::Vector3D simHitPos_local = VTXdigi_tools::GlobalToLocal(simHitPos_global, trafoMatrix);
  const dd4hep::rec::Vector3D residual_local = pos_local - simHitPos_local; // residual = observed - predicted

  ++(*m_hist1d.at(layer).at(hist1d_digiHitCharge))[pix.charge]; 
  ++(*m_hist1d.at(layer).at(hist1d_residualU))[residual_local.x() / dd4hep::um]; // convert to um
  ++(*m_hist1d.at(layer).at(hist1d_residualV))[residual_local.y() / dd4hep::um];
  ++(*m_hist1d.at(layer).at(hist1d_residualW))[residual_local.z() / dd4hep::um];
  ++(*m_hist1d.at(layer).at(hist1d_residualR))[residual_local.r() / dd4hep::um];
}

/* TODO: add */
// hist1d_digiHitsPerSimHit,
// hist1d_clusterSize,
// hist1d_clusterSize_createdInGenerator,
// hist1d_clusterSize_createdInSimulation,
// hist2d_clusterSize_vs_hit_z,
// hist2d_clusterSize_vs_hit_z_createdInGenerator,
// hist2d_clusterSize_vs_hit_z_createdInSim,
// hist2d_clusterSize_vs_module_z,