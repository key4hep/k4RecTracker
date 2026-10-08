// VTXdigi_Modular/src/VTXdigi_Modular.cpp
#include "VTXdigi_Modular.h"
#include "VTXdigi_tools.h"
#include <DDRec/CellIDPositionConverter.h>
#include <GaudiKernel/RndmGenerators.h>

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

  InitServicesAndGeometry();

  InitLayersAndSensors();

  // This needs to come in after the properties, geometry and services have all been initialized
  verbose() << "Initializing charge collection method: " << m_chargeCollectionMethod.value() << endmsg;
  m_chargeCollector = VTXdigi_tools::CreateChargeCollector(*this, m_chargeCollectionMethod);

  if (m_debugHistograms)
  // needs to run after charge collector is initialised
    InitHistograms();

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

  // seed RNG (done per-evt so it is thread-safe & reproducible)
  const auto uniqueID = m_uniqueIDService->getUniqueID(headers, name());
  TRandom3 randomGen(uniqueID);

  /* TODO: Implement fast digitization, where simHit pos is simply smeared by a Gaussian. But:
  * that would make a lot of the init etc unnecessary, so maybe keep it in a separate class? We cannot simply add a ChargeCollector implementation to do this, because ChargeCollectors act on pixels, but fast digitization acts on the position itself (ofc we could emulate this in clusters, but thats incredibly inefficient and stinks) */

  std::unordered_map<dd4hep::DDSegmentation::VolumeID, std::vector<VTXdigi_tools::SimHitWrapper>> sensorSimHits; // map from sensor (volumeID) to hits
  for (const edm4hep::SimTrackerHit& simTrackerHit : simTrackerHits) {
    if (CheckSimhitLayer(simTrackerHit)){
      const dd4hep::DDSegmentation::VolumeID volumeID = GetVolumeID(simTrackerHit.getCellID());
      if (volumeID == 0) {
        info() << "SimTrackerHit in event " << headers.at(0).getEventNumber() << " with cellID " << simTrackerHit.getCellID() << " is not in a valid sensor volume. Skipping this hit." << endmsg;
        continue;
      }

      sensorSimHits[volumeID].emplace_back(simTrackerHit, volumeID, *this); // simTrackerHits are copied here. Pointers to these are passed around (eg. in hit/pixel/cluster objects).

      switch (sensorSimHits[volumeID].back().mcParticleLevel()) {
        case VTXdigi_tools::MCParticleLevel::Primary:
          ++m_counter_accSimHitsFromPrimary;
          break;
        case VTXdigi_tools::MCParticleLevel::Secondary:
          ++m_counter_accSimHitsFromSecondary;
          break;
        case VTXdigi_tools::MCParticleLevel::Delta:
          ++m_counter_accSimHitsFromDelta;
      }
    }
  }
  debug() << " - Found " << simTrackerHits.size() << " simTrackerHits on " << sensorSimHits.size() << " individual sensors. Digitising sensor-wise..." << endmsg;

  auto digiHits = edm4hep::TrackerHitPlaneCollection();
  auto digiHitLinks = edm4hep::TrackerHitSimTrackerHitLinkCollection();
  std::vector<VTXdigi_tools::SimHitWrapper> hitsSensor;

  /* loop over sensors */
  for (const auto& [volumeID, simHits] : sensorSimHits) {
    debug() << "   - Processing sensor with volumeID " << volumeID << " (layer " << simHits.back().layer() << "). Has " << simHits.size() << " simHits." << endmsg;

    const TGeoHMatrix trafoMatrix = VTXdigi_tools::ComputeSensorTrafoMatrix(volumeID, m_volumeManager, m_sensorNormalRotation); // transformation from global detector to local sensor frame

    /* Fill hit map with charge from all simHits on this sensor, using the selected charge collection method */
    VTXdigi_tools::HitMap hitMap(m_pixelCount);

    for (const VTXdigi_tools::SimHitWrapper& simHit : simHits) {
      const dd4hep::rec::Vector3D pos_global = VTXdigi_tools::ConvertVector(simHit.hitPtr()->getPosition());
      simHit.SetTruthPos(VTXdigi_tools::Trafo_global_local(pos_global, trafoMatrix)); // do this only now to not compute the trafo matrix twice. TruthPos might be shifted by the charge collection algorithm later.

      if (!VTXdigi_tools::IsInsideVolume(simHit.truthPos(), ActiveVolumeDimensions(), VTXdigi_tools::kSensorBoundaryTolerance)) [[unlikely]] {
        warning() << "SimTrackerHit in event " << headers.at(0).getEventNumber() << " with cellID " << simHit.hitPtr()->getCellID() << " has truth position (" << simHit.truthPos().x() << ", " << simHit.truthPos().y() << ", " << simHit.truthPos().z() << ") local, which is outside the sensor volume. Skipping this hit." << endmsg;
        continue;
      }

      if (!VTXdigi_tools::IsInsideVolume(simHit.truthPos(), ActiveVolumeDimensions(), 0.0)) [[unlikely]] {
        const dd4hep::rec::Vector3D pos_clamped = VTXdigi_tools::ClampToVolume(simHit.truthPos(), ActiveVolumeDimensions());
        debug() << "     - Clamping simHit truth position from (" << simHit.truthPos().x() << ", " << simHit.truthPos().y() << ", " << simHit.truthPos().z() << ") local to sensor volume." << endmsg;
        simHit.SetTruthPos(pos_clamped); // clamp to avoid issues with binning (this also avoids checks in hot loops)
      }

      debug() << "     - Processing simHit, charge dep. " << simHit.charge() << " e at (" << simHit.truthPos().x() << ", " << simHit.truthPos().y() << ", " << simHit.truthPos().z() << ") local, (" << pos_global.x() << ", " << pos_global.y() << ", " << pos_global.z() << ") global" << endmsg;

      m_chargeCollector->FillHit(simHit, hitMap, trafoMatrix, randomGen); // uses the selected charge collection method

      if (m_debugHistograms.value())
        FillHistograms_perSimHit(simHit);
    }

    if (m_smearing_charge.value() > 0.f)
      hitMap.ApplyChargeSmearing(m_smearing_charge.value(), randomGen);
    if (m_threshold.value() > 0.f)
      hitMap.ApplyThreshold(m_threshold.value(), m_smearing_threshold.value(), randomGen);

    std::vector<VTXdigi_tools::Cluster> clusters = Clusterize(hitMap);

    CreateDigiHits(digiHits, digiHitLinks, volumeID, trafoMatrix, clusters, randomGen);
    if (m_debugHistograms.value())
      FillHistograms_perSensor(simHits, digiHits, trafoMatrix, volumeID);
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

  if (m_threshold.value() <= 5 * sqrt(m_smearing_charge.value()*m_smearing_charge.value() + m_smearing_threshold.value()*m_smearing_threshold.value()))
    warning() << "Threshold " << m_threshold.value() << " e- is less than 5 times the charge smearing and threshold dispersion sigma: sqrt(" << m_smearing_charge.value() << "^2 + " << m_smearing_threshold.value() << "^2) e-. This digitiser only applies smearing to pixels that have collected charge from simHits (a tiny bit is enough), so it does not simulate random firing of pixels. (doing this by drawing a noise for every pixel in the detector for every event would be INCREDIBLY slow. A work-around to simulate random pixels firing might be implemented in another algorith)." << endmsg;
  if (m_smearing_time.value() < 0.f)
    throw GaudiException("Time smearing sigma " + std::to_string(m_smearing_time.value()) + " ns is negative.", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  if (m_positionUncertainty.value().empty())
    info() << "Cluster position uncertainty not set. Using pitch/12 for every cluster." << endmsg;
  else if (m_positionUncertainty.value().size() == 2)
    info() << "Cluster position uncertainty set to (" << m_positionUncertainty.value().at(0) << " mm, " << m_positionUncertainty.value().at(1) << " mm) in u and v direction." << endmsg;
  else if (m_positionUncertainty.value().size() == 10) {
    std::string unc_u = "[";
    std::string unc_v = "[";
    for (int size = 0; size < 5; size++) {
      unc_u += std::to_string(m_positionUncertainty.value().at(size)) + ", ";
      unc_v += std::to_string(m_positionUncertainty.value().at(size + 5)) + ", ";
    }
    unc_u += "]";
    unc_v += "]";

    info() << "Cluster position uncertainty set to " << unc_u << " mm and " << unc_v << " mm for cluster lengths of (1, 2, 3, 4, 5+), in u and v direction, respectively." << endmsg;
  }
  else
    throw GaudiException("Property ClusterPositionUncertainty must be either empty (assign pitch/sqrt(12)), have exactly 2 values (for fixed uncertainty in u and v), or 10 values (for cluster-length based uncertainty).", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  m_uniqueIDService = service("uidSvc", false);
  if (!m_uniqueIDService)
    throw GaudiException("Unable to get UniqueIDGenSvc from name 'uidSvc'.", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  if (m_subDetName.value() == m_undefinedString)
    throw GaudiException("Property SubDetectorName is not set!", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  m_geoService = serviceLocator()->service(m_geoServiceName);
  if (!m_geoService)
    throw GaudiException("Unable to retrieve the GeoSvc", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  /* TODO: get volumeID without using geoService, a la Jessy (if this is advantageous) */
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
  if (!m_surfaceMap) {
    debug() << "Available surface maps: " << endmsg;
    debug() << surfaceManager->toString() << endmsg;
    throw GaudiException("Unable to retrieve the simSurface map for subdetector " + m_subDetName.value() + " (printed availabe surfaces above at debug level)", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
  }

  m_volumeManager = m_detector->volumeManager();
  if (!m_volumeManager.isValid())
    throw GaudiException("Unable to retrieve the VolumeManager from the DD4hep detector", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  m_cellIDPositionConverter = std::make_unique<dd4hep::rec::CellIDPositionConverter>(*m_detector);
  if (!m_cellIDPositionConverter)
    throw GaudiException("Unable to create CellIDPositionConverter", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

  { /* DD4hep has a transformation from the global detector coordinates to each sensors local system. The definition of the local system might change.
    * We define a rotation matrix to align the sensor surface axes (u,v,n) with the local system axes (x_local, y_local, z_local) for the first sensor we find, then apply this to all sensors.
    * Note the local system axes are usually referred to as (u, v, w) == (x_local, y_local, z_local) */

    auto surfaceMapIter = m_surfaceMap->begin();
    if (surfaceMapIter == m_surfaceMap->end())
      throw GaudiException("Surface map for subdetector " + m_subDetName.value() + " is empty.", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
    const dd4hep::rec::VolumeID volumeID = surfaceMapIter->first;
    const dd4hep::rec::ISurface* surface = surfaceMapIter->second;

    /* first: rotate sensor local coordinates st. u=x_local, v=y_local, n=z_local */
    const double epsilon = 1.0e-6; // reasonable for comparing to 1

    TGeoHMatrix sensorTrafoMatrix = m_volumeManager.lookupDetElement(volumeID).nominal().worldTransformation();
    const dd4hep::rec::Vector3D u = VTXdigi_tools::TrafoVec_global_local(surface->u(), sensorTrafoMatrix).unit();
    const dd4hep::rec::Vector3D v = VTXdigi_tools::TrafoVec_global_local(surface->v(), sensorTrafoMatrix).unit();
    const dd4hep::rec::Vector3D n = VTXdigi_tools::TrafoVec_global_local(surface->normal(), sensorTrafoMatrix).unit();

    /* check 1 - u,v,n are orthogonal and right-handed */
    if (std::abs(u.dot(v)) > epsilon || std::abs(u.dot(n)) > epsilon || std::abs(v.dot(n)) > epsilon)
      throw GaudiException("Sensor surface axes u " + VTXdigi_tools::VectorToString(u) + ", v " + VTXdigi_tools::VectorToString(v) + ", n " + VTXdigi_tools::VectorToString(n) + " are not orthogonal.", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
    if ((u.cross(v) - n).r() > epsilon)
      throw GaudiException("Sensor surface axes are left-handed (u x v != n). Cannot build a rotation for u " + VTXdigi_tools::VectorToString(u) + ", v " + VTXdigi_tools::VectorToString(v) + ", n " + VTXdigi_tools::VectorToString(n) + ".", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);

    /* assemble rotation matrix*/
    const double rot[9] = { u.x(), v.x(), n.x(),
                            u.y(), v.y(), n.y(),
                            u.z(), v.z(), n.z() };
    m_sensorNormalRotation.SetMatrix(rot);
    debug() << "   - Sensor local axes are u " << VTXdigi_tools::VectorToString(u) << ", v " << VTXdigi_tools::VectorToString(v) << ", n " << VTXdigi_tools::VectorToString(n) << ". Rotating them onto axes u-x, v-y, n-z in sensor local coordinates ." << endmsg;


    /* check 2 - rotation matrix works (this checks my maths) */
    const TGeoHMatrix sensorTrafoMatrix_corrected = VTXdigi_tools::ComputeSensorTrafoMatrix(volumeID, m_volumeManager, m_sensorNormalRotation);
    const dd4hep::rec::Vector3D u_corrected = VTXdigi_tools::TrafoVec_global_local(surface->u().unit(), sensorTrafoMatrix_corrected);
    const dd4hep::rec::Vector3D v_corrected = VTXdigi_tools::TrafoVec_global_local(surface->v().unit(), sensorTrafoMatrix_corrected);
    const dd4hep::rec::Vector3D n_corrected = VTXdigi_tools::TrafoVec_global_local(surface->normal().unit(), sensorTrafoMatrix_corrected);
    if ((u_corrected - dd4hep::rec::Vector3D(1., 0., 0.)).r() > epsilon || (v_corrected - dd4hep::rec::Vector3D(0., 1., 0.)).r() > epsilon || (n_corrected - dd4hep::rec::Vector3D(0., 0., 1.)).r() > epsilon)
      throw GaudiException("After rotation, sensor local axes are u " + VTXdigi_tools::VectorToString(u_corrected) + ", v " + VTXdigi_tools::VectorToString(v_corrected) + ", n " + VTXdigi_tools::VectorToString(n_corrected) + ". Expected u (1,0,0), v (0,1,0), n (0,0,1).", "VTXdigi_Modular::InitServicesAndGeometry()", StatusCode::FAILURE);
    debug() << "   - After rotation, sensor local axes are u " << VTXdigi_tools::VectorToString(u_corrected) << ", v " << VTXdigi_tools::VectorToString(v_corrected) << ", n " << VTXdigi_tools::VectorToString(n_corrected) << "." << endmsg;
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
  *  - m_sensorActiveThickness
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
  info() << " - Digitizing " << m_layers.value().size() << " layers: " << m_layers.value() << " in subdetector " << m_subDetName.value() << "." << endmsg;

  /* Get pixel pitch from the segmentation (from the readout that matches our simHitCollection)
  *  I don't know of a better way to do this (that also works for the IDEA detector model) */
  const auto simHitLocations = inputLocations("SimTrackHitCollectionName");
  if (simHitLocations.size() != 1)
    throw GaudiException("SimTrackHitCollectionName must contain exactly one collection name, got " + std::to_string(simHitLocations.size()), "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
  const std::string simHitCollectionName = simHitLocations.at(0);
  // doing this->getProperty("SimTrackHitCollectionName", simHitCollectionName) returns a list in string-format which is horrible to use

  const dd4hep::Detector::HandleMap readoutHandleMap = m_detector->readouts();
  if (readoutHandleMap.find(simHitCollectionName) == readoutHandleMap.end())
    throw GaudiException("Could not find readout matching SimTrackHitCollectionName \"" + simHitCollectionName + "\" in detector while checking geometry consistency.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
  debug() << "     - Found readout \"" << simHitCollectionName << "\". Getting segmentation." << endmsg;

  const dd4hep::Segmentation& segmentation = m_detector->readout(simHitCollectionName).segmentation();
  if (!segmentation.isValid())
    throw GaudiException("Segmentation for readout " + simHitCollectionName + " is not valid.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
  const auto cellDimensions = segmentation.cellDimensions(0); // this assumes all cells have the same dimensions (ie. only one sensor type in this readout)
  m_pixelPitch[0] = cellDimensions.at(0) * 10; // convert cm to mm
  m_pixelPitch[1] = cellDimensions.at(1) * 10;

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

        // FIRST: find thickness of the active sensor volume
        // using the active volume solid of this sensor
        dd4hep::VolumeID sensorVolumeID = sensorObj.volumeID();
        dd4hep::Volume sensorVolume = sensorObj.volume();
        std::array<double, 3> solidDimensions{}; // dimensions of active volume. indices must not correspond to local sensor (u,v,w) axes, but might be swapped

        // sensorVolume.solid() returns a dd4hep::Solid_type<T> object, generalised as dd4hep::Solid. This can either be a box or a trapezoid (for sensors in IDEA / ALLEGRO)
        try {
          dd4hep::Box sensorBox = sensorVolume.solid(); // directions 0,1,2 do not necessarily correspond to the local sensor u,v,w directions, but might be swapped.
          solidDimensions[0] = sensorBox.x() * 2 * 10;
          solidDimensions[1] = sensorBox.y() * 2 * 10;
          solidDimensions[2] = sensorBox.z() * 2 * 10;

          verbose() << "     - Sensor \"" << sensorKey << "\" (layer " << layerKey << ", module " << moduleKey << ", volumeID " << sensorVolumeID << ") has a dd4hep::Box solid." << endmsg;
        }
        catch (...) {
          try {
            dd4hep::Trd1 sensorTrd1 = sensorVolume.solid();
            solidDimensions[0] = (sensorTrd1.dX1() + sensorTrd1.dX2()) / 2 * 2 * 10; // average of upper and lower base of trapezoid. Convert half-length in cm to full length in mm
            solidDimensions[1] = sensorTrd1.dY() * 2 * 10; // Convert half-length in cm to full length in mm
            solidDimensions[2] = sensorTrd1.dZ() * 2 * 10;
            // Note: there is some weirdness in dX1() and dX2() with Trd1. I did not dig into this. Assume that these might be a bit funky.
            verbose() << "     - Sensor \"" << sensorKey << "\" (layer " << layerKey << ", module " << moduleKey << ", volumeID " << sensorVolumeID << ") has a dd4hep::Trd1 solid." << endmsg;
          }
          catch (...) {
            throw GaudiException("Unknown sensor solid type found (neither dd4hep::Box nor dd4hep::Trd1).", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
          }
        }

        const uint thicknessIndex = std::distance(solidDimensions.begin(), std::min_element(solidDimensions.begin(), solidDimensions.end()));
        const double solidThickness = solidDimensions.at(thicknessIndex);
        const double solidLength_0 = solidDimensions.at((thicknessIndex+1) % 3);
        const double solidLength_1 = solidDimensions.at((thicknessIndex+2) % 3);

        // SECOND: find lengths of the sensor (these will match the segmentation), and amount of inactive material above and below
        // using the sensitive surface of this sensor (which lies in the middle of the active volume)
        const auto surfaceIt = m_surfaceMap->find(sensorObj.volumeID());
        if (surfaceIt == m_surfaceMap->end()) {
          throw GaudiException("Could not find surface for sensor " + sensorKey + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layer) + " of subDetector " + m_subDetName.value() + " while checking geometry consistency.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
        }
        dd4hep::rec::ISurface* surface = surfaceIt->second;
        if (!surface) {
          throw GaudiException("Surface pointer for sensor " + sensorKey + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layer) + " of subDetector " + m_subDetName.value() + " is null while checking geometry consistency.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
        }

        const double surfaceLength_u = surface->length_along_u() * 10; // convert cm to mm
        const double surfaceLength_v = surface->length_along_v() * 10;
        const double surfaceThickness_above = surface->outerThickness() * 10; // sensor thickness measured from w=0 upwards, including inactive material above the active volume.
        const double surfaceThickness_below = surface->innerThickness() * 10; // same, but below
        // Note: the sensor local coordinate system (u,v,w) is centered on the active volume, so inactive material upper/lower might be assymetric

        // THIRD: consistency checks on sensor dimensions
        if (solidThickness > (surfaceThickness_above + surfaceThickness_below) + 1e-10) {
          throw GaudiException("Solid sensor thickness " + std::to_string(solidThickness) + " mm is larger than total sensor thickness " + std::to_string(surfaceThickness_above) + " + " + std::to_string(surfaceThickness_below) + " = " + std::to_string(surfaceThickness_above + surfaceThickness_below) + " mm (including inactive material) in sensor " + sensorKey + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layer) + " of subDetector " + m_subDetName.value() + ". This indicates an inconsistency in the geometry description.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
        }

        float tolerance = 0.01f; // 10 um should be small enough to find geometry inconsistencies, but does not trip on trapezoid weirdness in IDEA ultra light in O(2 um)
        if (
          !(( std::abs(solidLength_0-surfaceLength_u) < tolerance) && (std::abs(solidLength_1-surfaceLength_v) < tolerance)) &&
          !((std::abs(solidLength_1-surfaceLength_u) < tolerance) && (std::abs(solidLength_0-surfaceLength_v) < tolerance))
        )
        {
          throw GaudiException("Sensor solid (active volume) dimensions " + std::to_string(solidLength_0) + "mm and " + std::to_string(solidLength_1) + "mm could not be matched to sensitive surface dimensions (" + std::to_string(surfaceLength_u) + " x " + std::to_string(surfaceLength_v) + ") mm (including inactive material) within a tolerance of " + std::to_string(tolerance) + " mm. (sensor " + sensorKey + ", volumeID " + std::to_string(sensorVolumeID) + ", layer " + std::to_string(layer) + ", subDetector " + m_subDetName.value() + "). This indicates an inconsistency in the geometry description.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
        }

        // FOURTH: apply the parameters we found
        if (!membersDefined) {
          // For the first sensor we find, set the digitizer class members
          m_sensorActiveThickness = solidThickness;
          m_inactiveMaterialAbove = surfaceThickness_above - solidThickness/2.f;
          m_inactiveMaterialBelow = surfaceThickness_below - solidThickness/2.f;

          double pixelCountU = surfaceLength_u / m_pixelPitch[0];
          double pixelCountV = surfaceLength_v / m_pixelPitch[1];
          if (std::abs(pixelCountU - std::round(pixelCountU)) > 0.0001 || std::abs(pixelCountV - std::round(pixelCountV)) > 0.0001)
            throw GaudiException("Sensor side length (" + std::to_string(surfaceLength_u) + " x " + std::to_string(surfaceLength_v) + ") mm and pixel pitch (" + std::to_string(m_pixelPitch[0]) + " x " + std::to_string(m_pixelPitch[1]) + ") mm result in a non-integer pixel count (" + std::to_string(pixelCountU) + " x " + std::to_string(pixelCountV) + ") in subDetector " + m_subDetName.value() + ".", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
          m_pixelCount[0] = std::round(pixelCountU);
          m_pixelCount[1] = std::round(pixelCountV);
          membersDefined = true;
          debug() << "     - Setting expected dimensions according to first found sensor. Key: " << sensorKey << ", volumeID: " << sensorVolumeID << " (sensor " << sensorNumber << " in layer " << layer << ")" << endmsg;
          debug() << "       - Solid:  thickness " << solidThickness << "mm,  lengths " << solidLength_0 << " / " << solidLength_1 << "mm (arbitray ordering)" << endmsg;
          debug() << "       - Sensitive surface:  thickness " << surfaceThickness_above << "mm / " << surfaceThickness_below << "mm (above/below 0),  lengths " << surfaceLength_u << "mm / " << surfaceLength_v << "mm (u/v)" << endmsg;
        }
        else {
          // For all other sensors: check for consistency with first sensor
          if (std::abs(surfaceLength_u - m_pixelCount[0]*m_pixelPitch[0]) > 0.001 || std::abs(surfaceLength_v - m_pixelCount[1]*m_pixelPitch[1]) > 0.001 || std::abs(solidThickness - m_sensorActiveThickness) > 0.001)
            throw GaudiException("Sensor dimension mismatch found in sensor " + sensorKey + " (volumeID " + std::to_string(sensorVolumeID) + ") in layer " + std::to_string(layer) + " of subDetector " + m_subDetName.value() + ": expected dimensions of (" + std::to_string(m_pixelCount[0]*m_pixelPitch[0]) + " x " + std::to_string(m_pixelCount[1]*m_pixelPitch[1]) + " x " + std::to_string(m_sensorActiveThickness) + ") mm3, but found (" + std::to_string(surfaceLength_u) + " x " + std::to_string(surfaceLength_v) + " x " + std::to_string(solidThickness) + ") mm3. This algorithm expects exactly one type of sensor per subDetector. Use different instances of the algorithm if different layers consist of different sensors.", "VTXdigi_Modular::InitLayersAndSensors()", StatusCode::FAILURE);
        }
      } // loop over sensors
    } // loop over modules
  } // loop over layers

  info() << " - Retrieved sensor parameters: area (" << m_pixelCount[0]*m_pixelPitch[0] << " x " << m_pixelCount[1]*m_pixelPitch[1] << ") mm2, thickness " << m_sensorActiveThickness << " mm with (" << m_inactiveMaterialAbove << "/" << m_inactiveMaterialBelow << ") mm of inactive material above/below, pixel pitch (" << m_pixelPitch[0] << " x " << m_pixelPitch[1] << ") mm, pixel count (" << m_pixelCount[0] << " x " << m_pixelCount[1] << "). All " << sensorNumber << " sensors in the relevant layers share these parameters." << endmsg;
} // InitLayersAndSensors()

void VTXdigi_Modular::InitHistograms() {
  /* Define axes globally to make adjusting them easier
  * TODO: Make some of these adjustable via Gaudi Parameters? Might not be necessary.*/
  Gaudi::Accumulators::Axis<float> axis_xy{2000, -100, 100}; // global x/y in mm
  Gaudi::Accumulators::Axis<float> axis_z{4000, -200, 200};
  Gaudi::Accumulators::Axis<float> axis_cosTheta{100, 0, 1};
  Gaudi::Accumulators::Axis<float> axis_theta{4*180, 0, 180};
  Gaudi::Accumulators::Axis<float> axis_phi{4*180, -180, 180};

  Gaudi::Accumulators::Axis<float> axis_uv{1000, -20, 20};
  Gaudi::Accumulators::Axis<float> axis_w{200, -0.1, 0.1};

  Gaudi::Accumulators::Axis<float> axis_MomentumFraction{400, -1.0001f, 1.0001f};

  Gaudi::Accumulators::Axis<float> axis_moduleID{2000, -0.5f, 1999.5f};
  Gaudi::Accumulators::Axis<float> axis_clusterSize{60, 0.5f, 60.5f};
  Gaudi::Accumulators::Axis<float> axis_E{1000, 0, static_cast<float>(m_sensorActiveThickness)*2000.f};
  Gaudi::Accumulators::Axis<float> axis_charge{1000, 0, static_cast<float>(m_sensorActiveThickness)*500000.f};
  Gaudi::Accumulators::Axis<float> axis_particleE{1000, 0, 10.f}; // GeV
  Gaudi::Accumulators::Axis<float> axis_momentum_keV{10000, 0.f, 1000.f};
  Gaudi::Accumulators::Axis<float> axis_momentum_MeV{10000, 0.f, 1000.f};
  Gaudi::Accumulators::Axis<float> axis_momentum_GeV{10000, 0.f, 1000.f};
  Gaudi::Accumulators::Axis<float> axis_time{1000, 0.f, 1000.f};

  Gaudi::Accumulators::Axis<float> axis_residual{1600, -200.f, 200.f};
  Gaudi::Accumulators::Axis<float> axis_residual_abs{800, 0.f, 200.f};
  Gaudi::Accumulators::Axis<float> axis_residual_pixels{2020, -50.5f, 50.5f}; // 20 in-pix bins
  Gaudi::Accumulators::Axis<float> axis_pdg{1401, -700.5f, 700.5f};
  Gaudi::Accumulators::Axis<float> axis_mcParticleLevel{3, -0.5f, 2.5f};

  Gaudi::Accumulators::Axis<float> axis_pixels_u{
    static_cast<unsigned int>(m_pixelCount[0]),
    -0.5f,
    static_cast<float>(m_pixelCount[0]+0.5)};
  Gaudi::Accumulators::Axis<float> axis_pixels_v{
    static_cast<unsigned int>(m_pixelCount[1]),
    -0.5f,
    static_cast<float>(m_pixelCount[1]+0.5)};
  Gaudi::Accumulators::Axis<float> axis_inpix_u{
    100,
    -1.f * static_cast<float>(m_pixelPitch[0])/2.f * 1000.f,
    static_cast<float>(m_pixelPitch[0])/2.f * 1000.f};
  Gaudi::Accumulators::Axis<float> axis_inpix_v{
    100,
    -1.f * static_cast<float>(m_pixelPitch[1])/2.f * 1000.f,
    static_cast<float>(m_pixelPitch[1])/2.f * 1000.f};

  Gaudi::Accumulators::Axis<float> axis_pathLength{500, 0.f, static_cast<float>(m_sensorActiveThickness)*1000.f*10.f};
  Gaudi::Accumulators::Axis<float> axis_pathTravel{1000, -static_cast<float>(m_sensorActiveThickness)*1000.f*10.f, static_cast<float>(m_sensorActiveThickness)*1000.f*10.f};

  /* Fill histograms per layer */
  for (int layer : m_layers.value()) {
    if (layer == 0) {
      axis_z = Gaudi::Accumulators::Axis<float>{100, -96.5, 96.5}; // Want to cover layer 0 of IDEA vertex det perfectly to avoid binning-edge-effects, so we use a the correct length of 185 mm
    }

    std::array< std::unique_ptr< Gaudi::Accumulators::StaticHistogram< 1, Gaudi::Accumulators::atomicity::full, float > >, hist1dArrayLen > hist1d;

    hist1d.at(hist1d_simHit_depositedEnergy).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_energyDep",
        "SimHit deposited energy - Layer " + std::to_string(layer) + ";Energy [keV];Entries",
        axis_E
      }
    );
    hist1d.at(hist1d_simHit_depositedCharge).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_chargeDep",
        "SimHit deposited charge - Layer " + std::to_string(layer) + ";Charge [e-];Entries",
        axis_charge
      }
    );
    hist1d.at(hist1d_simHit_particleMomentum_keV).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_momentum_keV",
        "SimHit particle momentum at hit position - Layer " + std::to_string(layer) + ";Momentum [keV/c];Entries",
        axis_momentum_keV
      }
    );
    hist1d.at(hist1d_simHit_particleMomentum_MeV).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_momentum_MeV",
        "SimHit particle momentum at hit position - Layer " + std::to_string(layer) + ";Momentum [MeV/c];Entries",
        axis_momentum_MeV
      }
    );
    hist1d.at(hist1d_simHit_particleMomentum_GeV).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_momentum_GeV",
        "SimHit particle momentum at hit position - Layer " + std::to_string(layer) + ";Momentum [GeV/c];Entries",
        axis_momentum_GeV
      }
    );

    hist1d.at(hist1d_simHit_u).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/u_sensor_local",
        "SimHit position in sensor local u - Layer " + std::to_string(layer) + ";u [mm];Entries",
        axis_uv
      }
    );
    hist1d.at(hist1d_simHit_v).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/v_sensor_local",
        "SimHit position in sensor local v - Layer " + std::to_string(layer) + ";v [mm];Entries",
        axis_uv
      }
    );
    hist1d.at(hist1d_simHit_w).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/w_sensor_local",
        "SimHit position in sensor local w - Layer " + std::to_string(layer) + ";w [mm];Entries",
        axis_w
      }
    );

    hist1d.at(hist1d_simHit_x).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/x",
        "Global x-position of simHits - Layer " + std::to_string(layer) + ";x [mm];Entries",
        axis_xy
      }
    );
    hist1d.at(hist1d_simHit_y).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/y",
        "Global y-position of simHits - Layer " + std::to_string(layer) + ";y [mm];Entries",
        axis_xy
      }
    );
    hist1d.at(hist1d_simHit_z).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/z",
        "Global z-position of simHits - Layer " + std::to_string(layer) + ";z [mm];Entries",
        axis_z
      }
    );
    hist1d.at(hist1d_simHit_z_causedByPrimary).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/z_causedByPrimary",
        "Global z-position of simHits from particles created in the generator (accounts secondaries where MCParticle was deleted in ddsim) - Layer " + std::to_string(layer) + ";z [mm];Entries",
        axis_z
      }
    );
    hist1d.at(hist1d_simHit_z_causedBySecondary).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/z_causedBySecondary",
        "Global z-position of simHits from particles created in simulation (accounts secondaries where MCParticle was deleted in ddsim) - Layer " + std::to_string(layer) + ";z [mm];Entries",
        axis_z
      }
    );

    hist1d.at(hist1d_simHit_r).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/r",
        "Global simHit r position - Layer " + std::to_string(layer) + ";r [mm];Entries",
        axis_xy
      }
    );
    hist1d.at(hist1d_simHit_phi).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/phi",
        "Global simHit phi position - Layer " + std::to_string(layer) + ";phi [rad];Entries",
        axis_phi
      }
    );
    hist1d.at(hist1d_simHit_theta).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/theta",
        "Global simHit theta position - Layer " + std::to_string(layer) + ";theta [rad];Entries",
        axis_theta
      }
    );

    hist1d.at(hist1d_simHit_vertex_x).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/vertex_x",
        "X-position of the production vertex of simHits - Layer " + std::to_string(layer) + ";x [mm];Entries",
        axis_xy
      }
    );
    hist1d.at(hist1d_simHit_vertex_y).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/vertex_y",
        "Y-position of the production vertex of simHits - Layer " + std::to_string(layer) + ";y [mm];Entries",
        axis_xy
      }
    );
    hist1d.at(hist1d_simHit_vertex_z).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pos/vertex_z",
        "Z-position of the production vertex of simHits - Layer " + std::to_string(layer) + ";z [mm];Entries",
        axis_z
      }
    );


    hist1d.at(hist1d_simhit_particleMomentumDirection_x).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_momentumFraction_x",
        "X-component of the normalised momentum of simHits (at simHit position) - Layer " + std::to_string(layer) + ";Fraction of momentum in x direction [a.u.];Entries",
        axis_MomentumFraction
      }
    );
    hist1d.at(hist1d_simhit_particleMomentumDirection_y).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_momentumFraction_y",
        "Y-component of the normalised momentum of particles causing simHits (at simHit position) - Layer " + std::to_string(layer) + ";Fraction of momentum in y direction [a.u.];Entries",
        axis_MomentumFraction
      }
    );
    hist1d.at(hist1d_simhit_particleMomentumDirection_z).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_momentumFraction_z",
        "Z-component of the normalised momentum of particles causing simHits (at simHit position) - Layer " + std::to_string(layer) + ";Fraction of momentum in z direction [a.u.];Entries",
        axis_MomentumFraction
      }
    );

    hist1d.at(hist1d_simHit_particleMomentumInitialDirection_x).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_initialMomentumFraction_x",
        "X-component of the initial momentum of particles causing simHits (at their creation vertex) - Layer " + std::to_string(layer) + ";Fraction of momentum in x direction [a.u.];Entries",
        axis_MomentumFraction
      }
    );
    hist1d.at(hist1d_simHit_particleMomentumInitialDirection_y).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_initialMomentumFraction_y",
        "Y-component of the initial momentum of particles causing simHits (at their creation vertex) - Layer " + std::to_string(layer) + ";Fraction of momentum in y direction [a.u.];Entries",
        axis_MomentumFraction
      }
    );
    hist1d.at(hist1d_simHit_particleMomentumInitialDirection_z).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_initialMomentumFraction_z",
        "Z-component of the initial momentum of particles causing simHits (at their creation vertex) - Layer " + std::to_string(layer) + ";Fraction of momentum in z direction [a.u.];Entries",
        axis_MomentumFraction
      }
    );





    hist1d.at(hist1d_simHit_pdg).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pdg",
        "PDG codes of particles causing simHits - Layer " + std::to_string(layer) + ";PDG;Entries",
        axis_pdg
      }
    );

    hist1d.at(hist1d_simHit_mcParticleLevel).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_mcParticleLevel",
        "SimHits created by give MCParticle level (0 - primary, 1 - secondary, 2 - delta ray) - Layer " + std::to_string(layer) + ";MC particle level;Entries",
        axis_mcParticleLevel
      }
    );


    hist1d.at(hist1d_digiHit_collectedCharge).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_collectedCharge",
        "Charge collected per digiHit - Layer " + std::to_string(layer) + ";Charge [e-];Entries",
        axis_charge
      }
    );
    hist1d.at(hist1d_digiHit_collectedCharge_seedPixel).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_collectedCharge_seedPixel",
        "Charge collected in seed pixel (pixel with the highest charge per cluster) - Layer " + std::to_string(layer) + ";Charge [e-];Entries",
        axis_charge
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
    hist1d.at(hist1d_clusterSize_causedByPrimary).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_causedByPrimary",
        "Number of pixels per cluster (for simHits caused by particles created in the generator) - Layer " + std::to_string(layer) + ";Cluster size [pixels];Entries",
        axis_clusterSize
      }
    );
    hist1d.at(hist1d_clusterSize_causedBySecondary).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_causedBySecondary",
        "Number of pixels per cluster (for simHits caused by particles created in the Geant4 simulation) - Layer " + std::to_string(layer) + ";Cluster size [pixels];Entries",
        axis_clusterSize
      }
    );

    hist1d.at(hist1d_residual_u_toPrimaries).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_u_toPrimaries",
        "Residual (u_simHit - u_digiHit) wrt. all simHits from primary particles that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_u_toSecondaries).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_u_toSecondaries",
        "Residual (u_simHit - u_digiHit) wrt. all simHits from secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_u_toPrimariesSecondaries).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_u_toPrimariesSecondaries",
        "Residual (u_simHit - u_digiHit) wrt. all simHits from primary and secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_u_toPrimariesSecondariesDeltas).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_u_toPrimariesSecondariesDeltas",
        "Residual (u_simHit - u_digiHit) wrt. all simHits from primary and secondary particles (including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_u_maxEParticleOnSensor).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_u_maxEParticleOnSensor",
        "Residual (u_simHit - u_digiHit) in local u direction, to the simHit with the highest energy MCParticle on the sensor - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );

    hist1d.at(hist1d_residual_u_toPrimariesSecondaries_length1).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_u_toPrimariesSecondaries_length1",
        "Residual (u_simHit - u_digiHit), clusters with length 1 pix, in local u direction, wrt. all simHits from primary and secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_u_toPrimariesSecondaries_length2).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_u_toPrimariesSecondaries_length2",
        "Residual (u_simHit - u_digiHit), clusters with length 2 pix, in local u direction, wrt. all simHits from primary and secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_u_toPrimariesSecondaries_length3).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_u_toPrimariesSecondaries_length3",
        "Residual (u_simHit - u_digiHit), clusters with length 3 pix, in local u direction, wrt. all simHits from primary and secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_u_toPrimariesSecondaries_length4).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_u_toPrimariesSecondaries_length4",
        "Residual (u_simHit - u_digiHit), clusters with length 4 pix, in local u direction, wrt. all simHits from primary and secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_u_toPrimariesSecondaries_length5plus).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_u_toPrimariesSecondaries_length5plus",
        "Residual (u_simHit - u_digiHit), clusters with length 5+ pix, in local u direction, wrt. all simHits from primary and secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual u [um];Entries",
        axis_residual
      }
    );

    hist1d.at(hist1d_residual_v_toPrimaries).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_v_toPrimaries",
        "Residual (v_simHit - v_digiHit) wrt. all simHits from primary particles that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_v_toSecondaries).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_v_toSecondaries",
        "Residual (v_simHit - v_digiHit) wrt. all simHits from secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_v_toPrimariesSecondaries).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_v_toPrimariesSecondaries",
        "Residual (v_simHit - v_digiHit) wrt. all simHits from primary and secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_v_toPrimariesSecondariesDeltas).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_v_toPrimariesSecondariesDeltas",
        "Residual (v_simHit - v_digiHit) wrt. all simHits from primary and secondary particles (including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_v_maxEParticleOnSensor).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_v_maxEParticleOnSensor",
        "Residual (v_simHit - v_digiHit) in local v direction, to the simHit with the highest energy MCParticle on the sensor - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );

    hist1d.at(hist1d_residual_v_toPrimariesSecondaries_length1).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_v_toPrimariesSecondaries_length1",
        "Residual (v_simHit - v_digiHit), clusters with length 1 pix, in local v direction, wrt. all simHits from primary and secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_v_toPrimariesSecondaries_length2).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_v_toPrimariesSecondaries_length2",
        "Residual (v_simHit - v_digiHit), clusters with length 2 pix, in local v direction, wrt. all simHits from primary and secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_v_toPrimariesSecondaries_length3).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_v_toPrimariesSecondaries_length3",
        "Residual (v_simHit - v_digiHit), clusters with length 3 pix, in local v direction, wrt. all simHits from primary and secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_v_toPrimariesSecondaries_length4).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_v_toPrimariesSecondaries_length4",
        "Residual (v_simHit - v_digiHit), clusters with length 4 pix, in local v direction, wrt. all simHits from primary and secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
      }
    );
    hist1d.at(hist1d_residual_v_toPrimariesSecondaries_length5plus).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_v_toPrimariesSecondaries_length5plus",
        "Residual (v_simHit - v_digiHit), clusters with length 5 pix, in local v direction, wrt. all simHits from primary and secondary particles (not including delta rays) that contribute to this digiHit - Layer " + std::to_string(layer) + ";Residual v [um];Entries",
        axis_residual
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

    hist1d.at(hist1d_simHit_timeStamp).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_timeStamp",
        "SimHit time stamp - Layer " + std::to_string(layer) + ";Time [?s];Entries",
        axis_time
      }
    );
    hist1d.at(hist1d_digiHit_timeStamp).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_timeStamp",
        "Cluster time stamp - Layer " + std::to_string(layer) + ";Time [?s];Entries",
        axis_time
      }
    );


    hist1d.at(hist1d_simHitsPerDigiHit).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_numberOfSimHits",
        "Number of simHits per digiHit - Layer " + std::to_string(layer) + ";Number of SimHits;Entries",
        axis_clusterSize // not technically correct, but works
      }
    );

    hist1d.at(hist1d_highestEnergyParticleOnSensor_energy).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_highestEnergyParticleOnSensor_energy",
        "Energy of the particle with the highest energy on the sensor - Layer " + std::to_string(layer) + ";Energy [GeV];Entries",
        axis_particleE
      }
    );

    hist1d.at(hist1d_digiHit_u).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_pos/u_sensor_local",
        "DigiHit position on the sensor in u (eg. local x) - Layer " + std::to_string(layer) + ";Hit u position [mm];Entries",
        axis_uv
      }
    );
    hist1d.at(hist1d_digiHit_v).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_pos/v_sensor_local",
        "DigiHit position on the sensor in v (eg. local y) - Layer " + std::to_string(layer) + ";Hit v position [mm];Entries",
        axis_uv
      }
    );
    hist1d.at(hist1d_digiHit_w).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_pos/w_sensor_local",
        "DigiHit position on the sensor in w (eg. local z) - Layer " + std::to_string(layer) + ";Hit w position [mm];Entries",
        axis_w
      }
    );

    hist1d.at(hist1d_digiHit_x).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_pos/x",
        "DigiHit global x position - Layer " + std::to_string(layer) + ";Hit x position [mm];Entries",
        axis_xy
      }
    );
    hist1d.at(hist1d_digiHit_y).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_pos/y",
        "DigiHit global y position - Layer " + std::to_string(layer) + ";Hit y position [mm];Entries",
        axis_xy
      }
    );
    hist1d.at(hist1d_digiHit_z).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_pos/z",
        "DigiHit global z position - Layer " + std::to_string(layer) + ";Hit z position [mm];Entries",
        axis_z
      }
    );
    hist1d.at(hist1d_digiHit_r).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_pos/r",
        "DigiHit global r position - Layer " + std::to_string(layer) + ";Hit r position [mm];Entries",
        axis_xy
      }
    );
    hist1d.at(hist1d_digiHit_phi).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_pos/phi",
        "DigiHit global phi position - Layer " + std::to_string(layer) + ";Hit phi position [rad];Entries",
        axis_phi
      }
    );
    hist1d.at(hist1d_digiHit_theta).reset(
      new Gaudi::Accumulators::StaticHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_pos/theta",
        "DigiHit global theta position - Layer " + std::to_string(layer) + ";Hit theta position [rad];Entries",
        axis_theta
      }
    );





    m_hist1d.emplace(layer, std::move(hist1d));

    /* -- 1d-profile-hist -- */

    std::array< std::unique_ptr< Gaudi::Accumulators::StaticProfileHistogram<1,Gaudi::Accumulators::atomicity::full,float>>, histProfile1dArrayLen> histProfile1d;

    histProfile1d.at(histProfile1d_digiHit_collectedCharge_vs_global_z).reset(
      new Gaudi::Accumulators::StaticProfileHistogram<1, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_collectedCharge_vs_global_z",
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

    m_histProfile1d.emplace(layer, std::move(histProfile1d));

    /* -- 2d-hist -- */

    std::array< std::unique_ptr< Gaudi::Accumulators::StaticHistogram< 2, Gaudi::Accumulators::atomicity::full, float > >, hist2dArrayLen>  hist2d;

    hist2d.at(hist2d_digiHit_collectedCharge_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_collectedCharge_vs_global_z_2D",
        "DigiHit collected charge - Layer " + std::to_string(layer) + ";Global z position [mm];Collected charge [e-]",
        axis_z,
        axis_charge
      }
    );

    hist2d.at(hist2d_hitMap_simHits).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_hitMap",
        "Map of SimHits (pixel in which this simHit position lies) - Layer " + std::to_string(layer) + ";Pixel u;Pixel v;Entries",
        axis_pixels_u,
        axis_pixels_v
      }
    );
    hist2d.at(hist2d_hitMap_simHits_causedByPrimary).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_hitMap_causedByPrimary",
        "Map of SimHits from MCParticles created in the generator - Layer " + std::to_string(layer) + ";Pixel u;Pixel v;Entries",
        axis_pixels_u,
        axis_pixels_v
      }
    );
    hist2d.at(hist2d_hitMap_simHits_causedBySecondary).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_hitMap_causedBySecondary",
        "Map of SimHits from MCPartices created in simulation - Layer " + std::to_string(layer) + ";Pixel u;Pixel v;Entries",
        axis_pixels_u,
        axis_pixels_v
      }
    );
    hist2d.at(hist2d_hitMap_pixelHits).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_hitMap",
        "Map of pixel hits (pixels that collect charge over the threshold) - Layer " + std::to_string(layer) + ";Pixel u;Pixel v;Entries",
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
    hist2d.at(hist2d_clusterSize_vs_global_z_causedByPrimary).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_vs_global_z_causedByPrimary_2D",
        "Cluster size vs. global hit z (from MCParticles created in the generator) - Layer " + std::to_string(layer) + ";Global z position [mm];Pixels per cluster;Entries",
        axis_z,
        axis_clusterSize
      }
    );
    hist2d.at(hist2d_clusterSize_vs_global_z_causedBySecondary).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_clusterSize_vs_global_z_causedBySecondary_2D",
        "Cluster size vs. global hit z (from MCParticles created in simulation) - Layer " + std::to_string(layer) + ";Global z position [mm];Pixels per cluster;Entries",
        axis_z,
        axis_clusterSize
      }
    );

    hist2d.at(hist2d_residual_u_toPrimariesSecondaries_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_u_vs_global_z_2D",
        "Residual u (u_digiHit - u_simHit) vs. global z position - Layer " + std::to_string(layer) + ";Global z position [mm];Residual u [um];Entries",
        axis_z,
        axis_residual
      }
    );
    hist2d.at(hist2d_residual_v_toPrimariesSecondaries_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_v_vs_global_z_2D",
        "Residual v (v_digiHit - v_simHit) vs. global z position - Layer "+ std::to_string(layer) + ";Global z position [mm];Residual v [um];Entries",
        axis_z,
        axis_residual
      }
    );

    hist2d.at(hist2d_residual_u_toPrimariesSecondaries_vs_clusterPosUncertainty).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_vs_clusterPosUncertainty_u_2D",
        "Residual r (|r_digiHit - r_simHit|) vs. cluster position uncertainty in u direction - Layer " + std::to_string(layer) + ";Cluster position uncertainty u [um];Residual u [um];Entries",
        axis_residual_abs,
        axis_residual_abs
      }
    );
    hist2d.at(hist2d_residual_v_toPrimariesSecondaries_vs_clusterPosUncertainty).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/digiHit_residuals/residual_vs_clusterPosUncertainty_v_2D",
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
    hist2d.at(hist2d_pathTravel_v_vs_global_z).reset(
      new Gaudi::Accumulators::StaticHistogram<2, Gaudi::Accumulators::atomicity::full, float> {this,
        "Layer" + std::to_string(layer) + "/simHit_pathTravel_v_vs_global_z_2D",
        "Path travel length inside the sensor active volume in v direction vs. global z position - Layer " + std::to_string(layer) + ";Global z position [mm];Path length in v [um];Entries",
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
    hist2d.at(hist2d_pathTravel_vs_global_z).reset(
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

  /* global (eg not layer-specific) histograms */

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

bool VTXdigi_Modular::CheckSimhitLayer(const edm4hep::SimTrackerHit& simTrackerHit) const {
  ++m_counter_simHitsRead;
  const int layer = m_cellIdDecoder->get(simTrackerHit.getCellID(), "layer");
  if (m_layers.value().size()>0) {
    if (std::find(m_layers.value().begin(), m_layers.value().end(),  layer) == m_layers.value().end()) {
      ++m_counter_simHitsRejected_LayerNotToBeDigitized;
      return false;
    }
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

void VTXdigi_Modular::CreateDigiHits(edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitLinks, const dd4hep::DDSegmentation::VolumeID& volumeID, const TGeoHMatrix& trafoMatrix, const std::vector<VTXdigi_tools::Cluster>& clusters, TRandom3& randomGen) const {

  // Called for each cluster on a sensor, so volumeID and direction vectors are same among these clusters
  const dd4hep::rec::Vector3D direction_u_3d = VTXdigi_tools::Trafo_local_global(dd4hep::rec::Vector3D(1, 0, 0), trafoMatrix);
  const edm4hep::Vector2f direction_u = edm4hep::Vector2f(direction_u_3d.theta(), direction_u_3d.phi()); // edm4hep expects direction in to given as (theta, phi)

  const dd4hep::rec::Vector3D direction_v_3d = VTXdigi_tools::Trafo_local_global(dd4hep::rec::Vector3D(0, 1, 0), trafoMatrix);
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
    /* TODO: do we want cellID instead of VolumeID (cellID also encodes pixel segmentation instead of only the volume (ie. sensor) that the hit is in) */
    digiHit.setCellID(volumeID);
    digiHit.setEDep(cluster.charge / VTXdigi_tools::kChargePerkeV);

    // position
    VTXdigi_tools::PixelCoords clusterPos_pixC = cluster.ComputeCoG(m_clusterizeEndPixelsOnly.value());

    double cluster_pos_w;
    if (m_forceClusterPosToSensitiveSurface.value())
      cluster_pos_w = 0.f;
    else
      cluster_pos_w = m_chargeCollector->GetChargeCollectionDepthCenter();

    const dd4hep::rec::Vector3D clusterPos_local = VTXdigi_tools::Trafo_pixCoords_local(clusterPos_pixC, m_pixelPitch, m_pixelCount, cluster_pos_w);
    const dd4hep::rec::Vector3D clusterPos_global = VTXdigi_tools::Trafo_local_global(clusterPos_local, trafoMatrix);
    debug() << "     - Found cluster with " << cluster.pixels.size() << " pixels, charge " << cluster.charge << ", center at (" << clusterPos_pixC[0] << ", " << clusterPos_pixC[1] << "). Has " << cluster.simHits.size() << " contributing simHits." << endmsg;
    digiHit.setPosition(VTXdigi_tools::ConvertVector(clusterPos_global));

    // pos uncertainty
    digiHit.setU(direction_u);
    digiHit.setV(direction_v);
    if (m_positionUncertainty.value().empty()) {
      digiHit.setDu(m_pixelPitch[0] / std::sqrt(12));
      digiHit.setDv(m_pixelPitch[1] / std::sqrt(12));
    }
    else if (m_positionUncertainty.value().size() == 2) {
      digiHit.setDu(m_positionUncertainty.value().at(0));
      digiHit.setDv(m_positionUncertainty.value().at(1));
    }
    else if (m_positionUncertainty.value().size() == 10) {
      int clusterSize_u = cluster.GetSize(0);
      for (int size = 1; size < 5; size++) {
        if (clusterSize_u == size)
          digiHit.setDu(m_positionUncertainty.value().at(size - 1));
      }
      if (clusterSize_u >= 5)
        digiHit.setDu(m_positionUncertainty.value().at(4));

      int clusterSize_v = cluster.GetSize(1);
      for (int size = 1; size < 5; size++) {
        if (clusterSize_v == size)
          digiHit.setDv(m_positionUncertainty.value().at(size + 4));
      }
      if (clusterSize_v >= 5)
        digiHit.setDv(m_positionUncertainty.value().at(9));
    }
    debug() << "         - Set digiHit position uncertainty to (" << digiHit.getDu() << ", " << digiHit.getDv() << ") mm." << endmsg;

    // collect timestamp
    const VTXdigi_tools::Pixel* seedPixel = cluster.pixels.front();
    for (const VTXdigi_tools::Pixel* pix : cluster.pixels) {
      if (pix->charge > seedPixel->charge)
        seedPixel = pix; // use timestamp of pixel with the highest charge
    }
    float timeStamp = std::numeric_limits<float>::max();
    for (const auto& simHit : seedPixel->simHits) {
      timeStamp = std::min(timeStamp, simHit->hitPtr()->getTime()); // use earliest timestamp among simHits contributing to that pixel
    }
    if (m_smearing_time > 0.f)
      timeStamp += randomGen.Gaus(0., m_smearing_time.value());
    digiHit.setTime(timeStamp);

    // Create links to simHits
    for (const auto& simHit : cluster.simHits) {
      auto link = digiHitLinks.create();
      link.setFrom(digiHit);
      link.setTo(*(simHit->hitPtr()));
    }

    if (m_debugHistograms.value()) {
      FillHistograms_perDigiHit(cluster, digiHit, trafoMatrix);

      for (const VTXdigi_tools::Pixel* pix : cluster.pixels)
        FillHistograms_perPixel(volumeID, *pix);
    }
  } /* loop over clusters */
}

/* ---- Histogramming functions ---- */

void VTXdigi_Modular::FillHistograms_perSimHit(const VTXdigi_tools::SimHitWrapper& simHit) const {
  /* executed once for each simHit, no cuts applied before */

  const int layer = simHit.layer();
  const TGeoHMatrix trafoMatrix = VTXdigi_tools::ComputeSensorTrafoMatrix(simHit.volumeID(), m_volumeManager, m_sensorNormalRotation);

  const dd4hep::rec::Vector3D simHitPos_global = VTXdigi_tools::ConvertVector(simHit.hitPtr()->getPosition());
  const dd4hep::rec::Vector3D simHitPos_local = VTXdigi_tools::Trafo_global_local(simHitPos_global, trafoMatrix);
  const dd4hep::rec::Vector3D simHitProdVertex_global = VTXdigi_tools::ConvertVector(simHit.hitPtr()->getParticle().getVertex());
  const VTXdigi_tools::PixelIndex pixI = VTXdigi_tools::Trafo_local_pixI(simHitPos_local, m_pixelPitch, m_pixelCount);

  const dd4hep::rec::Vector3D simHitMomentum = VTXdigi_tools::ConvertVector(simHit.hitPtr()->getMomentum()); // in GeV
  const dd4hep::rec::Vector3D simHitMomentumInitial = VTXdigi_tools::ConvertVector(simHit.hitPtr()->getParticle().getMomentum()); // in GeV

  ++(*m_hist1d.at(layer).at(hist1d_simHit_pdg))[simHit.hitPtr()->getParticle().getPDG()];
  ++(*m_hist1d.at(layer).at(hist1d_simHit_mcParticleLevel))[static_cast<int>(simHit.mcParticleLevel())];

  ++(*m_hist1d.at(layer).at(hist1d_simHit_depositedEnergy))[simHit.hitPtr()->getEDep() * (dd4hep::GeV / dd4hep::keV)];
  ++(*m_hist1d.at(layer).at(hist1d_simHit_depositedCharge))[simHit.charge()];
  ++(*m_hist1d.at(layer).at(hist1d_simHit_particleMomentum_keV))[simHitMomentum.r() * 1.E6];
  ++(*m_hist1d.at(layer).at(hist1d_simHit_particleMomentum_MeV))[simHitMomentum.r() * 1.E3];
  ++(*m_hist1d.at(layer).at(hist1d_simHit_particleMomentum_GeV))[simHitMomentum.r()];

  ++(*m_hist1d.at(layer).at(hist1d_simHit_timeStamp))[simHit.hitPtr()->getTime()];

  ++(*m_hist1d.at(layer).at(hist1d_simHit_u))[simHitPos_local.x()];
  ++(*m_hist1d.at(layer).at(hist1d_simHit_v))[simHitPos_local.y()];
  ++(*m_hist1d.at(layer).at(hist1d_simHit_w))[simHitPos_local.z()];

  ++(*m_hist1d.at(layer).at(hist1d_simHit_x))[simHitPos_global.x()];
  ++(*m_hist1d.at(layer).at(hist1d_simHit_y))[simHitPos_global.y()];
  ++(*m_hist1d.at(layer).at(hist1d_simHit_z))[simHitPos_global.z()];

  const float r = std::sqrt(simHitPos_global.x()*simHitPos_global.x() + simHitPos_global.y()*simHitPos_global.y());
  ++(*m_hist1d.at(layer).at(hist1d_simHit_r))[r];
  const float phi = std::atan2(simHitPos_global.y(), simHitPos_global.x());
  ++(*m_hist1d.at(layer).at(hist1d_simHit_phi))[phi];
  const float theta = std::atan2(r, simHitPos_global.z());
  ++(*m_hist1d.at(layer).at(hist1d_simHit_theta))[theta];

  ++

  ++(*m_hist2d.at(layer).at(hist2d_simHit_xy))[{simHitPos_global.x(), simHitPos_global.y()}];
  ++(*m_hist2d.at(layer).at(hist2d_simHit_xz))[{simHitPos_global.x(), simHitPos_global.z()}];
  ++(*m_hist2d.at(layer).at(hist2d_simHit_yz))[{simHitPos_global.y(), simHitPos_global.z()}];

  ++(*m_hist1d.at(layer).at(hist1d_simHit_vertex_x))[simHitProdVertex_global.x()];
  ++(*m_hist1d.at(layer).at(hist1d_simHit_vertex_y))[simHitProdVertex_global.y()];
  ++(*m_hist1d.at(layer).at(hist1d_simHit_vertex_z))[simHitProdVertex_global.z()];

  ++(*m_hist1d.at(layer).at(hist1d_simhit_particleMomentumDirection_x))[simHitMomentum.x() / simHitMomentum.r()];
  ++(*m_hist1d.at(layer).at(hist1d_simhit_particleMomentumDirection_y))[simHitMomentum.y() / simHitMomentum.r()];
  ++(*m_hist1d.at(layer).at(hist1d_simhit_particleMomentumDirection_z))[simHitMomentum.z() / simHitMomentum.r()];

  ++(*m_hist1d.at(layer).at(hist1d_simHit_particleMomentumInitialDirection_x))[simHitMomentumInitial.x() / simHitMomentumInitial.r()];
  ++(*m_hist1d.at(layer).at(hist1d_simHit_particleMomentumInitialDirection_y))[simHitMomentumInitial.y() / simHitMomentumInitial.r()];
  ++(*m_hist1d.at(layer).at(hist1d_simHit_particleMomentumInitialDirection_z))[simHitMomentumInitial.z() / simHitMomentumInitial.r()];

  ++(*m_hist2d.at(layer).at(hist2d_hitMap_simHits))[{pixI[0], pixI[1]}];

  const VTXdigi_tools::MCParticleLevel mcParticleLevel = simHit.mcParticleLevel();

  if ( mcParticleLevel == VTXdigi_tools::MCParticleLevel::Primary ) {
    ++(*m_hist1d.at(layer).at(hist1d_simHit_z_causedByPrimary))[simHitPos_global.z()];

    ++(*m_hist2d.at(layer).at(hist2d_hitMap_simHits_causedByPrimary))[{pixI[0], pixI[1]}];
  }
  else if ( mcParticleLevel == VTXdigi_tools::MCParticleLevel::Secondary ) {
    ++(*m_hist1d.at(layer).at(hist1d_simHit_z_causedBySecondary))[simHitPos_global.z()];

    ++(*m_hist2d.at(layer).at(hist2d_hitMap_simHits_causedBySecondary))[{pixI[0], pixI[1]}];
  }
}

void VTXdigi_Modular::FillHistograms_perPixel(const dd4hep::DDSegmentation::VolumeID& volumeID, const VTXdigi_tools::Pixel& pix) const {
  /* executed once for each pixel */

  const int layer = GetLayer(volumeID);
  const std::array<int, 2> i_uv = pix.index;

  ++(*m_hist2d.at(layer).at(hist2d_hitMap_pixelHits))[{i_uv[0], i_uv[1]}];
}

void VTXdigi_Modular::FillHistograms_perDigiHit(const VTXdigi_tools::Cluster& cluster, const edm4hep::TrackerHitPlane& digiHit, const TGeoHMatrix& trafoMatrix) const {
  /* executed once for each digiHit */
  const int layer = GetLayer(digiHit.getCellID());
  const dd4hep::rec::Vector3D pos_global = VTXdigi_tools::ConvertVector(digiHit.getPosition());
  const dd4hep::rec::Vector3D pos_local = VTXdigi_tools::Trafo_global_local(pos_global, trafoMatrix);

  ++(*m_hist1d.at(layer).at(hist1d_digiHit_u))[ pos_local.x() ];
  ++(*m_hist1d.at(layer).at(hist1d_digiHit_v))[ pos_local.y() ];
  ++(*m_hist1d.at(layer).at(hist1d_digiHit_w))[ pos_local.z() ];

  ++(*m_hist1d.at(layer).at(hist1d_digiHit_x))[ pos_global.x() ];
  ++(*m_hist1d.at(layer).at(hist1d_digiHit_y))[ pos_global.y() ];
  ++(*m_hist1d.at(layer).at(hist1d_digiHit_z))[ pos_global.z() ];

  const float r = std::sqrt(pos_global.x()*pos_global.x() + pos_global.y()*pos_global.y());
  ++(*m_hist1d.at(layer).at(hist1d_digiHit_r))[ r ];
  const float phi = std::atan2(pos_global.y(), pos_global.x());
  ++(*m_hist1d.at(layer).at(hist1d_digiHit_phi))[ phi ];
  const float theta = std::atan2(r, pos_global.z());
  ++(*m_hist1d.at(layer).at(hist1d_digiHit_theta))[ theta ];

  ++(*m_hist1d.at(layer).at(hist1d_digiHit_collectedCharge))[ digiHit.getEDep() * VTXdigi_tools::kChargePerkeV ];
  ++(*m_hist1d.at(layer).at(hist1d_digiHit_collectedCharge_seedPixel))[ cluster.GetSeedPixelCharge() ];
  (*m_histProfile1d.at(layer).at(histProfile1d_digiHit_collectedCharge_vs_global_z))[ pos_global.z() ] += digiHit.getEDep() * VTXdigi_tools::kChargePerkeV;
  ++(*m_hist2d.at(layer).at(hist2d_digiHit_collectedCharge_vs_global_z))[ {pos_global.z(), digiHit.getEDep() * VTXdigi_tools::kChargePerkeV} ];

  ++(*m_hist1d.at(layer).at(hist1d_clusterSize))[ cluster.GetSize() ];
  ++(*m_hist1d.at(layer).at(hist1d_clusterSize_u))[ cluster.GetSize(0) ];
  ++(*m_hist1d.at(layer).at(hist1d_clusterSize_v))[ cluster.GetSize(1) ];

  (*m_histProfile1d.at(layer).at(histProfile1d_clusterSize_vs_global_z))[ pos_global.z() ] += cluster.GetSize();
  (*m_histProfile1d.at(layer).at(histProfile1d_clusterSize_u_vs_global_z))[ pos_global.z() ] += cluster.GetSize(0);
  (*m_histProfile1d.at(layer).at(histProfile1d_clusterSize_v_vs_global_z))[ pos_global.z() ] += cluster.GetSize(1);

  ++(*m_hist2d.at(layer).at(hist2d_clusterSize_vs_global_z))[ {pos_global.z(), cluster.GetSize()} ];
  ++(*m_hist2d.at(layer).at(hist2d_clusterSize_u_vs_global_z))[ {pos_global.z(), cluster.GetSize(0)} ];
  ++(*m_hist2d.at(layer).at(hist2d_clusterSize_v_vs_global_z))[ {pos_global.z(), cluster.GetSize(1)} ];

  ++(*m_hist1d.at(layer).at(hist1d_clusterPosUncertainty_u))[ digiHit.getDu() * 1000.f ]; // convert to um
  ++(*m_hist1d.at(layer).at(hist1d_clusterPosUncertainty_v))[ digiHit.getDv() * 1000.f ];

  ++(*m_hist1d.at(layer).at(hist1d_digiHit_timeStamp))[ digiHit.getTime() ];

  ++(*m_hist1d.at(layer).at(hist1d_simHitsPerDigiHit))[ cluster.simHits.size() ];

  // per simhit that is contributing to this digiHit (cluster)
  for (const auto& simHit : cluster.simHits) {
    const edm4hep::MCParticle mcParticle = simHit->hitPtr()->getParticle();

    const dd4hep::rec::Vector3D simHitPos_local = simHit->truthPos();
    const dd4hep::rec::Vector3D simHitPos_global = VTXdigi_tools::Trafo_local_global(simHitPos_local, trafoMatrix);
    const dd4hep::rec::Vector3D residual_local = simHitPos_local - pos_local; // residual = predicted - observed

    const float hit_z = simHitPos_global.z();

    ++(*m_hist1d.at(layer).at(hist1d_residual_u_toPrimariesSecondariesDeltas))[ residual_local.x()*1000.f ];
    ++(*m_hist1d.at(layer).at(hist1d_residual_v_toPrimariesSecondariesDeltas))[ residual_local.y()*1000.f ];

    const VTXdigi_tools::MCParticleLevel mcParticleLevel = simHit->mcParticleLevel();

    if (mcParticleLevel == VTXdigi_tools::MCParticleLevel::Primary) {
      ++(*m_hist1d.at(layer).at(hist1d_clusterSize_causedByPrimary))[ cluster.GetSize() ];
      ++(*m_hist2d.at(layer).at(hist2d_clusterSize_vs_global_z_causedByPrimary))[ {hit_z, cluster.GetSize()} ];

      ++(*m_hist1d.at(layer).at(hist1d_residual_u_toPrimaries))[ residual_local.x()*1000.f ];
      ++(*m_hist1d.at(layer).at(hist1d_residual_v_toPrimaries))[ residual_local.y()*1000.f ];
    }
    else if (mcParticleLevel == VTXdigi_tools::MCParticleLevel::Secondary) {
      ++(*m_hist1d.at(layer).at(hist1d_clusterSize_causedBySecondary))[ cluster.GetSize() ];
      ++(*m_hist2d.at(layer).at(hist2d_clusterSize_vs_global_z_causedBySecondary))[ {hit_z, cluster.GetSize()} ];

      ++(*m_hist1d.at(layer).at(hist1d_residual_u_toSecondaries))[ residual_local.x()*1000.f ];
      ++(*m_hist1d.at(layer).at(hist1d_residual_v_toSecondaries))[ residual_local.y()*1000.f ];
    }

    if ( mcParticleLevel == VTXdigi_tools::MCParticleLevel::Primary || mcParticleLevel == VTXdigi_tools::MCParticleLevel::Secondary ) {
      ++(*m_hist1d.at(layer).at(hist1d_residual_u_toPrimariesSecondaries))[ residual_local.x()*1000.f ];
      ++(*m_hist1d.at(layer).at(hist1d_residual_v_toPrimariesSecondaries))[ residual_local.y()*1000.f ];

      ++(*m_hist2d.at(layer).at(hist2d_residual_u_toPrimariesSecondaries_vs_global_z))[ {hit_z, residual_local.x()*1000.f} ];
      ++(*m_hist2d.at(layer).at(hist2d_residual_v_toPrimariesSecondaries_vs_global_z))[ {hit_z, residual_local.y()*1000.f} ];

      ++(*m_hist2d.at(layer).at(hist2d_residual_u_toPrimariesSecondaries_vs_clusterPosUncertainty))[ {std::abs(digiHit.getDu()) * 1000.f, std::abs(residual_local.x())*1000.f} ]; // convert to um
      ++(*m_hist2d.at(layer).at(hist2d_residual_v_toPrimariesSecondaries_vs_clusterPosUncertainty))[ {std::abs(digiHit.getDv()) * 1000.f, std::abs(residual_local.y())*1000.f} ];

      if (cluster.GetSize(0) == 1)
        ++(*m_hist1d.at(layer).at(hist1d_residual_u_toPrimariesSecondaries_length1))[ residual_local.x()*1000.f ];
      else if (cluster.GetSize(0) == 2)
        ++(*m_hist1d.at(layer).at(hist1d_residual_u_toPrimariesSecondaries_length2))[ residual_local.x()*1000.f ];
      else if (cluster.GetSize(0) == 3)
        ++(*m_hist1d.at(layer).at(hist1d_residual_u_toPrimariesSecondaries_length3))[ residual_local.x()*1000.f ];
      else if (cluster.GetSize(0) == 4)
        ++(*m_hist1d.at(layer).at(hist1d_residual_u_toPrimariesSecondaries_length4))[ residual_local.x()*1000.f ];
      else if (cluster.GetSize(0) >= 5)
        ++(*m_hist1d.at(layer).at(hist1d_residual_u_toPrimariesSecondaries_length5plus))[ residual_local.x()*1000.f ];

      if (cluster.GetSize(1) == 1)
        ++(*m_hist1d.at(layer).at(hist1d_residual_v_toPrimariesSecondaries_length1))[ residual_local.x()*1000.f ];
      else if (cluster.GetSize(1) == 2)
        ++(*m_hist1d.at(layer).at(hist1d_residual_v_toPrimariesSecondaries_length2))[ residual_local.x()*1000.f ];
      else if (cluster.GetSize(1) == 3)
        ++(*m_hist1d.at(layer).at(hist1d_residual_v_toPrimariesSecondaries_length3))[ residual_local.x()*1000.f ];
      else if (cluster.GetSize(1) == 4)
        ++(*m_hist1d.at(layer).at(hist1d_residual_v_toPrimariesSecondaries_length4))[ residual_local.x()*1000.f ];
      else if (cluster.GetSize(1) >= 5)
        ++(*m_hist1d.at(layer).at(hist1d_residual_v_toPrimariesSecondaries_length5plus))[ residual_local.x()*1000.f ];
    }
  } // loop over contributing simHits
}

void VTXdigi_Modular::FillHistograms_perSensor(const std::vector<VTXdigi_tools::SimHitWrapper>& simHits, const edm4hep::TrackerHitPlaneCollection& digiHits, const TGeoHMatrix& trafoMatrix, const dd4hep::DDSegmentation::VolumeID& volumeID) const {
  /* executed once for each sensor, after all clusters have been created */
  const int layer = GetLayer(volumeID);

  if (simHits.empty()) {
    debug() << " - No simHits found on this sensor." << endmsg;
    return;
  }

  // among simHits, find the one with the MCParticle with the highest energy
  // (needed for comparing to Allpix Squared residuals, when simulating single sensor with particle gun)
  // in case there are multiple simHits from the same MCParticle, any is fine
  size_t maxE_index = 0;
  float maxE = simHits.at(0).hitPtr()->getParticle().getEnergy(); // assume this is in GeV. Documentation is not clear

  for (size_t i = 1; i < simHits.size(); ++i) {
    const VTXdigi_tools::SimHitWrapper& simHit = simHits.at(i);
    const edm4hep::MCParticle mcParticle = simHit.hitPtr()->getParticle();
    if (mcParticle.getEnergy() > maxE && simHit.mcParticleLevel() == VTXdigi_tools::MCParticleLevel::Primary) {
      // we want to find the primary particle (assuming we only shot a single particle in simulation)
      maxE_index = i;
      maxE = mcParticle.getEnergy();
    }
  } // loop to find simHit with highest MCParticle energy.

  const VTXdigi_tools::SimHitWrapper& simHit = simHits.at(maxE_index);

  // extrapolate the simHit position to the depleted region depth centre (ie. the charge collection depth)
  // and compute the residuals to the digiHits
  // (necessary in case ddsim is ran with collectSingleDeposits=True)
  const dd4hep::rec::Vector3D simHit_pos_global = VTXdigi_tools::ConvertVector(simHit.hitPtr()->getPosition());
  const dd4hep::rec::Vector3D simHit_pos_local =VTXdigi_tools::Trafo_global_local(simHit_pos_global, trafoMatrix);

  // transform momentum to local coordinates
  double momentum_global[3] = {
    static_cast<double>(simHit.hitPtr()->getMomentum().x),
    static_cast<double>(simHit.hitPtr()->getMomentum().y),
    static_cast<double>(simHit.hitPtr()->getMomentum().z)
  };
  double momentum_local[3];
  trafoMatrix.MasterToLocalVect(momentum_global, momentum_local);
  dd4hep::rec::Vector3D simHit_dir_local = 1/std::abs(momentum_local[2]) * dd4hep::rec::Vector3D(momentum_local[0], momentum_local[1], momentum_local[2]); // normalised to w-component

  // const float targetDepth = m_chargeCollector->GetChargeCollectionDepthCenter();
  float targetDepth;
  if (m_LUT_shiftTruthPos.value())
    targetDepth = m_chargeCollector->GetChargeCollectionDepthCenter();
  else
    targetDepth = 0.f;
  const float t = (targetDepth - simHit_pos_local.z()) / simHit_dir_local.z();
  const dd4hep::rec::Vector3D simHit_pos_local_corr = simHit_pos_local + t * simHit_dir_local; // extrapolated position to the charge collection depth

  ++(*m_hist1d.at(layer).at(hist1d_highestEnergyParticleOnSensor_energy))[ maxE ]; // in GeV

  for (const auto& digiHit : digiHits) {
    const dd4hep::rec::Vector3D pos_global = VTXdigi_tools::ConvertVector(digiHit.getPosition());
    const dd4hep::rec::Vector3D pos_local = VTXdigi_tools::Trafo_global_local(pos_global, trafoMatrix);

    const dd4hep::rec::Vector3D residual_local = simHit_pos_local_corr - pos_local; // residual = predicted - observed

    ++(*m_hist1d.at(layer).at(hist1d_residual_u_maxEParticleOnSensor))[ residual_local.x()*1000.f ];
    ++(*m_hist1d.at(layer).at(hist1d_residual_v_maxEParticleOnSensor))[ residual_local.y()*1000.f ];
  }
}

void VTXdigi_Modular::FillHistograms_fromChargeCollector_perSimHit(const int layer, const dd4hep::rec::Vector3D& pathTravel, const float pathLength_Geant4, const dd4hep::rec::Vector3D& truthPos_local, const TGeoHMatrix& trafoMatrix) const {
  if (!m_debugHistograms.value()) return;

  const float pathLength = pathTravel.r(); // in mm
  const float factor_um_per_mm = 1000.f;
  const dd4hep::rec::Vector3D truthPos_global = VTXdigi_tools::Trafo_local_global(truthPos_local, trafoMatrix);

  ++(*m_hist1dglobal.at(hist1dglobal_pathTravel_r))[ pathLength*factor_um_per_mm ]; // convert from mm to um
  ++(*m_hist1dglobal.at(hist1dglobal_pathTravel_r_Geant4))[ pathLength_Geant4*factor_um_per_mm ];
  if (pathLength_Geant4 != 0.f) {
    ++(*m_hist1dglobal.at(hist1dglobal_pathTravel_r_ratio))[ pathLength / pathLength_Geant4 ];
  }

  ++(*m_hist1d.at(layer).at(hist1d_pathTravel_u))[ pathTravel.x()*factor_um_per_mm ];
  ++(*m_hist1d.at(layer).at(hist1d_pathTravel_v))[ pathTravel.y()*factor_um_per_mm ];
  ++(*m_hist1d.at(layer).at(hist1d_pathTravel_r))[ pathLength*factor_um_per_mm ];

  ++(*m_hist2d.at(layer).at(hist2d_pathTravel_u_vs_global_z))[ {truthPos_global.z(), pathTravel.x()*factor_um_per_mm} ];
  ++(*m_hist2d.at(layer).at(hist2d_pathTravel_v_vs_global_z))[ {truthPos_global.z(), pathTravel.y()*factor_um_per_mm} ];
  ++(*m_hist2d.at(layer).at(hist2d_pathTravel_w_vs_global_z))[ {truthPos_global.z(), pathTravel.z()*factor_um_per_mm} ];
  ++(*m_hist2d.at(layer).at(hist2d_pathTravel_vs_global_z))[ {truthPos_global.z(), pathLength*factor_um_per_mm} ];
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
  info() << " | " << std::setw(colWidths[0]) << std::left << "( accepted SimTrackerHits from particles created in generator )"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_accSimHitsFromPrimary.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "( accepted SimTrackerHits from particles created in Geant4 simulation, excl. deltas )"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_accSimHitsFromSecondary.value() << " |" << endmsg;
  info() << " | " << std::setw(colWidths[0]) << std::left << "( accepted SimTrackerHits from delta rays )"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_accSimHitsFromDelta.value() << " |" << endmsg;

  info()	<<	" | "	<<	std::setw(colWidths[0])	<<	std::left	<<	"Digi hits created"
         << " | " << std::setw(colWidths[1]) << std::right << m_counter_digiHitsCreated.value() << " |" << endmsg;
}

/* ---- Public helpers ---- */

dd4hep::DDSegmentation::CellID VTXdigi_Modular::GetCellID(const dd4hep::rec::Vector3D& pos_global) const {
  const dd4hep::Position pos = 0.1 * dd4hep::Position(pos_global.x(), pos_global.y(), pos_global.z()); // convert our natively used mm -> dd4hep's cm
  return m_cellIDPositionConverter->cellID(pos); // returns 0 if the position is outside of any sensitive volume
}

dd4hep::DDSegmentation::VolumeID VTXdigi_Modular::GetVolumeID(const dd4hep::DDSegmentation::CellID& cellID) const {
  /* - volumeID identifies the detector element volume (eg. sensor)
   * - cellID identifies a cell inside a volume (eg. pixel)
   * - we could convert cellID->volumeID by masking bits of the cellID, but this breaks if the cellID encoding changes later on
   *  -> use the volumeManager to convert correctly, only handle cellID's where absolutely necessary
   * - Note: this throws if a unknown cellID is passed. Catch & return 0 instead */
  if (cellID == 0) {
    // cellID 0 is used to indicate that the position is outside of any sensitive volume. Preserve this in volumeID
    return 0;
  }

  try{
    return m_volumeManager.lookupContext(cellID)->element.volumeID();  // throws if unknown
  }
  catch (const std::exception& e) {
    warning() << "VTXdigi_Modular::GetVolumeID(): Failed to convert cellID to volumeID, returning volumeID 0: " << e.what() << endmsg;
    return 0;
  }
}

int VTXdigi_Modular::GetLayer(const dd4hep::DDSegmentation::VolumeID& volumeID) const {
  return static_cast<int>(m_cellIdDecoder->get(volumeID, "layer"));
}
