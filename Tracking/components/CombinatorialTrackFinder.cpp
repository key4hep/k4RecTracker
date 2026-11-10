/**
 * CombinatorialTrackFinder: track finding on the planar hits of one tracking system.
 * See CombinatorialTrackFinder.h for an overview of the algorithm.
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

#include "CombinatorialTrackFinder.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/SurfaceHelper.h"
#include "DDRec/SurfaceManager.h"
#include "ExhaustiveSearch.h"
#include "k4FWCore/DataHandle.h"
#include "podio/Frame.h"
#include <iomanip>
#include <iostream>

CombinatorialTrackFinder::CombinatorialTrackFinder(const std::string& name, ISvcLocator* pSvcLocator)
    : MultiTransformer(
          name, pSvcLocator,
          {KeyValues{"InputHitCollection", {"TrackerHits"}}, KeyValues{"HeaderCollectionName", {"EventHeader"}},
           KeyValues{"InputRecoSimLinkCollection", {"TrackerHitSimTrackerHitLink"}}},
          {KeyValues{"OutputTrackCollection", {"DisplacedTracks"}}, KeyValues{"OutputHitCollection", {"TrackHits"}},
           KeyValues{"OutputRecoSimLinkCollection", {"DisplacedTrackHitSimLinks"}}}) {}

StatusCode CombinatorialTrackFinder::initialize() {
  StatusCode sc = Algorithm::initialize();
  if (!sc.isSuccess())
    return sc;

  if (!m_geoSvc.retrieve()) {
    error() << "Failed to retrieve GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  m_detector = m_geoSvc->getDetector();
  if (!m_detector) {
    error() << "Failed to access detector" << endmsg;
    return StatusCode::FAILURE;
  }

  // cellID decoder
  std::string cellIDEncodingString = m_geoSvc->constantAsString(m_encodingStringParameter);
  m_bitFieldCoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(cellIDEncodingString);

  // Check that the configured cellID field names exist in this readout.
  {
    const auto& desc = m_bitFieldCoder->fieldDescription();
    if (!readoutHasField(m_layerFieldName.value())) {
      error() << "Readout '" << m_encodingStringParameter.value() << "' (" << desc << ") has no field '"
              << m_layerFieldName.value() << "'. Set LayerFieldName to match this readout." << endmsg;
      return StatusCode::FAILURE;
    }
    m_hasRegionField = readoutHasField(m_regionFieldName.value());
    if (!m_hasRegionField) {
      warning() << "Readout '" << m_encodingStringParameter.value() << "' (" << desc << ") has no field '"
                << m_regionFieldName.value() << "'; treating every hit as barrel. Set RegionFieldName "
                << "(\"type\" for MuonSystemReadoutID, \"side\" for GlobalTrackerReadoutID)." << endmsg;
    }
    // When one readout spans several subdetectors, "system" goes into the layer key,
    // so that equally numbered layers of different subdetectors stay apart.
    // (automatic: on exactly when DetectorNames lists several detectors)
    const bool multiDet = m_detectorNames.value().size() > 1;
    m_useSystemInCompositeID = multiDet && readoutHasField("system");
    if (multiDet && !m_useSystemInCompositeID) {
      warning() << m_detectorNames.value().size() << " detectors share this readout but "
                << "'system' is NOT folded into the layer key -- layers with the same "
                << "index in different subdetectors will be merged." << endmsg;
    }

    info() << "CellID decoding: layer field '" << m_layerFieldName.value() << "', region field '"
           << m_regionFieldName.value() << (m_hasRegionField ? "'" : "' (absent, barrel-only)")
           << ", system in layer key: " << (m_useSystemInCompositeID ? "yes" : "no") << endmsg;
  }

  // Surface map: the one of the detector, or merged over several detectors
  auto surfaceMan = m_detector->extension<dd4hep::rec::SurfaceManager>();
  std::vector<std::string> detNames = m_detectorNames.value();
  if (detNames.empty()) {
    error() << "DetectorNames is empty: give the name(s) of the detector(s) to track in" << endmsg;
    return StatusCode::FAILURE;
  }

  if (detNames.size() == 1) {
    m_surfaceMap = surfaceMan->map(detNames.front());
    if (!m_surfaceMap) {
      error() << "Could not find surface map for detector: " << detNames.front() << endmsg;
      return StatusCode::FAILURE;
    }
  } else {
    m_mergedSurfaceMap.clear();
    for (const auto& name : detNames) {
      const auto* m = surfaceMan->map(name);
      if (!m) {
        error() << "Could not find surface map for detector: " << name << endmsg;
        return StatusCode::FAILURE;
      }
      m_mergedSurfaceMap.insert(m->begin(), m->end());
      info() << "  merged " << m->size() << " surfaces from '" << name << "'" << endmsg;
    }
    m_surfaceMap = &m_mergedSurfaceMap;
    info() << "Merged surface map over " << detNames.size() << " detectors: " << m_surfaceMap->size() << " surfaces"
           << endmsg;
  }

  // Magnetic field (the solenoid boundary is found from it below, by findSolenoidBoundary)
  m_field = m_detector->field();
  const dd4hep::Direction fieldAtCentre = m_field.magneticField(dd4hep::Position(0, 0, 0));
  info() << "Magnetic field at the detector centre: (" << fieldAtCentre.x() / dd4hep::tesla << ", "
         << fieldAtCentre.y() / dd4hep::tesla << ", " << fieldAtCentre.z() / dd4hep::tesla << ") T" << endmsg;

  // the settings shared by both search strategies
  m_commonConfig.minTrackHits = m_minTrackHits.value();
  m_commonConfig.minPtGeV = m_minComboPT.value();
  m_commonConfig.maxPairHitDistanceCm = m_maxPairHitDist.value();
  m_commonConfig.bFieldTesla = fieldAtCentre.z() / dd4hep::tesla;

  dd4hep::DetElement world = m_detector->world();
  dd4hep::rec::SurfaceHelper surfHelper(world);

  // Surfaces of the selected detector(s), grouped by layer
  m_surfaces.clear();
  for (const auto& [cellID, surface] : *m_surfaceMap) {
    const dd4hep::rec::Surface* concreteSurface = dynamic_cast<const dd4hep::rec::Surface*>(surface);
    if (concreteSurface) {
      m_surfaces.push_back(concreteSurface);
    }
  }

  m_surfacesByLayer = getSurfacesByLayer();

  // Layer key -> region (barrel / endcap), taken from the readout's region field.
  // Two different layers must not share a key.
  m_cidRegion.clear();
  {
    int nMismatch = 0;
    for (const auto* s : m_surfaces) {
      dd4hep::DetElement det = s->detElement();
      if (!det.isValid())
        continue;
      const uint64_t vid = det.volumeID();
      const int cid = getCompositeID(vid);
      const int reg = getTypeID(vid);
      auto [it, inserted] = m_cidRegion.emplace(cid, reg);
      if (!inserted && it->second != reg) {
        error() << "compositeID " << cid << " is produced by surfaces in different regions (" << it->second << " and "
                << reg << "): the 100*system+layer encoding collides "
                << "for this detector" << endmsg;
        return StatusCode::FAILURE;
      }
      if (inserted && reg != k4rt::regionOfCid(cid, nullptr))
        ++nMismatch;
    }
    if (nMismatch > 0)
      info() << nMismatch << " compositeID(s) whose region differs from the cid-range rule "
             << "(barrel system ID >= 10); using the readout region field" << endmsg;
  }

  // Per-layer geometry (radius / z extents, barrel or endcap from the surface normals).
  m_layerGeo.build(m_surfacesByLayer);
  info() << "Layer geometry (inside-out by closest approach to the IP):" << endmsg;
  for (const auto& d : m_layerGeo.layers()) {
    info() << "  cid " << std::setw(6) << d.compositeID << "  " << std::setw(7) << k4rt::kindName(d.kind) << "  r=["
           << std::fixed << std::setprecision(2) << d.rMin_cm << "," << d.rMax_cm << "]"
           << "  z=[" << d.zMin_cm << "," << d.zMax_cm << "] cm"
           << "  dMin=" << d.dMin_cm << "  nSurf=" << d.nSurfaces << std::defaultfloat << std::setprecision(6)
           << endmsg;
  }
  for (const auto& p : m_layerGeo.validate())
    warning() << "LayerGeometry: " << p << endmsg;

  if (m_ipStateMode.value() != "None" && m_ipStateMode.value() != "Helix" &&
      m_ipStateMode.value() != "SolenoidCrossing") {
    error() << "Unknown IPStateMode '" << m_ipStateMode.value() << "' (expected None, Helix or SolenoidCrossing)"
            << endmsg;
    return StatusCode::FAILURE;
  }
  // --- Strategy ---
  if (m_strategy.value() == "SeedAndExtend") {
    m_useSeedExtend = true;
  } else if (m_strategy.value() != "Exhaustive") {
    error() << "Unknown Strategy '" << m_strategy.value() << "' (expected Exhaustive or SeedAndExtend)" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_useSeedExtend) {
    k4rt::SeedConfig seedConfig;
    seedConfig.policy = m_seedLayerPolicy.value();
    seedConfig.maxCandPerLayer = m_maxCandPerLayer.value();
    seedConfig.shortTrackMaxHits = m_shortTrackMaxHits.value();
    seedConfig.shadowHitWindow_cm = m_shadowHitWindow.value();
    seedConfig.kinkMaxChi2Triplet = m_kinkMaxChi2Triplet.value();
    seedConfig.seedKinkMaxChi2NDF = m_seedKinkMaxChi2NDF.value();
    seedConfig.trackKinkMaxChi2NDF = m_trackKinkMaxChi2NDF.value();
    seedConfig.common = m_commonConfig; // the settings shared with the exhaustive search
    if (seedConfig.policy != "innermost" && seedConfig.policy != "outermost" && seedConfig.policy != "consecutive") {
      error() << "Unknown SeedLayerPolicy '" << seedConfig.policy << "' (expected innermost, outermost or consecutive)"
              << endmsg;
      return StatusCode::FAILURE;
    }
    {
      // material of each layer, for the multiple-scattering term of the kink chi2
      // x/X0 of a layer: integrated from halfway to its inner neighbour to halfway to its
      // outer neighbour, each side at most 150 cm deep
      const double layerMaterialHalfDepthCm = 150.0;
      m_layerGeo.measureMaterial(*m_detector, m_surfacesByLayer, 24, layerMaterialHalfDepthCm);
      auto& m = info();
      m << "SeedAndExtend: B = " << seedConfig.common.bFieldTesla << " T; layer x/X0 [%]:";
      for (const auto& d : m_layerGeo.layers())
        m << " " << d.compositeID << ":" << std::fixed << std::setprecision(2) << 100.0 * d.xOverX0;
      m << std::defaultfloat << std::setprecision(6) << endmsg;
    }
    // the seed finder uses the same circle fit and detector surfaces
    auto fit = [this](const std::vector<edm4hep::TrackerHitPlane>& hv, double& x0, double& y0, double& R, double& chi2,
                      Eigen::Matrix3d& cov, std::vector<double>& res) {
      try {
        return fitCircleNHits(hv, x0, y0, R, chi2, cov, res);
      } catch (...) {
        return false;
      }
    };
    auto plane = [this](uint64_t cellID, Eigen::Vector3d& o, Eigen::Vector3d& n) {
      const dd4hep::rec::Surface* s = findSurfaceByID(cellID);
      if (!s)
        return false;
      const auto so = s->origin();
      const auto sn = s->normal(so);
      o = Eigen::Vector3d(so.x(), so.y(), so.z());
      n = Eigen::Vector3d(sn.x(), sn.y(), sn.z()).normalized();
      return true;
    };
    m_seedFinder = std::make_unique<k4rt::SeedFinder>(m_layerGeo, seedConfig, fit, plane);
    for (const auto& p : m_seedFinder->prepare())
      warning() << "SeedFinder: " << p << endmsg;
    {
      auto& m = info();
      m << "Strategy SeedAndExtend: passes";
      for (int n : m_seedFinder->passes())
        m << " [" << n << " hits: " << m_seedFinder->numSets(n) << " layer sets]";
      m << " (policy " << seedConfig.policy << ")" << endmsg;
    }
    if (m_seedFinder->seedLayerSets().empty()) {
      error() << "SeedAndExtend: no seed layer sets - check SeedLayerPolicy" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  info() << "Found " << m_surfaces.size() << " detector surfaces in " << m_surfacesByLayer.size() << " layers"
         << endmsg;

  // Solenoid boundary, from where Bz changes sign
  m_solenoidBoundaryMm = findSolenoidBoundary(false);
  info() << "Solenoid Bz boundary (radial): R = " << std::fixed << std::setprecision(4) << m_solenoidBoundaryMm
         << " mm  (" << m_solenoidBoundaryMm / 10.0 << " cm)" << endmsg;

  m_solenoidHalfZmm = findSolenoidBoundary(true);
  info() << "Solenoid Bz boundary (endcap): |Z| = " << std::fixed << std::setprecision(4) << m_solenoidHalfZmm
         << " mm  (" << m_solenoidHalfZmm / 10.0 << " cm)" << endmsg;

  return StatusCode::SUCCESS;
}

// Operator method
std::tuple<edm4hep::TrackCollection, edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection>
CombinatorialTrackFinder::operator()(const edm4hep::TrackerHitPlaneCollection& hits,
                                     const edm4hep::EventHeaderCollection& headers,
                                     const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks) const {

  // Find tracks using the direct EDM4hep approach
  edm4hep::TrackCollection trackCollection;
  edm4hep::TrackerHitPlaneCollection outputHits;
  edm4hep::TrackerHitSimTrackerHitLinkCollection outputLinks;

  unsigned int eventNumber = 0;
  if (!headers.empty())
    eventNumber = headers[0].getEventNumber();
  debug() << "Event " << eventNumber << ": " << hits.size() << " hits" << endmsg;

  if (hits.size() < 3) {
    m_tracksPerEvent += 0;
    return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks));
  }

  try {
    findTracks(&hits, trackCollection, outputHits, recoSimLinks, outputLinks);
  } catch (const std::exception& ex) {
    error() << "Exception during track finding: " << ex.what() << endmsg;
    m_tracksPerEvent += 0;
    return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks));
  }

  m_tracksPerEvent += trackCollection.size();

  // one line per found track (DEBUG only)
  if (msgLevel(MSG::DEBUG))
    printTrackSummary(trackCollection);

  return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks));
}

// DEBUG printout after each event: one line per track, with its state at the first hit
// (the comparison with the MC truth belongs in the analysis).
void CombinatorialTrackFinder::printTrackSummary(const edm4hep::TrackCollection& trackCollection) const {
  debug() << "  " << trackCollection.size() << " track(s)" << endmsg;
  for (size_t i = 0; i < trackCollection.size(); ++i) {
    const auto& track = trackCollection[i];
    for (const auto& state : track.getTrackStates()) {
      if (state.location != edm4hep::TrackState::AtFirstHit)
        continue;
      // field at the first hit (positions in mm, DD4hep in cm)
      double bField = m_field.magneticField(dd4hep::Position(0, 0, 0)).z() / dd4hep::tesla;
      if (track.trackerHits_size() > 0) {
        const auto position = track.getTrackerHits(0).getPosition();
        bField =
            m_field.magneticField(dd4hep::Position(position[0] / 10.0, position[1] / 10.0, position[2] / 10.0)).z() /
            dd4hep::tesla;
      }
      const double pT = (std::abs(state.omega) > 1e-12) ? k4rt::ptFromOmega(bField, state.omega) : 0.0;
      debug() << "  track " << i + 1 << ": " << track.trackerHits_size() << " hits, chi2/ndf "
              << track.getChi2() / std::max(1, track.getNdf()) << ", pT " << pT << " GeV, tanLambda " << state.tanLambda
              << endmsg;
      break;
    }
  }
}

const dd4hep::rec::Surface* CombinatorialTrackFinder::findSurface(const edm4hep::TrackerHitPlane& hit) const {
  return findSurfaceByID(hit.getCellID());
}

const dd4hep::rec::Surface* CombinatorialTrackFinder::findSurfaceByID(uint64_t cellID) const {
  dd4hep::rec::SurfaceMap::const_iterator sI = m_surfaceMap->find(cellID);

  if (sI != m_surfaceMap->end()) {
    return dynamic_cast<const dd4hep::rec::Surface*>(sI->second);
  }

  return nullptr;
}

// Region of a cellID: 0 barrel, +1 / -1 endcaps (field RegionFieldName).
int CombinatorialTrackFinder::getTypeID(uint64_t cellID) const {
  if (!m_hasRegionField)
    return 0; // single-region readout: treat everything as barrel
  return m_bitFieldCoder->get(cellID, m_regionFieldName.value());
}

int CombinatorialTrackFinder::getLayerID(uint64_t cellID) const {
  return m_bitFieldCoder->get(cellID, m_layerFieldName.value());
}

// Layer key of a cellID (layout in LayerGeometry.h):
//   barrel:    100*system + layer
//   +endcap:   100*system + layer + 10000
//   -endcap: -(100*system + layer + 10000)
// The system term is only added when several subdetectors share one readout.
int CombinatorialTrackFinder::getCompositeID(uint64_t cellID) const {
  int layerID = getLayerID(cellID);
  int typeID = getTypeID(cellID);
  int base = layerID;
  if (m_useSystemInCompositeID)
    base += 100 * m_bitFieldCoder->get(cellID, "system");
  if (typeID == 0)
    return base;
  else if (typeID == 1)
    return k4rt::kEndcapKeyOffset + base;
  else
    return -k4rt::kEndcapKeyOffset - base;
}

// Detector surfaces grouped by layer key.
std::map<int, std::vector<const dd4hep::rec::Surface*>> CombinatorialTrackFinder::getSurfacesByLayer() const {
  std::map<int, std::vector<const dd4hep::rec::Surface*>> result;

  for (const auto& surface : m_surfaces) {
    dd4hep::DetElement det = surface->detElement();
    if (!det.isValid())
      continue;

    uint64_t volID = det.volumeID();

    // same key as for the hits
    const int compositeID = getCompositeID(volID);
    result[compositeID].push_back(surface);
  }

  if (msgLevel(MSG::INFO)) {
    info() << "Layer distribution:" << endmsg;
    for (const auto& [layer, surfaces] : result) {
      std::string typeStr;
      int actualLayer;

      const int region = k4rt::regionOfKey(layer);
      typeStr = (region == 0) ? "Barrel" : (region > 0 ? "Positive Endcap" : "Negative Endcap");
      actualLayer = k4rt::keyWithoutRegion(layer);

      info() << "  " << typeStr << " Layer " << actualLayer << " (ID=" << layer << "): " << surfaces.size()
             << " surfaces" << endmsg;
    }
  }

  return result;
}

// pT from omega, with Bz at the state's reference point.
double CombinatorialTrackFinder::getPT(const edm4hep::TrackState& state) const {
  double omega = state.omega; // 1/mm

  double x = state.referencePoint[0];
  double y = state.referencePoint[1];
  double z = state.referencePoint[2];

  dd4hep::Position fieldPos(x / 10.0, y / 10.0, z / 10.0); // mm -> cm
  double bField = m_field.magneticField(fieldPos).z() / dd4hep::tesla;

  return k4rt::ptFromOmega(bField, omega); // GeV/c
}

// Where along the fitted helix a stored track state is placed.
enum class StatePoint { FirstHit, LastHit, Origin };

// Describe the fitted helix as an EDM4hep track state at one point of the track.
//
// The helix is the fitted circle (centre, radius) plus a straight-line fit of z versus
// the transverse arc length. Everything is geometric, so it works for either field sign:
//   - the particle moves outward, from the hit nearest the origin to the farthest one;
//   - omega > 0 means clockwise seen from +z (EDM4hep); in a field along +z its sign is
//     the charge;
//   - the reference point is the innermost hit, the outermost hit or the origin, and
//     d0, phi, z0 are given at the circle point closest to it (EDM4hep convention:
//     pca = ref + (-d0 sin(phi), d0 cos(phi), z0), centre = pca + (1/omega)(sin(phi), -cos(phi)),
//     z(s) = pca_z + tanLambda * s).
// Units: centre and radius in cm, circleCov is the 3x3 covariance of (x0, y0, R) in cm^2;
// the state is written in mm. Returns false (state unchanged) for degenerate input.
static bool makeEdm4hepStateFromCircle(double centreXcm, double centreYcm, double radiusCm,
                                       const Eigen::Matrix3d& circleCov,
                                       const std::vector<edm4hep::TrackerHitPlane>& trackHits, StatePoint where,
                                       edm4hep::TrackState& state) {
  const double centreX = centreXcm * 10.0; // mm
  const double centreY = centreYcm * 10.0;
  const double radius = radiusCm * 10.0;
  if (radius <= 0.0 || trackHits.size() < 2)
    return false;

  // 1. Hits sorted by distance from the origin: the order an outward-moving particle
  //    crosses them.
  struct HitPoint {
    double x, y, z;   // mm
    double distance;  // 3D distance from the origin, mm
    double varianceZ; // mm^2
  };
  std::vector<HitPoint> points;
  for (const auto& hit : trackHits) {
    const auto& pos = hit.getPosition();
    HitPoint point;
    point.x = pos[0];
    point.y = pos[1];
    point.z = pos[2];
    point.distance = std::sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
    point.varianceZ = hit.getCovMatrix(edm4hep::Cartesian::z, edm4hep::Cartesian::z);
    points.push_back(point);
  }
  std::sort(points.begin(), points.end(), [](const HitPoint& a, const HitPoint& b) { return a.distance < b.distance; });

  // 2. Sense of rotation seen from +z (+1 = counter-clockwise): sum of the cross
  //    products of consecutive hits, measured from the circle centre.
  double crossSum = 0.0;
  for (size_t i = 0; i + 1 < points.size(); ++i) {
    const double ax = points[i].x - centreX;
    const double ay = points[i].y - centreY;
    const double bx = points[i + 1].x - centreX;
    const double by = points[i + 1].y - centreY;
    crossSum += ax * by - ay * bx;
  }
  if (crossSum == 0.0)
    return false;
  const double sense = (crossSum > 0.0) ? 1.0 : -1.0;

  // 3. Reference point and the circle point closest to it (pca).
  double refX = 0.0, refY = 0.0, refZ = 0.0;
  size_t nearestHit = 0; // hit next to the reference point, used to count the arc
  if (where == StatePoint::FirstHit) {
    refX = points.front().x;
    refY = points.front().y;
    refZ = points.front().z;
  } else if (where == StatePoint::LastHit) {
    nearestHit = points.size() - 1;
    refX = points.back().x;
    refY = points.back().y;
    refZ = points.back().z;
  }
  const double toRefX = refX - centreX;
  const double toRefY = refY - centreY;
  const double centreToRef = std::sqrt(toRefX * toRefX + toRefY * toRefY);
  if (centreToRef < 1e-6)
    return false;                                     // reference point at the circle centre: no closest point
  const double anglePca = std::atan2(toRefY, toRefX); // angle of the pca around the centre
  const double pcaX = centreX + radius * toRefX / centreToRef;
  const double pcaY = centreY + radius * toRefY / centreToRef;

  // 4. phi = direction of travel at the pca (radius vector turned by +-90 degrees),
  //    omega > 0 for clockwise motion, d0 from pca - ref = (-d0 sin(phi), d0 cos(phi)).
  double phi = anglePca + sense * M_PI / 2.0;
  phi = k4rt::wrapAngle(phi);
  const double omega = -sense / radius; // 1/mm
  const double d0 = -(pcaX - refX) * std::sin(phi) + (pcaY - refY) * std::cos(phi);

  // 5. Angle turned along the direction of travel, hit by hit (no 2 pi jumps), then the
  //    arc length of every hit measured from the pca.
  std::vector<double> turned(points.size());
  double previousAngle = std::atan2(points[0].y - centreY, points[0].x - centreX);
  turned[0] = 0.0;
  for (size_t i = 1; i < points.size(); ++i) {
    const double angle = std::atan2(points[i].y - centreY, points[i].x - centreX);
    double step = sense * (angle - previousAngle);
    step = k4rt::wrapAngle(step);
    turned[i] = turned[i - 1] + step;
    previousAngle = angle;
  }
  const double angleNearest = std::atan2(points[nearestHit].y - centreY, points[nearestHit].x - centreX);
  double pcaStep = sense * (anglePca - angleNearest);
  pcaStep = k4rt::wrapAngle(pcaStep);
  const double turnedPca = turned[nearestHit] + pcaStep;

  // 6. Weighted straight-line fit z = zPca + tanLambda * s, weights 1/sigma_z^2 from the hit
  //    covariance. A disk hit measures z exactly (z is along its normal), so its sigma_z^2 is 0
  //    or rounding noise (1e-34 mm^2): it gets the floor of (1 um)^2, the largest weight.
  const double minimumVarianceZ = 1e-6; // mm^2 = (1 um)^2
  std::vector<double> weight(points.size());
  std::vector<double> arcLength(points.size());
  double sumW = 0.0, sumWS = 0.0, sumWZ = 0.0;
  for (size_t i = 0; i < points.size(); ++i) {
    weight[i] = 1.0 / std::max(points[i].varianceZ, minimumVarianceZ);
    arcLength[i] = radius * (turned[i] - turnedPca);
    sumW += weight[i];
    sumWS += weight[i] * arcLength[i];
    sumWZ += weight[i] * points[i].z;
  }
  // the fit about the weighted means: no cancellation even when one weight is much larger
  const double meanS = sumWS / sumW;
  const double meanZ = sumWZ / sumW;
  double sumWDSDS = 0.0, sumWDSDZ = 0.0;
  for (size_t i = 0; i < points.size(); ++i) {
    const double ds = arcLength[i] - meanS;
    sumWDSDS += weight[i] * ds * ds;
    sumWDSDZ += weight[i] * ds * (points[i].z - meanZ);
  }
  // degenerate when the hits do not spread along the arc (weighted spread below 1 um)
  if (!(sumWDSDS / sumW > 1e-6))
    return false;
  const double tanLambda = sumWDSDZ / sumWDSDS;
  const double zPca = meanZ - tanLambda * meanS;

  // 7. Covariance. d0, phi and omega from the circle-fit covariance (gradients with
  //    respect to x0, y0, R in cm), with lower limits on the errors;
  //    z0 and tanLambda from the line fit.
  const double D2 = centreToRef * centreToRef / 100.0;                              // cm^2
  const Eigen::Vector3d gradD0(-toRefX / centreToRef, -toRefY / centreToRef, -1.0); // cm/cm
  const Eigen::Vector3d gradPhi(toRefY / 10.0 / D2, -toRefX / 10.0 / D2, 0.0);      // rad/cm
  const Eigen::Vector3d gradOmega(0.0, 0.0, 1.0 / (radiusCm * radiusCm * 10.0));    // (1/mm)/cm
  const double varD0 = std::max(gradD0.dot(circleCov * gradD0) * 100.0, 0.01);      // mm^2
  const double varPhi = std::max(gradPhi.dot(circleCov * gradPhi), 0.01);
  const double varOmega = std::max(gradOmega.dot(circleCov * gradOmega), std::pow(0.1 * omega, 2));

  using TP = edm4hep::TrackParams;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j <= i; ++j)
      state.setCovMatrix(0.0f, static_cast<TP>(i), static_cast<TP>(j));
  state.setCovMatrix(static_cast<float>(varD0), TP::d0, TP::d0);
  state.setCovMatrix(static_cast<float>(varPhi), TP::phi, TP::phi);
  state.setCovMatrix(static_cast<float>(varOmega), TP::omega, TP::omega);
  // line-fit errors: var(tanLambda) = 1 / S, var(z0) = 1 / sumW + meanS^2 / S,
  // cov(z0, tanLambda) = -meanS / S, with S = sum of w (s - meanS)^2
  state.setCovMatrix(static_cast<float>(1.0 / sumW + meanS * meanS / sumWDSDS), TP::z0, TP::z0);
  state.setCovMatrix(static_cast<float>(1.0 / sumWDSDS), TP::tanLambda, TP::tanLambda);
  state.setCovMatrix(static_cast<float>(-meanS / sumWDSDS), TP::z0, TP::tanLambda);

  state.D0 = static_cast<float>(d0);
  state.phi = static_cast<float>(phi);
  state.omega = static_cast<float>(omega);
  state.Z0 = static_cast<float>(zPca - refZ);
  state.tanLambda = static_cast<float>(tanLambda);
  state.referencePoint = edm4hep::Vector3f(refX, refY, refZ);
  return true;
}

// Track state with the given parameters and a default diagonal covariance.
edm4hep::TrackState CombinatorialTrackFinder::createTrackState(double d0, double phi, double omega, double z0,
                                                               double tanLambda, int location) const {

  edm4hep::TrackState state;

  state.D0 = d0;
  state.phi = phi;
  state.omega = omega;
  state.Z0 = z0;
  state.tanLambda = tanLambda;
  state.location = location;

  // callers that know the fit uncertainties overwrite these
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j <= i; ++j)
      state.setCovMatrix(0.0, static_cast<edm4hep::TrackParams>(i), static_cast<edm4hep::TrackParams>(j));

  // D0/Z0 in mm^2, phi/tanLambda in rad^2, omega in (1/mm)^2
  state.setCovMatrix(1.0, edm4hep::TrackParams::d0, edm4hep::TrackParams::d0);
  state.setCovMatrix(0.01, edm4hep::TrackParams::phi, edm4hep::TrackParams::phi);
  state.setCovMatrix(1e-10, edm4hep::TrackParams::omega, edm4hep::TrackParams::omega);
  state.setCovMatrix(1.0, edm4hep::TrackParams::z0, edm4hep::TrackParams::z0);
  state.setCovMatrix(0.01, edm4hep::TrackParams::tanLambda, edm4hep::TrackParams::tanLambda);

  return state;
}

//------------------------------------------------------------------------------
// Track finding
//------------------------------------------------------------------------------

namespace {
/// One group of hits searched on its own: its layers (inner to outer) and the hits of each.
struct ClusterDef {
  std::vector<int> cids;
  std::vector<std::vector<size_t>> hitsForComp;
};

/// Union-find (disjoint sets) over hit indices, used by the hit clustering: unite() joins the
/// groups of two hits, find() gives the representative of a hit's group.
struct UnionFind {
  std::vector<size_t> parent;
  explicit UnionFind(size_t n) : parent(n) {
    for (size_t i = 0; i < n; ++i)
      parent[i] = i;
  }
  size_t find(size_t x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]]; // path halving
      x = parent[x];
    }
    return x;
  }
  void unite(size_t a, size_t b) {
    a = find(a);
    b = find(b);
    if (a != b)
      parent[a] = b;
  }
};
} // namespace

// The track search in one event: one track per iteration, until no candidate is left.
// Each iteration
//   1. collects the unused hits (collectHits),
//   2. chooses a combination of hits: the next SeedAndExtend candidate, or the best combination
//      of the exhaustive search (searchCombinations),
//   3. (Exhaustive) refines it: guided hit choice in crowded layers, neighbour-track cuts,
//      outlier removal,
//   4. stores it as a track with its states, hits and reco-sim links (buildTrack).
class CombinatorialTrackFinder::EventSearch {
public:
  EventSearch(const CombinatorialTrackFinder& finder, const edm4hep::TrackerHitPlaneCollection& hits,
              edm4hep::TrackCollection& finalTracks, edm4hep::TrackerHitPlaneCollection& outputHits,
              const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks,
              edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinks);

  void run();

private:
  // ---- the steps of one iteration ----
  bool collectHits();
  bool nextSeedAndExtendCandidate(k4rt::CombResult& best);
  bool searchCombinations(k4rt::CombResult& best);
  std::vector<ClusterDef> buildClusters(const std::vector<std::vector<size_t>>& hitsForCompAll) const;
  void joinLayers(UnionFind& groups, const std::vector<std::vector<size_t>>& hitsForCompAll, int step, double phiLimit,
                  double distLimit) const;
  void searchCluster(const ClusterDef& cluster, const k4rt::RoadCuts& roadCuts, k4rt::CombResult& best,
                     k4rt::CombResult& fallback, bool& found);
  void selectCrowdedLayerHits(k4rt::CombResult& best);
  void collectCleanHits(std::vector<edm4hep::TrackerHitPlane>& cleanHits, std::vector<size_t>& cleanHIIs,
                        std::vector<size_t>& cleanGIdxs) const;
  void applyGuidedChoice(const std::vector<edm4hep::TrackerHitPlane>& guidedHits, const std::vector<size_t>& guidedHIIs,
                         const std::vector<size_t>& guidedGIdxs, k4rt::CombResult& best) const;
  bool rejectedAsNeighbour(const k4rt::CombResult& best);
  bool nearAcceptedTrack(size_t hii, double distanceSquaredMm2) const;
  void removeOutliers(k4rt::CombResult& best) const;
  void buildTrack(const k4rt::CombResult& best);

  // ---- hit tests of the exhaustive search ----
  bool isolated(const std::vector<size_t>& gIdxVec) const;
  bool coincident(const edm4hep::TrackerHitPlane& a, const edm4hep::TrackerHitPlane& b) const;
  int completeDoublets(const std::vector<size_t>& gIdxVec) const;

  // ---- reco-sim links ----
  // Key of a hit: (collectionID, index). Merged input collections all start their indices at 0.
  template <typename Hit>
  static uint64_t hitKey(const Hit& h) {
    const auto id = h.id();
    return (static_cast<uint64_t>(id.collectionID) << 32) | static_cast<uint32_t>(id.index);
  }
  void propagateLink(const edm4hep::MutableTrackerHitPlane& outHit, const edm4hep::TrackerHitPlane& inHit);

  // ---- the event ----
  const CombinatorialTrackFinder& finder;
  const edm4hep::TrackerHitPlaneCollection& hits;
  edm4hep::TrackCollection& finalTracks;
  edm4hep::TrackerHitPlaneCollection& outputHits;
  edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinks;
  /// input hit -> all SimTrackerHits that contributed to it (one digi hit can have several)
  std::unordered_map<uint64_t, std::vector<edm4hep::SimTrackerHit>> recoToSimMap;
  /// hits already used by accepted tracks
  std::vector<bool> globalUsedHits;
  /// SeedAndExtend: all candidates are built in the first iteration, then taken one per
  /// iteration, best first, skipping those that share a hit with an accepted track
  std::vector<k4rt::SECandidate> seCands;
  bool seBuilt = false;
  size_t seNext = 0;
  /// hit positions (mm) of the accepted tracks, for the neighbour-track cuts
  std::vector<std::array<double, 3>> acceptedTrackHitPositions;
  int trackNumber = 0;
  /// hard pair coincidence (RequirePairCoincidence)
  bool pcHard = false;

  // ---- the current iteration (filled by collectHits) ----
  std::vector<k4rt::HitInfo> allHitInfo;                         ///< the unused hits with a surface
  std::map<int, std::vector<size_t>> hitIndicesByCompositeLayer; ///< layer key -> allHitInfo indices
  std::vector<int> compositeLayerIDs;                            ///< layer keys, inner -> outer
  std::unordered_map<size_t, int> gIdxToComposite;               ///< hit index -> layer key
  std::unordered_map<size_t, double> hitEdepScore;               ///< allHitInfo index -> edep score (keV)
  std::vector<bool> usedHits;                                    ///< hits used before this iteration
  double absBRefForK3 = 0.0; ///< |Bz| at the centroid of the hits (T), for the pT of a combination
  double rCapK3 = 0.0;       ///< radius (cm) of a track with pT = MaxComboPT
  int maxK = 0;              ///< most hits of a combination
  int minK = 0;              ///< fewest hits of a combination
};

CombinatorialTrackFinder::EventSearch::EventSearch(const CombinatorialTrackFinder& finderIn,
                                                   const edm4hep::TrackerHitPlaneCollection& hitsIn,
                                                   edm4hep::TrackCollection& finalTracksIn,
                                                   edm4hep::TrackerHitPlaneCollection& outputHitsIn,
                                                   const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks,
                                                   edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinksIn)
    : finder(finderIn), hits(hitsIn), finalTracks(finalTracksIn), outputHits(outputHitsIn), outputLinks(outputLinksIn),
      globalUsedHits(hitsIn.size(), false), pcHard(finderIn.m_requirePairCoinc.value()) {
  for (const auto& link : recoSimLinks) {
    auto recoHit = link.getFrom();
    recoToSimMap[hitKey(recoHit)].push_back(link.getTo());
  }
}

// Given a newly created output hit whose data was copied from an input hit, every reco-sim link
// of the input hit is copied to the output link collection (cluster digitisation can give several
// SimHits per hit, e.g. a muon and a delta ray in the same readout cell).
void CombinatorialTrackFinder::EventSearch::propagateLink(const edm4hep::MutableTrackerHitPlane& outHit,
                                                          const edm4hep::TrackerHitPlane& inHit) {
  auto it = recoToSimMap.find(hitKey(inHit));
  if (it == recoToSimMap.end())
    return;
  for (const auto& sim : it->second) {
    auto outLink = outputLinks.create();
    outLink.setFrom(outHit);
    outLink.setTo(sim);
    outLink.setWeight(1.0f);
  }
}

void CombinatorialTrackFinder::EventSearch::run() {
  if (hits.size() < 3)
    return;
  while (true) {
    trackNumber++;
    if (!collectHits())
      break;

    // the combination of hits of this iteration
    k4rt::CombResult best;
    const bool found = finder.m_useSeedExtend ? nextSeedAndExtendCandidate(best) : searchCombinations(best);
    if (!found)
      break;

    if (!finder.m_useSeedExtend) {
      if (finder.m_doCrowdedLayerHitSelection)
        selectCrowdedLayerHits(best);
      if (rejectedAsNeighbour(best))
        continue;
      removeOutliers(best);
    }
    buildTrack(best);
  }
}

// The unused hits of this iteration (with a known surface, below MaxHitEdepKeV), their layers in
// inner -> outer order, the energy-deposit score of every hit and the field at their centroid.
// false: fewer than 3 hits or layers left, the search stops.
bool CombinatorialTrackFinder::EventSearch::collectHits() {
  int remainingHits = 0;
  for (size_t i = 0; i < globalUsedHits.size(); ++i) {
    if (!globalUsedHits[i])
      remainingHits++;
  }
  if (remainingHits < 3)
    return false;

  allHitInfo.clear();
  hitIndicesByCompositeLayer.clear();
  allHitInfo.reserve(remainingHits);

  // MaxHitEdepKeV: skip hits with a large energy deposit (EDM4hep edep is in GeV)
  const double edepThreshGeV = (finder.m_maxHitEdepKeV.value() > 0.0)
                                   ? finder.m_maxHitEdepKeV.value() * 1e-6 // keV -> GeV
                                   : -1.0;                                 // disabled
  for (size_t i = 0; i < hits.size(); ++i) {
    if (globalUsedHits[i])
      continue;
    const auto& hit = hits[i];
    if (edepThreshGeV > 0.0 && hit.getEDep() > edepThreshGeV)
      continue;
    if (finder.findSurface(hit)) {
      const int layerID = finder.getLayerID(hit.getCellID());
      const int typeID = finder.getTypeID(hit.getCellID());
      const int compositeID = finder.getCompositeID(hit.getCellID());
      allHitInfo.emplace_back(i, hit, layerID, typeID, compositeID);
      hitIndicesByCompositeLayer[compositeID].push_back(allHitInfo.size() - 1);
    }
  }
  if (allHitInfo.size() < 3)
    return false;

  // layer keys sorted by layer number (innermost first), then by key
  compositeLayerIDs.clear();
  for (const auto& [compositeID, indices] : hitIndicesByCompositeLayer)
    compositeLayerIDs.push_back(compositeID);
  std::sort(compositeLayerIDs.begin(), compositeLayerIDs.end(), [](int a, int b) {
    const int la = k4rt::keyWithoutRegion(a);
    const int lb = k4rt::keyWithoutRegion(b);
    if (la != lb)
      return la < lb; // smaller layer number => nearer to IP
    return a < b;     // tie-break by compositeID numeric order
  });
  if (compositeLayerIDs.size() < 3)
    return false;

  // Field at the centroid of the hits, to turn MaxComboPT into a radius cap for the
  // k = 3 score: the radius (cm) of a track with pT = MaxComboPT
  double bRefForK3 = 0.0;
  {
    double sx = 0.0, sy = 0.0, sz = 0.0;
    for (const auto& hi : allHitInfo) {
      const auto& p = hi.hit.getPosition(); // mm
      sx += p.x;
      sy += p.y;
      sz += p.z;
    }
    const double n = static_cast<double>(allHitInfo.size());
    dd4hep::Position centroid(sx / n / 10.0, sy / n / 10.0, sz / n / 10.0); // cm
    bRefForK3 = finder.m_field.magneticField(centroid).z() / dd4hep::tesla;
  }
  absBRefForK3 = std::max(std::abs(bRefForK3), 1e-6);
  rCapK3 = k4rt::radiusFromPt(absBRefForK3, finder.m_maxComboPT.value());

  // --- Energy-deposit score ---
  // In a layer with several hits, the right one usually has an energy deposit like
  // the hits in layers with a single hit (the event baseline); delta rays and
  // secondaries deposit much more or much less.
  //   hit score  = distance of edep from [baseline, 2 x baseline]  (keV, lower = better)
  //   combo score: k = 3 : w * avgScore / EdepNormKeV - min(R, R_cap) / R_cap
  //                k >= 4: chi2/ndf + w * avgScore / EdepNormKeV
  // MaxHitsPerCompositeID keeps only the best-scoring hits of crowded layers.
  std::unordered_map<size_t, double> hitEdepKeV; // allHitInfo index -> edep (keV)
  hitEdepKeV.reserve(allHitInfo.size());
  for (size_t idx = 0; idx < allHitInfo.size(); ++idx)
    hitEdepKeV[idx] = allHitInfo[idx].hit.getEDep() * 1e6; // GeV -> keV

  // baseline: mean edep of the layers with exactly one hit
  double meanEdepSingle = 0.0;
  int nSingle = 0;
  for (const auto& [cid, cidIndices] : hitIndicesByCompositeLayer) {
    if (cidIndices.size() == 1) {
      meanEdepSingle += hitEdepKeV.at(cidIndices[0]);
      ++nSingle;
    }
  }
  if (nSingle > 0)
    meanEdepSingle /= static_cast<double>(nSingle);

  // no single-hit layer: use the median edep of all hits
  if (nSingle == 0 && !hitEdepKeV.empty()) {
    std::vector<double> allEdeps;
    allEdeps.reserve(hitEdepKeV.size());
    for (const auto& [idx, e] : hitEdepKeV)
      allEdeps.push_back(e);
    std::sort(allEdeps.begin(), allEdeps.end());
    const size_t mid = allEdeps.size() / 2;
    meanEdepSingle = (allEdeps.size() % 2 == 0) ? 0.5 * (allEdeps[mid - 1] + allEdeps[mid]) : allEdeps[mid];
  }

  hitEdepScore.clear();
  hitEdepScore.reserve(allHitInfo.size());
  if (finder.m_edepScoreWeight.value() > 0.0) {
    for (size_t idx = 0; idx < allHitInfo.size(); ++idx) {
      const double edep_ = hitEdepKeV.at(idx);
      const double lowDev_ = std::max(0.0, meanEdepSingle - edep_);
      const double highDev_ = std::max(0.0, edep_ - 2.0 * meanEdepSingle);
      hitEdepScore[idx] = lowDev_ + highDev_;
    }

    const int maxHPC = finder.m_maxHitsPerComposite.value();
    if (maxHPC > 0) {
      for (auto& [cid, cidIndices] : hitIndicesByCompositeLayer) {
        if (static_cast<int>(cidIndices.size()) <= maxHPC)
          continue;
        std::sort(cidIndices.begin(), cidIndices.end(),
                  [this](size_t a, size_t b) { return hitEdepScore[a] < hitEdepScore[b]; });
        cidIndices.resize(maxHPC);
      }
    }
  }

  usedHits = globalUsedHits;

  // number of hits searched (after the trimming), and the size range of a combination
  int nAvail = 0;
  for (const auto& cid : compositeLayerIDs)
    nAvail += static_cast<int>(hitIndicesByCompositeLayer[cid].size());
  maxK = std::min(nAvail, finder.m_maxCombinatorialHits.value());
  minK = std::max(3, finder.m_commonConfig.minTrackHits);

  // hit index -> layer key
  gIdxToComposite.clear();
  gIdxToComposite.reserve(allHitInfo.size());
  for (const auto& hi : allHitInfo)
    gIdxToComposite[hi.index] = hi.compositeID;
  return true;
}

// SeedAndExtend: the next candidate (best first) whose hits are all still free.
bool CombinatorialTrackFinder::EventSearch::nextSeedAndExtendCandidate(k4rt::CombResult& best) {
  if (!seBuilt) {
    std::vector<std::pair<size_t, int>> pool;
    pool.reserve(allHitInfo.size());
    for (const auto& hi : allHitInfo)
      pool.emplace_back(hi.index, hi.compositeID);
    k4rt::SeedStats ss; // statistics of the search (filled by run, not used here)
    seCands = finder.m_seedFinder->run(hits, pool, ss);
    seBuilt = true;
  }
  std::unordered_map<size_t, size_t> gToHii;
  for (size_t i = 0; i < allHitInfo.size(); ++i)
    gToHii[allHitInfo[i].index] = i;
  while (seNext < seCands.size()) {
    const auto& c = seCands[seNext++];
    bool free = true;
    for (size_t g : c.gIdx)
      if (usedHits[g] || !gToHii.count(g)) {
        free = false;
        break;
      }
    if (!free)
      continue;
    best = k4rt::CombResult{};
    best.gIdxVec = c.gIdx;
    for (size_t g : c.gIdx)
      best.hiiVec.push_back(gToHii[g]);
    best.chi2ndf = c.chi2ndf;
    best.avgProxCm = 0.0;
    best.x0 = c.x0;
    best.y0 = c.y0;
    best.radius = c.R;
    best.fitCov = c.cov;
    best.residuals = c.residuals;
    best.nHits = static_cast<int>(c.gIdx.size());
    best.kinkChi2 = c.kinkChi2;
    best.kinkNdf = c.kinkNdf;
    best.looper = c.looper;
    return true;
  }
  return false;
}

// Isolation: each hit of the combination must be HitIsolationCut (cm) away from
// every other unused hit in a different layer. Hits in the same layer (secondaries
// next to the track) are not counted.
bool CombinatorialTrackFinder::EventSearch::isolated(const std::vector<size_t>& gIdxVec) const {
  const double cut = finder.m_hitIsolationCut.value();
  if (cut <= 0.0)
    return true;
  for (size_t gi : gIdxVec) {
    const auto refPos = hits[gi].getPosition();
    Eigen::Vector3d ref(refPos[0] / 10.0, refPos[1] / 10.0, refPos[2] / 10.0);
    int giComposite = -999;
    auto cgit = gIdxToComposite.find(gi);
    if (cgit != gIdxToComposite.end())
      giComposite = cgit->second;
    for (size_t j = 0; j < hits.size(); ++j) {
      if (j == gi)
        continue;
      if (globalUsedHits[j])
        continue;
      bool inCombo = false;
      for (size_t si : gIdxVec) {
        if (si == j) {
          inCombo = true;
          break;
        }
      }
      if (inCombo)
        continue;
      auto cjit = gIdxToComposite.find(j);
      if (cjit == gIdxToComposite.end())
        continue; // not in this iteration's hit pool
      if (cjit->second == giComposite)
        continue;
      const auto p = hits[j].getPosition();
      Eigen::Vector3d cand(p[0] / 10.0, p[1] / 10.0, p[2] / 10.0);
      if ((cand - ref).norm() < cut)
        return false;
    }
  }
  return true;
}

// --- Pair coincidence (RequirePairCoincidence), for detectors built of layer ---
// doublets (0,1), (2,3), (4,5): two hits of a doublet are coincident when
// |dphi| < 0.15 rad, |deta| < 0.10 and their distance is below PairCoincidenceMaxDist.
bool CombinatorialTrackFinder::EventSearch::coincident(const edm4hep::TrackerHitPlane& a,
                                                       const edm4hep::TrackerHitPlane& b) const {
  constexpr double kPairTolPhi = 0.15;
  constexpr double kPairTolEta = 0.10;
  const double kPairMaxDistMM = finder.m_pairCoincMaxDist.value() * 10.0; // cm -> mm (hit positions are in mm)
  const auto& pa = a.getPosition();
  const auto& pb = b.getPosition();
  double dx = pb[0] - pa[0];
  double dy = pb[1] - pa[1];
  double dz = pb[2] - pa[2];
  if ((dx * dx + dy * dy + dz * dz) > (kPairMaxDistMM * kPairMaxDistMM))
    return false;
  double dPhi = std::atan2(pb[1], pb[0]) - std::atan2(pa[1], pa[0]);
  dPhi = k4rt::wrapAngle(dPhi);
  if (std::abs(dPhi) > kPairTolPhi)
    return false;
  double rA = std::sqrt(pa[0] * pa[0] + pa[1] * pa[1]);
  double rB = std::sqrt(pb[0] * pb[0] + pb[1] * pb[1]);
  if (rA < 1e-6 || rB < 1e-6)
    return false;
  double etaA = std::asinh(pa[2] / rA);
  double etaB = std::asinh(pb[2] / rB);
  return std::abs(etaB - etaA) <= kPairTolEta;
}

// Complete doublets of a combination: hits in both layers of one of (0,1)/(2,3)/(4,5), in the
// same region, and mutually coincident.
int CombinatorialTrackFinder::EventSearch::completeDoublets(const std::vector<size_t>& gIdxVec) const {
  int n_complete = 0;
  const int pairs[3][2] = {{0, 1}, {2, 3}, {4, 5}};
  for (int p = 0; p < 3; ++p) {
    bool complete = false;
    for (size_t gi_a : gIdxVec) {
      auto itA = gIdxToComposite.find(gi_a);
      if (itA == gIdxToComposite.end())
        continue;
      int cidA = itA->second;
      if (k4rt::layerOfKey(cidA) != pairs[p][0])
        continue;
      int regionA = k4rt::regionOfKey(cidA);
      for (size_t gi_b : gIdxVec) {
        if (gi_b == gi_a)
          continue;
        auto itB = gIdxToComposite.find(gi_b);
        if (itB == gIdxToComposite.end())
          continue;
        int cidB = itB->second;
        if (k4rt::layerOfKey(cidB) != pairs[p][1])
          continue;
        if (k4rt::regionOfKey(cidB) != regionA)
          continue;
        if (coincident(hits[gi_a], hits[gi_b])) {
          complete = true;
          break;
        }
      }
      if (complete)
        break;
    }
    if (complete)
      ++n_complete;
  }
  return n_complete;
}

// --- Hit clustering (DoHitClustering) ---
// Hits of neighbouring layers that pass the road cuts are joined (union-find); each cluster is
// searched on its own. Without clustering: one cluster, all hits.
std::vector<ClusterDef>
CombinatorialTrackFinder::EventSearch::buildClusters(const std::vector<std::vector<size_t>>& hitsForCompAll) const {
  std::vector<ClusterDef> activeClusters;
  const int nDistinctAll = static_cast<int>(compositeLayerIDs.size());
  if (!(finder.m_doHitClustering.value() && nDistinctAll >= 2)) {
    ClusterDef cd;
    cd.cids = compositeLayerIDs;
    cd.hitsForComp = hitsForCompAll;
    activeClusters.push_back(std::move(cd));
    return activeClusters;
  }
  UnionFind groups(allHitInfo.size());
  const double phiCut = finder.m_maxConsecDeltaPhi.value();
  const double distCut = finder.m_maxConsecHitDist.value(); // cm

  joinLayers(groups, hitsForCompAll, 1, phiCut, distCut);
  // also across one missing layer, with doubled cuts
  joinLayers(groups, hitsForCompAll, 2, (phiCut > 0.0) ? phiCut * 2.0 : -1.0, (distCut > 0.0) ? distCut * 2.0 : -1.0);

  std::unordered_map<size_t, std::vector<size_t>> rootToHIIs;
  for (int ci = 0; ci < nDistinctAll; ++ci) {
    for (size_t hii : hitsForCompAll[ci])
      rootToHIIs[groups.find(hii)].push_back(hii);
  }

  // keep clusters with at least minK layers, layers in inner -> outer order
  for (auto& [root, hiis] : rootToHIIs) {
    std::unordered_map<int, std::vector<size_t>> cidMap;
    for (size_t hii : hiis)
      cidMap[allHitInfo[hii].compositeID].push_back(hii);
    if (static_cast<int>(cidMap.size()) < minK)
      continue;
    ClusterDef cd;
    for (int cid : compositeLayerIDs) {
      auto it = cidMap.find(cid);
      if (it == cidMap.end())
        continue;
      cd.cids.push_back(cid);
      cd.hitsForComp.push_back(std::move(it->second));
    }
    activeClusters.push_back(std::move(cd));
  }
  return activeClusters;
}

// Joins the hits of layers `step` apart (1: neighbouring layers; 2: across one missing layer)
// that pass the phi and distance limits into one cluster.
void CombinatorialTrackFinder::EventSearch::joinLayers(UnionFind& groups,
                                                       const std::vector<std::vector<size_t>>& hitsForCompAll, int step,
                                                       double phiLimit, double distLimit) const {
  const int nDistinctAll = static_cast<int>(hitsForCompAll.size());
  for (int ci = 0; ci + step < nDistinctAll; ++ci) {
    for (size_t hii_a : hitsForCompAll[ci]) {
      const auto pa = hits[allHitInfo[hii_a].index].getPosition();
      const double phiA = std::atan2(pa[1], pa[0]);
      for (size_t hii_b : hitsForCompAll[ci + step]) {
        const auto pb = hits[allHitInfo[hii_b].index].getPosition();
        if (phiLimit > 0.0) {
          double dPhi = std::atan2(pb[1], pb[0]) - phiA;
          dPhi = k4rt::wrapAngle(dPhi);
          if (std::abs(dPhi) > phiLimit)
            continue;
        }
        if (distLimit > 0.0) {
          double dx = pb[0] - pa[0], dy = pb[1] - pa[1], dz = pb[2] - pa[2];
          if (std::sqrt(dx * dx + dy * dy + dz * dz) * 0.1 > distLimit)
            continue;
        }
        groups.unite(hii_a, hii_b);
      }
    }
  }
}

// --- N-hit combinatorial search (Strategy = Exhaustive) ---
// All C(N,k) combinations for k = MaxCombinatorialHits down to MinTrackHits, each scored by the
// chi2/ndf of its circle fit, after road, isolation and quality cuts.
bool CombinatorialTrackFinder::EventSearch::searchCombinations(k4rt::CombResult& best) {
  // Road cuts (k4rt::RoadCuts, ExhaustiveSearch.h): quick geometric tests that reject a
  // combination of hits before its circle fit.
  k4rt::RoadCutConfig roadCfg;
  roadCfg.maxConsecDeltaPhi = finder.m_maxConsecDeltaPhi.value();
  roadCfg.maxComboPhiSpread = finder.m_maxComboPhiSpread.value();
  roadCfg.maxConsecHitDist = finder.m_maxConsecHitDist.value();
  roadCfg.maxPairHitDist = finder.m_commonConfig.maxPairHitDistanceCm;
  roadCfg.maxComboTanLambdaDev = finder.m_maxComboTanLambdaDev.value();
  const k4rt::RoadCuts roadCuts(roadCfg, hits, gIdxToComposite, &finder.m_cidRegion);

  // unused hits of each layer: hitsForCompAll[i] belongs to compositeLayerIDs[i]
  const int nDistinctAll = static_cast<int>(compositeLayerIDs.size());
  std::vector<std::vector<size_t>> hitsForCompAll(nDistinctAll);
  for (int ci = 0; ci < nDistinctAll; ++ci) {
    for (size_t hii : hitIndicesByCompositeLayer.at(compositeLayerIDs[ci])) {
      if (!usedHits[allHitInfo[hii].index])
        hitsForCompAll[ci].push_back(hii);
    }
  }

  // Pair coincidence: first pass only with hits that have a coincident partner
  std::vector<ClusterDef> activeClustersSeed;
  if (pcHard) {
    ClusterDef cd;
    cd.cids = compositeLayerIDs;
    cd.hitsForComp.assign(nDistinctAll, std::vector<size_t>{});
    for (int ci = 0; ci < nDistinctAll; ++ci) {
      for (size_t hii : hitsForCompAll[ci]) {
        // the partner layer of the hit's doublet: 0 <-> 1, 2 <-> 3, 4 <-> 5 (same region and subdetector)
        const int composite = allHitInfo[hii].compositeID;
        const int layer = k4rt::layerOfKey(composite);
        const int absPartner = std::abs(composite) - layer + (layer ^ 1);
        const int partnerKey = (composite < 0) ? -absPartner : absPartner;
        // the hit is kept if a hit of the partner layer is coincident with it
        auto partnerHits = hitIndicesByCompositeLayer.find(partnerKey);
        if (partnerHits == hitIndicesByCompositeLayer.end())
          continue;
        for (size_t partnerHii : partnerHits->second) {
          if (coincident(allHitInfo[hii].hit, allHitInfo[partnerHii].hit)) {
            cd.hitsForComp[ci].push_back(hii);
            break;
          }
        }
      }
    }
    activeClustersSeed.push_back(std::move(cd));
  }
  const std::vector<ClusterDef> activeClusters = buildClusters(hitsForCompAll);

  // passes: hard pair coincidence = coincident hits first, then all hits if the best
  // combination has < 3 complete doublets; otherwise one pass over all hits
  k4rt::CombResult fallback; // best over-threshold combo; promoted if no valid combo exists
  bool found = false;
  const int nPairPasses = pcHard ? 2 : 1;
  for (int pairPassIdx = 0; pairPassIdx < nPairPasses; ++pairPassIdx) {
    // second pass only if the first one's best has fewer than 3 complete doublets
    if (pairPassIdx == 1 && found && completeDoublets(best.gIdxVec) >= 3)
      break;
    const std::vector<ClusterDef>& clustersThisPass =
        (pcHard && pairPassIdx == 0) ? activeClustersSeed : activeClusters;
    for (const auto& cluster : clustersThisPass)
      searchCluster(cluster, roadCuts, best, fallback, found);
  }

  if (!found && fallback.nHits > 0) {
    // every combination was above MaxComboPT: use the fallback
    best = fallback;
    found = true;
  }
  return found;
}

// All combinations of one cluster: for k = most hits down to minK, every choice of k layers and,
// for each, every combination of one hit per chosen layer. The best combination so far is kept
// in `best`, the best one above MaxComboPT in `fallback`.
void CombinatorialTrackFinder::EventSearch::searchCluster(const ClusterDef& cluster, const k4rt::RoadCuts& roadCuts,
                                                          k4rt::CombResult& best, k4rt::CombResult& fallback,
                                                          bool& found) {
  const std::vector<std::vector<size_t>>& hitsForComp = cluster.hitsForComp;
  const int nDistinct = static_cast<int>(cluster.cids.size());
  const int maxKeff = std::min(maxK, nDistinct);

  for (int k = maxKeff; k >= minK; --k) {
    // the chosen layers: indices into the cluster's layers, the first choice is 0, 1, ..., k-1
    std::vector<int> cidxCombo(k);
    for (int i = 0; i < k; ++i)
      cidxCombo[i] = i;

    while (true) {
      // the hits of each chosen layer
      bool anySlotEmpty = false;
      std::vector<std::vector<size_t>> slotHits(k);
      for (int i = 0; i < k; ++i) {
        slotHits[i] = hitsForComp[cidxCombo[i]];
        if (slotHits[i].empty()) {
          anySlotEmpty = true;
          break;
        }
      }

      if (!anySlotEmpty) {
        // every combination of one hit per layer (pos[i]: the hit taken in layer i)
        std::vector<int> pos(k, 0);
        while (true) {
          std::vector<size_t> hiiVec(k), gIdxVec(k);
          for (int i = 0; i < k; ++i) {
            hiiVec[i] = slotHits[i][pos[i]];
            gIdxVec[i] = allHitInfo[hiiVec[i]].index;
          }

          bool anyUsed = false;
          for (size_t gi : gIdxVec) {
            if (usedHits[gi]) {
              anyUsed = true;
              break;
            }
          }

          if (!anyUsed) {
            if (!roadCuts.passesRoadCuts(gIdxVec)) {
              ++finder.m_rejectedRoadCuts;
            } else if (!roadCuts.passesDistanceCuts(gIdxVec)) {
              ++finder.m_rejectedDistanceCuts;
            } else if (!isolated(gIdxVec)) {
              // rejected: a hit is too close to another unused hit (HitIsolationCut)
            } else if (!roadCuts.passesTanLambdaCut(gIdxVec)) {
              ++finder.m_rejectedTanLambdaCut;
            } else if (pcHard && completeDoublets(gIdxVec) < 2) {
              // rejected: hard pair coincidence needs at least 2 complete doublets
            } else {
              ++finder.m_combinationsTried;
              std::vector<edm4hep::TrackerHitPlane> comboHits;
              comboHits.reserve(k);
              for (size_t gi : gIdxVec)
                comboHits.push_back(hits[gi]);

              // circle parameters start at 0, so the pT below is defined even if the fit fails
              double cx0 = 0.0, cy0 = 0.0, cr = 0.0, cchi2ndf = 0.0;
              Eigen::Matrix3d cCov = Eigen::Matrix3d::Zero();
              std::vector<double> cRes;
              bool fitOK = false;
              try {
                fitOK = finder.fitCircleNHits(comboHits, cx0, cy0, cr, cchi2ndf, cCov, cRes);
              } catch (...) {
              }

              // pT of the combination, from its circle radius
              const double comboPT = k4rt::ptFromRadius(absBRefForK3, cr);
              if (!fitOK) {
                // the circle fit failed
              } else if (cchi2ndf >= finder.m_nHitMaxChi2NDF.value()) {
                ++finder.m_rejectedChi2Cut;
              } else if (comboPT > finder.m_maxComboPT.value()) {
                // pT above MaxComboPT: not a candidate, but kept as fallback in
                // case every combination is above it (near-straight tracks)
                ++finder.m_rejectedPtCuts;
                // best fallback: more hits; then lower chi2/ndf (k >= 4) or
                // lower pT (k = 3)
                bool isBetterFallback = false;
                if (k != fallback.nHits) {
                  isBetterFallback = (k > fallback.nHits);
                } else if (k == 3) {
                  isBetterFallback = (cr < fallback.radius);
                } else {
                  isBetterFallback = (cchi2ndf < fallback.chi2ndf);
                }
                if (isBetterFallback) {
                  fallback.hiiVec = hiiVec;
                  fallback.gIdxVec = gIdxVec;
                  fallback.chi2ndf = cchi2ndf;
                  fallback.avgProxCm = 0.0;
                  fallback.x0 = cx0;
                  fallback.y0 = cy0;
                  fallback.radius = cr;
                  fallback.fitCov = cCov;
                  fallback.residuals = cRes;
                  fallback.nHits = k;
                }
              } else if (finder.m_commonConfig.minPtGeV > 0.0 && comboPT < finder.m_commonConfig.minPtGeV) {
                // pT below MinComboPT: rejected (no fallback)
                ++finder.m_rejectedPtCuts;
              } else {
                // Ranking of the combinations (lower score is better):
                //  1. more hits always wins;
                //  2. k >= 4: chi2/ndf + w * edep score / EdepNormKeV;
                //  3. k = 3 (chi2 is always 0 for 3 points):
                //     w * edep score / EdepNormKeV - min(R, R_cap) / R_cap,
                //     i.e. straighter (higher pT) and MIP-like is better.
                // The edep score is the mean hit score against the event
                // baseline (see collectHits).
                double cEdepDev = 0.0;
                const double w = finder.m_edepScoreWeight.value();
                if (w > 0.0 && !hitEdepScore.empty()) {
                  for (size_t hii_i : hiiVec)
                    cEdepDev += hitEdepScore.at(hii_i);
                  cEdepDev /= static_cast<double>(k);
                }
                const double edepNorm = std::max(finder.m_edepNormKeV.value(), 1e-9);
                bool isBetter;
                if (k != best.nHits) {
                  isBetter = (k > best.nHits);
                } else if (k == 3) {
                  // k = 3 (radius term saturates at R_cap; w = 0: radius only)
                  if (w > 0.0) {
                    const double cRadNorm = std::min(cr, rCapK3) / rCapK3;
                    const double bRadNorm = std::min(best.radius, rCapK3) / rCapK3;
                    const double cScore3 = w * cEdepDev / edepNorm - cRadNorm;
                    const double bScore3 = w * best.avgProxCm / edepNorm - bRadNorm;
                    isBetter = (cScore3 < bScore3);
                  } else {
                    isBetter = (cr > best.radius);
                  }
                } else {
                  const double cScore = cchi2ndf + (w > 0.0 ? w * cEdepDev / edepNorm : 0.0);
                  const double bScore = best.chi2ndf + (w > 0.0 ? w * best.avgProxCm / edepNorm : 0.0);
                  isBetter = (cScore < bScore);
                }
                if (isBetter) {
                  best.hiiVec = hiiVec;
                  best.gIdxVec = gIdxVec;
                  best.chi2ndf = cchi2ndf;
                  best.avgProxCm = cEdepDev;
                  best.x0 = cx0;
                  best.y0 = cy0;
                  best.radius = cr;
                  best.fitCov = cCov;
                  best.residuals = cRes;
                  best.nHits = k;
                  found = true;
                }
              }
            }
          }

          // next combination of hits (like an odometer: the last layer turns first)
          int carry = k - 1;
          while (carry >= 0) {
            pos[carry]++;
            if (pos[carry] < static_cast<int>(slotHits[carry].size()))
              break;
            pos[carry] = 0;
            carry--;
          }
          if (carry < 0)
            break;
        }
      }

      // next choice of k layers (lexicographic order of the index sets)
      int i = k - 1;
      while (i >= 0 && cidxCombo[i] == nDistinct - k + i)
        --i;
      if (i < 0)
        break;
      ++cidxCombo[i];
      for (int j = i + 1; j < k; ++j)
        cidxCombo[j] = cidxCombo[j - 1] + 1;
    }
  }
}

// clean hits (layers with exactly one unused hit), ordered by layer key
void CombinatorialTrackFinder::EventSearch::collectCleanHits(std::vector<edm4hep::TrackerHitPlane>& cleanHits,
                                                             std::vector<size_t>& cleanHIIs,
                                                             std::vector<size_t>& cleanGIdxs) const {
  for (const auto& [cid, cidHIIs] : hitIndicesByCompositeLayer) {
    if (static_cast<int>(cidHIIs.size()) != 1)
      continue;
    size_t hii = cidHIIs[0];
    if (globalUsedHits[allHitInfo[hii].index])
      continue;
    cleanHIIs.push_back(hii);
    cleanGIdxs.push_back(allHitInfo[hii].index);
    cleanHits.push_back(allHitInfo[hii].hit);
  }
}

// refit the guided choice; keep it if it is different and better
void CombinatorialTrackFinder::EventSearch::applyGuidedChoice(const std::vector<edm4hep::TrackerHitPlane>& guidedHits,
                                                              const std::vector<size_t>& guidedHIIs,
                                                              const std::vector<size_t>& guidedGIdxs,
                                                              k4rt::CombResult& best) const {
  double gx0 = 0, gy0 = 0, gr = 0, gchi2 = 0;
  Eigen::Matrix3d gCov;
  std::vector<double> gRes;
  bool guidedOK = false;
  try {
    guidedOK = finder.fitCircleNHits(guidedHits, gx0, gy0, gr, gchi2, gCov, gRes);
  } catch (...) {
  }
  if (!guidedOK || gchi2 > finder.m_nHitMaxChi2NDF.value())
    return;
  bool differentChoice = false;
  for (size_t gHII : guidedHIIs) {
    bool inGeneral = false;
    for (size_t bHII : best.hiiVec)
      if (bHII == gHII) {
        inGeneral = true;
        break;
      }
    if (!inGeneral) {
      differentChoice = true;
      break;
    }
  }
  if (differentChoice && gchi2 < best.chi2ndf) {
    best.hiiVec = guidedHIIs;
    best.gIdxVec = guidedGIdxs;
    best.chi2ndf = gchi2;
    best.x0 = gx0;
    best.y0 = gy0;
    best.radius = gr;
    best.fitCov = gCov;
    best.residuals = gRes;
    best.nHits = static_cast<int>(guidedHIIs.size());
  }
}

// --- Guided hit choice in crowded layers (doCrowdedLayerHitSelection) ---
// Layers with one hit are "clean". In each crowded layer the hit is chosen by
//   A (>= 3 clean layers): distance to the circle through the clean hits,
//      score 0.75 * dr / SigmaHitDefault + 0.25 * edep score / EdepNormKeV;
//   B (2 clean layers): kink angle to the line through the two clean hits,
//      score 0.75 * kink / MaxConsecDeltaPhi + 0.25 * edep score / EdepNormKeV.
// The result replaces the best combination if its refit has a lower chi2/ndf.
void CombinatorialTrackFinder::EventSearch::selectCrowdedLayerHits(k4rt::CombResult& best) {
  int nCleanCIDs = 0;
  int nCrowdedCIDs = 0;
  for (const auto& [cid, cidHIIs] : hitIndicesByCompositeLayer) {
    if (static_cast<int>(cidHIIs.size()) == 1)
      ++nCleanCIDs;
    else if (static_cast<int>(cidHIIs.size()) > 1)
      ++nCrowdedCIDs;
  }
  const int totalCIDs = nCleanCIDs + nCrowdedCIDs;
  const double edepNorm = std::max(finder.m_edepNormKeV.value(), 1e-9); // keV

  // --- A: circle through the clean hits ---
  if (nCleanCIDs >= 3 && nCrowdedCIDs >= 1) {
    std::vector<edm4hep::TrackerHitPlane> cleanHits;
    std::vector<size_t> cleanHIIs, cleanGIdxs;
    collectCleanHits(cleanHits, cleanHIIs, cleanGIdxs);
    if (static_cast<int>(cleanHits.size()) < 3)
      return;
    double cx0 = 0, cy0 = 0, cr = 0, cchi2 = 0;
    Eigen::Matrix3d cCov;
    std::vector<double> cRes;
    bool cleanOK = false;
    try {
      cleanOK = finder.fitCircleNHits(cleanHits, cx0, cy0, cr, cchi2, cCov, cRes);
    } catch (...) {
    }
    if (!cleanOK || !(cr > 0))
      return;
    const double drNorm = std::max(finder.m_sigmaHitDefault.value(), 1e-9); // cm
    std::vector<edm4hep::TrackerHitPlane> guidedHits = cleanHits;
    std::vector<size_t> guidedHIIs = cleanHIIs;
    std::vector<size_t> guidedGIdxs = cleanGIdxs;
    for (const auto& [cid, cidHIIs] : hitIndicesByCompositeLayer) {
      if (static_cast<int>(cidHIIs.size()) <= 1)
        continue;
      size_t bestHII = SIZE_MAX;
      double bestScore = std::numeric_limits<double>::max();
      for (size_t hii : cidHIIs) {
        if (globalUsedHits[allHitInfo[hii].index])
          continue;
        const auto& hp = allHitInfo[hii].hit.getPosition(); // mm
        double px = hp.x * 0.1, py = hp.y * 0.1;            // -> cm
        double dr = std::abs(std::sqrt((px - cx0) * (px - cx0) + (py - cy0) * (py - cy0)) - cr);
        double edepDev = hitEdepScore.count(hii) ? hitEdepScore.at(hii) : 0.0;
        double score = 0.75 * (dr / drNorm) + 0.25 * (edepDev / edepNorm);
        if (score < bestScore) {
          bestScore = score;
          bestHII = hii;
        }
      }
      if (bestHII == SIZE_MAX)
        return; // a crowded layer without a free hit
      guidedHIIs.push_back(bestHII);
      guidedGIdxs.push_back(allHitInfo[bestHII].index);
      guidedHits.push_back(allHitInfo[bestHII].hit);
    }
    applyGuidedChoice(guidedHits, guidedHIIs, guidedGIdxs, best);
  }
  // --- B: line through the two clean hits ---
  else if (nCleanCIDs == 2 && nCrowdedCIDs >= 1 && totalCIDs >= 3) {
    std::vector<edm4hep::TrackerHitPlane> cleanHits;
    std::vector<size_t> cleanHIIs, cleanGIdxs;
    collectCleanHits(cleanHits, cleanHIIs, cleanGIdxs);
    if (static_cast<int>(cleanHits.size()) != 2)
      return;
    // direction from the inner (A) to the outer (B) clean hit, cm
    const auto& posA_mm = cleanHits[0].getPosition();
    const auto& posB_mm = cleanHits[1].getPosition();
    double ax = posA_mm.x * 0.1, ay = posA_mm.y * 0.1, az = posA_mm.z * 0.1;
    double bx = posB_mm.x * 0.1, by = posB_mm.y * 0.1, bz = posB_mm.z * 0.1;
    double dabx = bx - ax, daby = by - ay, dabz = bz - az;
    double lenAB = std::sqrt(dabx * dabx + daby * daby + dabz * dabz);
    if (!(lenAB > 1e-9))
      return;                                                                     // degenerate if A == B
    const double sigmaKappa = std::max(finder.m_maxConsecDeltaPhi.value(), 1e-4); // rad
    std::vector<edm4hep::TrackerHitPlane> guidedHits = cleanHits;
    std::vector<size_t> guidedHIIs = cleanHIIs;
    std::vector<size_t> guidedGIdxs = cleanGIdxs;
    for (const auto& [cid, cidHIIs] : hitIndicesByCompositeLayer) {
      if (static_cast<int>(cidHIIs.size()) <= 1)
        continue;
      size_t bestHII = SIZE_MAX;
      double bestScore = std::numeric_limits<double>::max();
      for (size_t hii : cidHIIs) {
        if (globalUsedHits[allHitInfo[hii].index])
          continue;
        const auto& hp = allHitInfo[hii].hit.getPosition(); // mm
        double cx = hp.x * 0.1, cy = hp.y * 0.1, cz = hp.z * 0.1;
        double dacx = cx - ax, dacy = cy - ay, dacz = cz - az;
        double lenAC = std::sqrt(dacx * dacx + dacy * dacy + dacz * dacz);
        // angle at A between AB and AC
        double kappa = 0.0;
        if (lenAC > 1e-9) {
          double cosK = (dabx * dacx + daby * dacy + dabz * dacz) / (lenAB * lenAC);
          cosK = std::max(-1.0, std::min(1.0, cosK));
          kappa = std::acos(cosK); // rad
        }
        double edepDev = hitEdepScore.count(hii) ? hitEdepScore.at(hii) : 0.0;
        double score = 0.75 * (kappa / sigmaKappa) + 0.25 * (edepDev / edepNorm);
        if (score < bestScore) {
          bestScore = score;
          bestHII = hii;
        }
      }
      if (bestHII == SIZE_MAX)
        return; // a crowded layer without a free hit
      guidedHIIs.push_back(bestHII);
      guidedGIdxs.push_back(allHitInfo[bestHII].index);
      guidedHits.push_back(allHitInfo[bestHII].hit);
    }
    applyGuidedChoice(guidedHits, guidedHIIs, guidedGIdxs, best);
  }
}

// --- Neighbour-track cuts ---
// A combination with a hit within NeighbourTrackMaxDist (cm) of an accepted track may be a
// delta-ray shadow. It is kept only with >= NeighbourTrackMinHits hits and chi2/ndf <=
// NeighbourTrackMaxChi2NDF; otherwise its hits near the accepted track are marked as used and
// true is returned (the search goes on with the next iteration).
bool CombinatorialTrackFinder::EventSearch::rejectedAsNeighbour(const k4rt::CombResult& best) {
  if (!(finder.m_neighbourTrackMaxDist.value() > 0.0) || acceptedTrackHitPositions.empty())
    return false;
  const double distThreshMm = finder.m_neighbourTrackMaxDist.value() * 10.0; // cm -> mm
  const double distThreshSq = distThreshMm * distThreshMm;
  bool isNeighbour = false;
  for (size_t hii : best.hiiVec) {
    if (nearAcceptedTrack(hii, distThreshSq)) {
      isNeighbour = true;
      break;
    }
  }
  if (!isNeighbour)
    return false;
  const bool passHits = (static_cast<int>(best.hiiVec.size()) >= finder.m_neighbourTrackMinHits.value());
  const bool passChi2 = (best.chi2ndf <= finder.m_neighbourTrackMaxChi2NDF.value());
  if (passHits && passChi2)
    return false;
  // only the hits near the accepted track are used up
  for (size_t hii : best.hiiVec) {
    if (nearAcceptedTrack(hii, distThreshSq))
      globalUsedHits[allHitInfo[hii].index] = true;
  }
  return true;
}

// is the hit allHitInfo[hii] closer than sqrt(distanceSquaredMm2) to a hit of an accepted track?
bool CombinatorialTrackFinder::EventSearch::nearAcceptedTrack(size_t hii, double distanceSquaredMm2) const {
  const auto& hp = allHitInfo[hii].hit.getPosition(); // mm
  for (const auto& ap : acceptedTrackHitPositions) {
    double dx = hp.x - ap[0], dy = hp.y - ap[1], dz = hp.z - ap[2];
    if (dx * dx + dy * dy + dz * dz < distanceSquaredMm2)
      return true;
  }
  return false;
}

// --- Outlier removal: drop the worst hit while its residual is above ---
//    OutlierSigma * SigmaHitDefault, and refit
void CombinatorialTrackFinder::EventSearch::removeOutliers(k4rt::CombResult& best) const {
  const double outlierCutCm = finder.m_outlierSigma.value() * finder.m_sigmaHitDefault.value();
  for (int iter = 0; iter < finder.m_maxOutlierIterations.value(); ++iter) {
    if (static_cast<int>(best.gIdxVec.size()) <= minK)
      break;

    size_t worstIdx = 0;
    double worstRes = 0.0;
    for (size_t ri = 0; ri < best.residuals.size(); ++ri) {
      if (best.residuals[ri] > worstRes) {
        worstRes = best.residuals[ri];
        worstIdx = ri;
      }
    }
    if (worstRes <= outlierCutCm)
      break;

    best.gIdxVec.erase(best.gIdxVec.begin() + worstIdx);
    best.hiiVec.erase(best.hiiVec.begin() + worstIdx);
    best.residuals.erase(best.residuals.begin() + worstIdx);
    best.nHits = static_cast<int>(best.gIdxVec.size());

    if (best.nHits < minK)
      break;

    std::vector<edm4hep::TrackerHitPlane> refitHits;
    for (size_t gi : best.gIdxVec)
      refitHits.push_back(hits[gi]);
    double rx0, ry0, rr, rchi2;
    Eigen::Matrix3d rCov;
    std::vector<double> rRes;
    bool refitOK = false;
    try {
      refitOK = finder.fitCircleNHits(refitHits, rx0, ry0, rr, rchi2, rCov, rRes);
    } catch (...) {
    }
    if (refitOK) {
      best.x0 = rx0;
      best.y0 = ry0;
      best.radius = rr;
      best.chi2ndf = rchi2;
      best.fitCov = rCov;
      best.residuals = rRes;
    }
  }
}

// The chosen combination as an EDM4hep track: its hits are marked as used, the fitted helix is
// stored as track states, the hits and their reco-sim links are copied to the output.
void CombinatorialTrackFinder::EventSearch::buildTrack(const k4rt::CombResult& best) {
  // pT of the combination: from its circle radius and the field at the centroid of its hits
  // (inside or outside the solenoid); used by the SolenoidCrossing AtIP state
  double pT = 0.0;
  {
    const int Nh = best.nHits;
    double sumXf = 0, sumYf = 0, sumZf = 0;
    for (size_t gi : best.gIdxVec) {
      const auto p = hits[gi].getPosition();
      sumXf += p[0] / 10.0;
      sumYf += p[1] / 10.0;
      sumZf += p[2] / 10.0;
    }
    dd4hep::Position avgPosN(sumXf / Nh, sumYf / Nh, sumZf / Nh);
    double bFieldN = finder.m_field.magneticField(avgPosN).z() / dd4hep::tesla;
    pT = k4rt::ptFromRadius(bFieldN, best.radius); // GeV/c
  }

  // the track's hits, marked as used
  std::vector<edm4hep::TrackerHitPlane> trackHits;
  trackHits.reserve(best.gIdxVec.size());
  for (size_t idx : best.gIdxVec) {
    globalUsedHits[idx] = true;
    trackHits.push_back(hits[idx]);
  }

  // for the neighbour-track cuts of the next iterations
  if (finder.m_neighbourTrackMaxDist.value() > 0.0) {
    for (const auto& th : trackHits) {
      const auto& p = th.getPosition(); // mm
      acceptedTrackHitPositions.push_back({p.x, p.y, p.z});
    }
  }

  auto finalTrack = finalTracks.create();

  // The fitted helix as EDM4hep track states (makeEdm4hepStateFromCircle): at the first
  // hit, at the last hit and, with IPStateMode = Helix, at the closest approach to the origin.
  if (trackHits.size() >= 3) {
    if (best.radius > 0) {
      struct StoredState {
        StatePoint where;
        int location;
        bool wanted;
      };
      const StoredState stored[] = {
          {StatePoint::FirstHit, edm4hep::TrackState::AtFirstHit, true},
          {StatePoint::LastHit, edm4hep::TrackState::AtLastHit, true},
          {StatePoint::Origin, edm4hep::TrackState::AtIP, finder.m_ipStateMode.value() == "Helix"},
      };
      for (const auto& entry : stored) {
        if (!entry.wanted)
          continue;
        edm4hep::TrackState state;
        if (!makeEdm4hepStateFromCircle(best.x0, best.y0, best.radius, best.fitCov, trackHits, entry.where, state)) {
          finder.warning() << "Could not build track state " << entry.location << " of track " << trackNumber
                           << " (degenerate hits)" << endmsg;
          continue;
        }
        state.location = entry.location;
        finalTrack.addToTrackStates(state);
      }
      // the circle fit has one radial residual per hit and 3 parameters
      const int fitNdf = std::max(1, static_cast<int>(trackHits.size()) - 3);
      finalTrack.setChi2(static_cast<float>(best.chi2ndf * fitNdf)); // raw chi2
      finalTrack.setNdf(fitNdf);
    } else {
      finder.warning() << "Circle fit failed for track " << trackNumber << endmsg;
    }
  }

  // circle-fit ndf (n_hits - 3)
  finalTrack.setNdf(std::max(1, static_cast<int>(trackHits.size()) - 3));
  // SeedAndExtend: the stored chi2 is the kink chi2 (with multiple scattering);
  // loopers get track type 1
  if (best.kinkNdf > 0) {
    finalTrack.setChi2(static_cast<float>(best.kinkChi2));
    finalTrack.setNdf(best.kinkNdf);
  }
  if (best.looper)
    finalTrack.setType(1);

  // copy the hits to the output collection
  for (const auto& hit : trackHits) {
    auto outHit = outputHits.create();
    outHit.setCellID(hit.getCellID());
    outHit.setTime(hit.getTime());
    outHit.setEDep(hit.getEDep());
    outHit.setEDepError(hit.getEDepError());
    outHit.setPosition(hit.getPosition());
    outHit.setCovMatrix(hit.getCovMatrix());
    outHit.setDu(hit.getDu());
    outHit.setDv(hit.getDv());
    outHit.setU(hit.getU());
    outHit.setV(hit.getV());
    propagateLink(outHit, hit);
    finalTrack.addToTrackerHits(outHit);
  }

  // ===== IPStateMode = SolenoidCrossing: AtIP from the helix matched at the solenoid =====
  // Only for tracks with all hits outside the solenoid. The rotation sense and the dip
  // (tanLambda) of the outer helix come from its stored state at the first hit.
  edm4hep::TrackState firstHitState;
  bool haveFirstHitState = false;
  for (const auto& state : finalTrack.getTrackStates()) {
    if (state.location == edm4hep::TrackState::AtFirstHit) {
      firstHitState = state;
      haveFirstHitState = true;
    }
  }
  if (finder.m_ipStateMode.value() == "SolenoidCrossing" && finder.allHitsOutsideSolenoid(trackHits) &&
      haveFirstHitState && firstHitState.omega != 0.0f && finalTrack.trackerHits_size() > 0) {
    // the innermost hit anchors the crossing (cm)
    auto refHit = finalTrack.getTrackerHits(0);
    double xr = refHit.getPosition()[0] / 10.0; // mm -> cm
    double yr = refHit.getPosition()[1] / 10.0;
    double zr = refHit.getPosition()[2] / 10.0;
    // barrel hit: crossing of the solenoid barrel; endcap hit: of its end face
    bool isEndcapHit = (finder.getTypeID(refHit.getCellID()) != 0);
    // EDM4hep: omega > 0 means clockwise seen from +z
    const bool clockwise = (firstHitState.omega > 0.0f);
    bool atIPok = finder.analyticalInnerPropagation(best.x0, best.y0, best.radius, clockwise, pT,
                                                    firstHitState.tanLambda, xr, yr, zr, isEndcapHit, finalTrack);
    finder.m_solenoidCrossingBuilt += atIPok;
  }

  // ---- per-track run statistics ----
  finder.m_hitsPerTrack += trackHits.size();
}

void CombinatorialTrackFinder::findTracks(const edm4hep::TrackerHitPlaneCollection* hits,
                                          edm4hep::TrackCollection& finalTracks,
                                          edm4hep::TrackerHitPlaneCollection& outputHits,
                                          const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks,
                                          edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinks) const {
  EventSearch search(*this, *hits, finalTracks, outputHits, recoSimLinks, outputLinks);
  search.run();
}

// The two points where a circle (centre (cx, cy) at the distance d from the origin, radius R)
// crosses the solenoid cylinder (radius rSolenoid around the beam): an intersection of two
// circles. On the line from the origin towards the centre (unit vector e) the crossings lie at
//   h = (rSolenoid^2 + d^2 - R^2) / (2 d)   from the origin,
// and +-l = sqrt(rSolenoid^2 - h^2) away from that line along the perpendicular unit vector p.
static void circleCylinderCrossings(double cx, double cy, double R, double d, double rSolenoid, double& P1x,
                                    double& P1y, double& P2x, double& P2y) {
  double h = (rSolenoid * rSolenoid + d * d - R * R) / (2.0 * d);
  double l = std::sqrt(std::max(0.0, rSolenoid * rSolenoid - h * h));
  double ex = cx / d, ey = cy / d;
  double px = -ey, py = ex;
  P1x = h * ex + l * px;
  P1y = h * ey + l * py;
  P2x = h * ex - l * px;
  P2y = h * ey - l * py;
}

// Angle on a circle (centre (cx, cy)) from the reference hit (at the angle thetaHit about the
// centre) back to the point P, in [0, 2pi), going against the direction of travel (the track came
// from the IP side).
static double backwardArcAngle(double cx, double cy, double thetaHit, bool clockwise, double Px, double Py) {
  double tc = std::atan2(Py - cy, Px - cx);
  return clockwise ? k4rt::wrapAngleZeroToTwoPi(tc - thetaHit) : k4rt::wrapAngleZeroToTwoPi(thetaHit - tc);
}

//
// IPStateMode = SolenoidCrossing: AtIP for a track measured outside the solenoid.
//
// The solenoid boundary is a cylinder at R = R_sol closed by two disks at
// |z| = Z_sol. Going backward from the innermost hit, the outer circle is followed
// to the boundary crossing (barrel or end disk); there an inner helix with the
// inner field and the same pT and tanLambda is attached, and its perigee is stored.
//
bool CombinatorialTrackFinder::analyticalInnerPropagation(double cx_o, double cy_o, double R_outer, bool clockwise,
                                                          double pT, double tanLambda, double x_ref, double y_ref,
                                                          double z_ref, bool isEndcap,
                                                          edm4hep::MutableTrack& finalTrack) const {

  const double R_sol = m_solenoidBoundaryMm / 10.0; // mm -> cm
  const double Z_sol = m_solenoidHalfZmm / 10.0;    // mm -> cm

  double d_oc = std::sqrt(cx_o * cx_o + cy_o * cy_o); // outer center dist from IP (cm)

  // Backward direction on the outer circle (+1: increasing angle, for a clockwise track);
  // the same sign gives the reversed curvature of the inner helix.
  const double inner_sign = clockwise ? 1.0 : -1.0;

  const double theta_hit = std::atan2(y_ref - cy_o, x_ref - cx_o);

  // --- Solenoid boundary crossing ---
  double xc_sol, yc_sol, z_sol;

  if (!isEndcap) {
    // --- barrel: outer circle x solenoid cylinder ---
    if (d_oc < 1e-6) {
      return false;
    }
    if (R_outer <= std::abs(d_oc - R_sol) || R_outer >= d_oc + R_sol) {
      return false;
    }

    double P1x, P1y, P2x, P2y;
    circleCylinderCrossings(cx_o, cy_o, R_outer, d_oc, R_sol, P1x, P1y, P2x, P2y);

    // the crossing reached first going backward from the hit
    double arc1 = backwardArcAngle(cx_o, cy_o, theta_hit, clockwise, P1x, P1y);
    double arc2 = backwardArcAngle(cx_o, cy_o, theta_hit, clockwise, P2x, P2y);
    double arc_back;
    if (arc1 <= arc2) {
      xc_sol = P1x;
      yc_sol = P1y;
      arc_back = arc1;
    } else {
      xc_sol = P2x;
      yc_sol = P2y;
      arc_back = arc2;
    }

    double s_bwd = R_outer * arc_back;
    z_sol = z_ref - tanLambda * s_bwd;

    // The crossing must be on the IP side of the hit in z; if not, try the other one.
    if (std::abs(z_sol) > std::abs(z_ref) + 1e-3) {
      double arc_alt, s_bwd_alt, z_sol_alt;
      double xc_alt, yc_alt;
      if (arc1 <= arc2) {
        xc_alt = P2x;
        yc_alt = P2y;
        arc_alt = arc2;
      } else {
        xc_alt = P1x;
        yc_alt = P1y;
        arc_alt = arc1;
      }
      s_bwd_alt = R_outer * arc_alt;
      z_sol_alt = z_ref - tanLambda * s_bwd_alt;
      if (std::abs(z_sol_alt) > std::abs(z_ref) + 1e-3) {
        return false;
      }
      xc_sol = xc_alt;
      yc_sol = yc_alt;
      s_bwd = s_bwd_alt;
      z_sol = z_sol_alt;
    }

  } else {
    // --- end disk: arc s backward from the hit with z(s) = +-Z_sol ---
    if (std::abs(tanLambda) < 1e-6) {
      return false;
    }

    double z_cap = (z_ref < 0.0) ? -Z_sol : +Z_sol;

    // z_ref - tanLambda * s = z_cap
    double s_endcap = (z_ref - z_cap) / tanLambda;
    if (s_endcap < 0.0) {
      // hit already at |z| < Z_sol: use the barrel crossing instead
      if (d_oc < 1e-6 || R_outer <= std::abs(d_oc - R_sol) || R_outer >= d_oc + R_sol) {
        return false;
      }
      double P1x, P1y, P2x, P2y;
      circleCylinderCrossings(cx_o, cy_o, R_outer, d_oc, R_sol, P1x, P1y, P2x, P2y);
      double arc1 = backwardArcAngle(cx_o, cy_o, theta_hit, clockwise, P1x, P1y),
             arc2 = backwardArcAngle(cx_o, cy_o, theta_hit, clockwise, P2x, P2y);
      double s_bwd_fb;
      if (arc1 <= arc2) {
        xc_sol = P1x;
        yc_sol = P1y;
        s_bwd_fb = arc1 * R_outer;
      } else {
        xc_sol = P2x;
        yc_sol = P2y;
        s_bwd_fb = arc2 * R_outer;
      }
      z_sol = z_ref - tanLambda * s_bwd_fb;
      if (std::abs(z_sol) > std::abs(z_ref) + 1e-3) {
        return false;
      }
    } else {
      if (s_endcap > 2.0 * M_PI * R_outer) {
        return false;
      }

      // the point on the outer circle reached after the arc s_endcap backward from the hit
      double theta_ref = std::atan2(y_ref - cy_o, x_ref - cx_o);
      double d_theta = s_endcap / R_outer;
      double theta_cap = theta_ref + inner_sign * d_theta;

      xc_sol = cx_o + R_outer * std::cos(theta_cap);
      yc_sol = cy_o + R_outer * std::sin(theta_cap);
      z_sol = z_cap;

    } // end if (s_endcap < 0.0)
  } // end if (!isEndcap)

  // --- outward radial direction at the crossing ---
  double rx = (xc_sol - cx_o) / R_outer;
  double ry = (yc_sol - cy_o) / R_outer;

  // --- inner field, at the origin ---
  double B_inner = m_field.magneticField(dd4hep::Position(0.0, 0.0, 0.0)).z() / dd4hep::tesla;
  if (std::abs(B_inner) < 1e-6) {
    return false;
  }

  // --- inner circle: same pT, field reversed, so the centre is on the other side ---
  // Both circles go through the crossing point with the same tangent; with the field reversed the
  // track bends the other way, so the inner centre lies on the outward radial line of the outer
  // circle, at the inner radius R_inner = pT / (0.3 B_inner) beyond the crossing.
  double R_inner = k4rt::radiusFromPt(B_inner, pT); // cm
  double cx_i = xc_sol + R_inner * rx;
  double cy_i = yc_sol + R_inner * ry;

  // --- inner helix at its perigee (closest approach to the beam line) ---
  // d0 = distance of the inner centre from the origin minus R_inner; the direction at the
  // perigee is perpendicular to the line from the origin to the centre (+-90 degrees).
  double d0_inner_cm = std::sqrt(cx_i * cx_i + cy_i * cy_i) - R_inner;
  double phi0_inner = k4rt::wrapAngle(std::atan2(cy_i, cx_i) - inner_sign * M_PI / 2.0);
  double omega_inner_mm = inner_sign / (R_inner * 10.0); // 1/mm, signed

  // --- z at the perigee: backward arc from the crossing, the shorter way (< pi) ---
  double theta_c = std::atan2(yc_sol - cy_i, xc_sol - cx_i);
  double theta_ip = std::atan2(-cy_i, -cx_i);
  double bwd_arc_inner = (omega_inner_mm > 0) ? k4rt::wrapAngleZeroToTwoPi(theta_c - theta_ip)
                                              : k4rt::wrapAngleZeroToTwoPi(theta_ip - theta_c);
  if (bwd_arc_inner > M_PI)
    bwd_arc_inner = 2.0 * M_PI - bwd_arc_inner;

  double s_sol_to_ip = R_inner * bwd_arc_inner;
  double z0_inner_cm = z_sol - tanLambda * s_sol_to_ip;

  // --- stored as AtIP ---
  edm4hep::TrackState atIPState = createTrackState(d0_inner_cm * 10.0, // cm -> mm
                                                   phi0_inner, omega_inner_mm,
                                                   z0_inner_cm * 10.0, // cm -> mm
                                                   tanLambda, edm4hep::TrackState::AtIP);
  finalTrack.addToTrackStates(atIPState);

  return true;
}

bool CombinatorialTrackFinder::fitCircleNHits(const std::vector<edm4hep::TrackerHitPlane>& hits, double& x0, double& y0,
                                              double& radius, double& chi2ndf, Eigen::Matrix3d& fitCov3x3,
                                              std::vector<double>& residualsCm) const {

  const int N = static_cast<int>(hits.size());
  if (N < 3) {
    warning() << "fitCircleNHits: need at least 3 hits, got " << N << endmsg;
    return false;
  }

  // positions (cm) and weights: one radial residual per hit, so one resolution (du)
  std::vector<Eigen::Vector2d> pts(N);
  std::vector<double> sigma(N), w(N);
  for (int i = 0; i < N; ++i) {
    const auto& pos = hits[i].getPosition();
    pts[i] = {pos[0] / 10.0, pos[1] / 10.0};
    double su = hits[i].getDu();
    if (su <= 0)
      su = m_sigmaHitDefault.value() * 10.0; // property is in cm, du is in mm
    sigma[i] = su / 10.0;                    // mm -> cm
    w[i] = 1.0 / (sigma[i] * sigma[i]);
  }

  // --- Step 1: algebraic weighted least-squares circle fit (a first guess) ---
  // A circle can be written x^2 + y^2 + D*x + E*y + F = 0, which is LINEAR in (D, E, F): every
  // hit gives one equation D*x + E*y + F = -(x^2 + y^2), and the over-determined system A p = b
  // is solved by least squares through the normal equations (A^T A) p = A^T b. Completing the
  // squares gives the centre (-D/2, -E/2) and the radius sqrt(D^2 + E^2 - 4F) / 2.
  Eigen::MatrixXd A(N, 3);
  Eigen::VectorXd b(N);
  for (int i = 0; i < N; ++i) {
    double x = pts[i].x(), y = pts[i].y();
    A(i, 0) = w[i] * x;
    A(i, 1) = w[i] * y;
    A(i, 2) = w[i] * 1.0;
    b(i) = -w[i] * (x * x + y * y);
  }
  Eigen::Vector3d sol = (A.transpose() * A).ldlt().solve(A.transpose() * b);
  double D = sol(0), E = sol(1), F = sol(2);
  x0 = -D / 2.0;
  y0 = -E / 2.0;
  double disc = D * D + E * E - 4 * F;
  if (disc <= 0) {
    return false;
  }
  radius = std::sqrt(disc) / 2.0;

  // --- Step 2: Gauss-Newton geometric refinement ---
  // The algebraic fit does not minimise the real distances of the hits from the circle. Here the
  // residual of hit i is r_i = (|p_i - c| - R) / sigma_i (its distance from the circle, in units
  // of its error), and the parameters (cx, cy, R) are improved by Gauss-Newton steps
  //   delta = -(J^T J)^-1 J^T r,   J = d r / d(cx, cy, R) = -(dx/dist, dy/dist, 1) / sigma
  // until the step is below 1e-6 cm (or after 10 steps).
  Eigen::Vector3d params(x0, y0, radius);
  const int maxIter = 10;
  const double tol = 1e-6;
  for (int iter = 0; iter < maxIter; ++iter) {
    double cx = params(0), cy = params(1), r = params(2);
    Eigen::MatrixXd J(N, 3);
    Eigen::VectorXd res(N);
    for (int i = 0; i < N; ++i) {
      double dx = pts[i].x() - cx, dy = pts[i].y() - cy;
      double dist = std::sqrt(dx * dx + dy * dy);
      double sqw = std::sqrt(w[i]);
      res(i) = sqw * (dist - r);
      if (dist > 1e-12) {
        J(i, 0) = -sqw * dx / dist;
        J(i, 1) = -sqw * dy / dist;
        J(i, 2) = -sqw;
      } else {
        J(i, 0) = J(i, 1) = J(i, 2) = 0.0;
      }
    }
    Eigen::Vector3d delta = -(J.transpose() * J).ldlt().solve(J.transpose() * res);
    params += delta;
    if (delta.norm() < tol)
      break;
  }
  x0 = params(0);
  y0 = params(1);
  radius = params(2);
  if (radius <= 0) {
    return false;
  }

  // --- Step 3: covariance of (x0, y0, R): (J^T J)^-1 with the Jacobian at the final parameters ---
  {
    Eigen::MatrixXd Jf(N, 3);
    for (int i = 0; i < N; ++i) {
      double dx = pts[i].x() - x0, dy = pts[i].y() - y0;
      double dist = std::sqrt(dx * dx + dy * dy);
      double sqw = std::sqrt(w[i]);
      if (dist > 1e-12) {
        Jf(i, 0) = -sqw * dx / dist;
        Jf(i, 1) = -sqw * dy / dist;
        Jf(i, 2) = -sqw;
      } else {
        Jf(i, 0) = Jf(i, 1) = Jf(i, 2) = 0.0;
      }
    }
    fitCov3x3 = (Jf.transpose() * Jf).ldlt().solve(Eigen::Matrix3d::Identity());
  }

  // --- Step 4: per-hit residuals and chi2/ndf (3 fitted parameters: ndf = N - 3, at least 1) ---
  residualsCm.resize(N);
  double chi2 = 0.0;
  for (int i = 0; i < N; ++i) {
    double dx = pts[i].x() - x0, dy = pts[i].y() - y0;
    double res = std::abs(std::sqrt(dx * dx + dy * dy) - radius);
    residualsCm[i] = res;
    chi2 += (res / sigma[i]) * (res / sigma[i]);
  }
  int ndf = N - 3;
  if (ndf < 1)
    ndf = 1;
  chi2ndf = chi2 / ndf;

  return true;
}
// ============================================================================
// SOLENOID BOUNDARY SCAN
// ============================================================================
double CombinatorialTrackFinder::findSolenoidBoundary(bool alongZ) const {
  // The field is followed outward from near the IP in 1 mm steps until it is no longer the
  // central solenoid field; a binary search then finds that point to 0.01 mm.
  //   alongZ = false: along +x in the plane z = 0 (the radius of the solenoid)
  //   alongZ = true:  along z at half the solenoid radius (its half-length); not on the beam
  //                   axis, where compensating solenoids can change the field near the IP
  // "No longer the central field": Bz has moved more than half-way from its value at the IP
  // towards zero or beyond. This catches a change of sign (return field) and a step to no
  // field at all (some field maps have a 0 T gap at the end of the coil).
  const std::string tag = alongZ ? "[SolScanZ] |Z|" : "[SolScan] R";
  const double stepMm = 1.0;
  const double startMm = 10.0;
  const double scanLimitMm = 10000.0; // 10 m: beyond any FCC-ee solenoid

  const double bzCentral = fieldAlongScan(startMm, alongZ);
  if (std::abs(bzCentral) < 1e-6) {
    warning() << "  " << tag << ": no field near the IP; the solenoid boundary is set to " << scanLimitMm << " mm"
              << endmsg;
    return scanLimitMm;
  }
  // 1. Coarse scan: the 1 mm interval [lowMm, highMm] where the central field ends
  double lowMm = -1.0;
  double highMm = -1.0;
  for (double distanceMm = startMm + stepMm; distanceMm <= scanLimitMm; distanceMm += stepMm) {
    if (outsideCentralField(distanceMm, alongZ, bzCentral)) {
      lowMm = distanceMm - stepMm;
      highMm = distanceMm;
      break;
    }
  }
  if (lowMm < 0.0) {
    warning() << "  " << tag << ": the central field does not end within " << scanLimitMm
              << " mm; this value is used as the boundary" << endmsg;
    return scanLimitMm;
  }

  // 2. Binary search inside the interval, down to 0.01 mm
  for (int iteration = 0; iteration < 40; ++iteration) {
    const double middleMm = 0.5 * (lowMm + highMm);
    if (outsideCentralField(middleMm, alongZ, bzCentral))
      highMm = middleMm;
    else
      lowMm = middleMm;
    if (highMm - lowMm < 0.01)
      break;
  }

  const double boundaryMm = 0.5 * (lowMm + highMm);
  return boundaryMm;
}

double CombinatorialTrackFinder::fieldAlongScan(double distanceMm, bool alongZ) const {
  const double distanceCm = distanceMm / 10.0; // DD4hep positions are in cm
  dd4hep::Position position(distanceCm, 0.0, 0.0);
  // the z scan runs at half the solenoid radius (found by the radial scan, which runs first)
  if (alongZ)
    position = dd4hep::Position(0.5 * m_solenoidBoundaryMm / 10.0, 0.0, distanceCm);
  return m_field.magneticField(position).z() / dd4hep::tesla;
}

bool CombinatorialTrackFinder::outsideCentralField(double distanceMm, bool alongZ, double bzCentral) const {
  return std::abs(fieldAlongScan(distanceMm, alongZ) - bzCentral) > 0.5 * std::abs(bzCentral);
}

bool CombinatorialTrackFinder::readoutHasField(const std::string& name) const {
  for (const auto& field : m_bitFieldCoder->fields())
    if (field.name() == name)
      return true;
  return false;
}

// ============================================================================
// SOLENOID REGION
// ============================================================================
bool CombinatorialTrackFinder::isInsideSolenoid(const Eigen::Vector3d& pos) const {
  double r = std::sqrt(pos.x() * pos.x() + pos.y() * pos.y());
  double z = std::abs(pos.z());
  return (r < m_solenoidBoundaryMm / 10.0) && (z < m_solenoidHalfZmm / 10.0);
}

// True when no hit of the track lies inside the solenoid.
bool CombinatorialTrackFinder::allHitsOutsideSolenoid(const std::vector<edm4hep::TrackerHitPlane>& trackHits) const {
  for (const auto& hit : trackHits) {
    const auto& pos = hit.getPosition();
    const Eigen::Vector3d posCm(pos[0] / 10.0, pos[1] / 10.0, pos[2] / 10.0);
    if (isInsideSolenoid(posCm))
      return false;
  }
  return true;
}
