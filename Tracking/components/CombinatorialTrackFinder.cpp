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
#include <chrono>
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
    auto hasField = [&](const std::string& n) {
      for (const auto& f : m_bitFieldCoder->fields())
        if (f.name() == n)
          return true;
      return false;
    };
    if (!hasField(m_layerFieldName.value())) {
      error() << "Readout '" << m_encodingStringParameter.value() << "' (" << desc << ") has no field '"
              << m_layerFieldName.value() << "'. Set LayerFieldName to match this readout." << endmsg;
      return StatusCode::FAILURE;
    }
    m_hasRegionField = hasField(m_regionFieldName.value());
    if (!m_hasRegionField) {
      warning() << "Readout '" << m_encodingStringParameter.value() << "' (" << desc << ") has no field '"
                << m_regionFieldName.value() << "'; treating every hit as barrel. Set RegionFieldName "
                << "(\"type\" for MuonSystemReadoutID, \"side\" for GlobalTrackerReadoutID)." << endmsg;
    }
    // When one readout spans several subdetectors, "system" goes into the layer key,
    // so that equally numbered layers of different subdetectors stay apart.
    // (automatic: on exactly when DetectorNames lists several detectors)
    const bool multiDet = m_detectorNames.value().size() > 1;
    m_useSystemInCompositeID = multiDet && hasField("system");
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

  // Magnetic field (printed at a few points for reference)
  m_field = m_detector->field();

  dd4hep::Position center(0, 0, 0);
  dd4hep::Direction bfield = m_field.magneticField(center);
  info() << "Magnetic field at detector center: (" << bfield.x() / dd4hep::tesla << ", " << bfield.y() / dd4hep::tesla
         << ", " << bfield.z() / dd4hep::tesla << ") Tesla" << endmsg;

  const double solenoidRadius = 230.0; // cm
  const std::vector<dd4hep::Position> testPoints = {
      // in cm
      dd4hep::Position(0, 0, 0),
      dd4hep::Position(100, 100, 100),
      dd4hep::Position(solenoidRadius, 0, 0),
      dd4hep::Position(0, solenoidRadius, 0),
      dd4hep::Position(solenoidRadius / sqrt(2), solenoidRadius / sqrt(2), 0),
      dd4hep::Position(249, 0, 0),
      dd4hep::Position(0, 251, 0),
      dd4hep::Position(0, 0, 249),
      dd4hep::Position(0, 0, 218),
      dd4hep::Position(0, 0, -217),
      dd4hep::Position(0, 500, 0),
      dd4hep::Position(460, 0, 0),
      dd4hep::Position(200, 200, 460),
      dd4hep::Position(460, 0, 460)};

  info() << "Magnetic field at different positions:" << endmsg;
  for (const auto& point : testPoints) {
    dd4hep::Direction bfieldVec = m_field.magneticField(point);
    info() << "  Position (" << point.x() << ", " << point.y() << ", " << point.z() << ") cm: B = ("
           << bfieldVec.x() / dd4hep::tesla << ", " << bfieldVec.y() / dd4hep::tesla << ", "
           << bfieldVec.z() / dd4hep::tesla << ") Tesla" << endmsg;
  }

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
    seedConfig.maxHitDist_cm = m_maxPairHitDist.value(); // the same distance as the exhaustive search
    seedConfig.maxCandPerLayer = m_maxCandPerLayer.value();
    seedConfig.minTrackHits = m_minTrackHits.value();
    seedConfig.shortTrackMaxHits = m_shortTrackMaxHits.value();
    seedConfig.shadowHitWindow_cm = m_shadowHitWindow.value();
    seedConfig.kinkMaxChi2Triplet = m_kinkMaxChi2Triplet.value();
    seedConfig.seedKinkMaxChi2NDF = m_seedKinkMaxChi2NDF.value();
    seedConfig.trackKinkMaxChi2NDF = m_trackKinkMaxChi2NDF.value();
    seedConfig.shortTrackMinPt = m_minComboPT.value(); // one minimum pT for both searches
    seedConfig.bField_T = m_field.magneticField(dd4hep::Position(0, 0, 0)).z() / dd4hep::tesla;
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
      m << "SeedAndExtend: B = " << seedConfig.bField_T << " T; layer x/X0 [%]:";
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

  m_statTotalEvents++;

  debug() << "\n" << std::string(80, '=') << endmsg;

  unsigned int eventNumber = 0;
  if (!headers.empty()) {
    eventNumber = headers[0].getEventNumber();
    debug() << "Processing Event #" << eventNumber << endmsg;
  } else {
    debug() << "Processing new event (unknown number)" << endmsg;
  }

  debug() << std::string(80, '-') << endmsg;
  debug() << "Processing " << hits.size() << " tracker hits" << endmsg;

  if (msgLevel(MSG::DEBUG)) {
    debug() << "Hit details:" << endmsg;
    for (size_t i = 0; i < hits.size(); ++i) {
      const auto& hit = hits[i];
      const auto& pos = hit.getPosition();
      debug() << "  Hit " << i << " at (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ") mm" << endmsg;
    }
  }

  if (hits.size() < 3) {
    warning() << "Not enough hits to create tracks. Need at least 3." << endmsg;
    debug() << std::string(80, '=') << "\n" << endmsg;
    return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks));
  }

  try {
    findTracks(&hits, trackCollection, outputHits, recoSimLinks, outputLinks);
  } catch (const std::exception& ex) {
    error() << "Exception during track finding: " << ex.what() << endmsg;
    debug() << std::string(80, '=') << "\n" << endmsg;
    return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks));
  }

  debug() << "Found " << trackCollection.size() << " tracks" << endmsg;

  // per-track summary with the matched MC particles (DEBUG only)
  if (msgLevel(MSG::DEBUG))
    printTrackSummary(trackCollection, recoSimLinks);

  debug() << std::string(80, '=') << "\n" << endmsg;

  return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks));
}

// DEBUG printout after each event: every track's state at the first hit and AtIP,
// next to the MC particle of its first hit.
void CombinatorialTrackFinder::printTrackSummary(
    const edm4hep::TrackCollection& trackCollection,
    const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks) const {
  // Per-track printout, with the matched MC particle when the sim links are there
  std::unordered_map<uint32_t, edm4hep::SimTrackerHit> summarySimMap;
  for (const auto& link : recoSimLinks)
    summarySimMap[link.getFrom().id().index] = link.getTo();

  // --- Per-track summary ---
  for (size_t i = 0; i < trackCollection.size(); ++i) {
    const auto& track = trackCollection[i];
    debug() << "--- Track " << (i + 1) << "  (" << track.trackerHits_size() << " hits) ---" << endmsg;

    // state at the first hit
    for (int j = 0; j < static_cast<int>(track.trackStates_size()); ++j) {
      auto ts = track.getTrackStates(j);
      if (ts.location != edm4hep::TrackState::AtFirstHit)
        continue;
      double B_outer = m_field.magneticField(dd4hep::Position(0, 0, 0)).z() / dd4hep::tesla;
      // field at the first hit
      if (track.trackerHits_size() > 0) {
        auto h0 = track.getTrackerHits(0);
        B_outer = m_field
                      .magneticField(dd4hep::Position(h0.getPosition()[0] / 10.0, h0.getPosition()[1] / 10.0,
                                                      h0.getPosition()[2] / 10.0))
                      .z() /
                  dd4hep::tesla;
      }
      double pT_out = (std::abs(ts.omega) > 1e-12) ? k4rt::ptFromOmega(B_outer, ts.omega) : 0.0;
      double eta_out = std::asinh(ts.tanLambda);
      debug() << "  [AtFirstHit]       pT=" << std::fixed << std::setprecision(2) << pT_out
              << " GeV  phi=" << std::setprecision(3) << ts.phi << " rad  eta=" << std::setprecision(3) << eta_out
              << "  d0=" << std::setprecision(2) << ts.D0 / 10.0 << " cm"
              << "  z0=" << ts.Z0 / 10.0 << " cm  (w.r.t. the first hit)" << endmsg;
      break;
    }

    // state at the IP, if stored (IPStateMode)
    bool hasIP = false;
    for (int j = 0; j < static_cast<int>(track.trackStates_size()); ++j) {
      auto ts = track.getTrackStates(j);
      if (ts.location != edm4hep::TrackState::AtIP)
        continue;
      double B_inner = m_field.magneticField(dd4hep::Position(0, 0, 0)).z() / dd4hep::tesla;
      double pT_ip = (std::abs(ts.omega) > 1e-12) ? k4rt::ptFromOmega(B_inner, ts.omega) : 0.0;
      double eta_ip = std::asinh(ts.tanLambda);
      debug() << "  [Inner/AtIP]       pT=" << std::fixed << std::setprecision(2) << pT_ip
              << " GeV  phi=" << std::setprecision(3) << ts.phi << " rad  eta=" << std::setprecision(3) << eta_ip
              << "  d0=" << std::setprecision(2) << ts.D0 / 10.0 << " cm"
              << "  z0=" << ts.Z0 / 10.0 << " cm" << endmsg;
      hasIP = true;
      break;
    }

    // MC particle with most hits on the track
    std::unordered_map<int, std::pair<edm4hep::MCParticle, int>> mcCount;
    for (int j = 0; j < static_cast<int>(track.trackerHits_size()); ++j) {
      auto hit = track.getTrackerHits(j);
      auto it = summarySimMap.find(hit.id().index);
      if (it == summarySimMap.end())
        continue;
      auto mc = it->second.getParticle();
      if (!mc.isAvailable())
        continue;
      mcCount[mc.id().index].first = mc;
      mcCount[mc.id().index].second++;
    }
    int bestN = 0;
    edm4hep::MCParticle bestMC;
    bool hasMC = false;
    for (const auto& [idx, pr] : mcCount)
      if (pr.second > bestN) {
        bestN = pr.second;
        bestMC = pr.first;
        hasMC = true;
      }

    if (hasMC) {
      const auto& mom = bestMC.getMomentum();
      double mcPT = std::sqrt(mom.x * mom.x + mom.y * mom.y);
      double mcEta = (mcPT > 0) ? std::asinh(mom.z / mcPT) : 0.0;
      double mcPhi = std::atan2(mom.y, mom.x);
      const auto& vtx = bestMC.getVertex();
      double mcD0 = (-vtx.x * std::sin(mcPhi) + vtx.y * std::cos(mcPhi)) / 10.0; // mm->cm
      double mcZ0 = vtx.z / 10.0;                                                // mm -> cm

      debug() << "  [Truth PDG=" << bestMC.getPDG() << " q=" << bestMC.getCharge() << "]        pT=" << std::fixed
              << std::setprecision(2) << mcPT << " GeV  phi=" << std::setprecision(3) << mcPhi
              << " rad  eta=" << std::setprecision(3) << mcEta << "  d0=" << std::setprecision(2) << mcD0 << " cm"
              << "  z0=" << mcZ0 << " cm" << endmsg;

      if (hasIP) {
        for (int j = 0; j < static_cast<int>(track.trackStates_size()); ++j) {
          auto ts = track.getTrackStates(j);
          if (ts.location != edm4hep::TrackState::AtIP)
            continue;
          double B_inner = m_field.magneticField(dd4hep::Position(0, 0, 0)).z() / dd4hep::tesla;
          double pT_ip = (std::abs(ts.omega) > 1e-12) ? k4rt::ptFromOmega(B_inner, ts.omega) : 0.0;
          double dphi = ts.phi - mcPhi;
          dphi = k4rt::wrapAngle(dphi);
          debug() << "  [Resolution AtIP]  dpT/pT=" << std::fixed << std::setprecision(1)
                  << (mcPT > 0 ? (pT_ip - mcPT) / mcPT * 100.0 : 0.0) << "%"
                  << "  dphi=" << std::setprecision(4) << dphi << " rad"
                  << "  deta=" << std::setprecision(4) << std::asinh(ts.tanLambda) - mcEta
                  << "  dd0=" << std::setprecision(2) << ts.D0 / 10.0 - mcD0 << " cm"
                  << "  dz0=" << ts.Z0 / 10.0 - mcZ0 << " cm" << endmsg;
          break;
        }
      }
    } else {
      debug() << "  [Truth] no MC link found" << endmsg;
    }
  }
}

// Prints the run statistics.
StatusCode CombinatorialTrackFinder::finalize() {
  int nEvt = m_statTotalEvents.load();
  int nTrk = m_statTracksReconstructed.load();
  int nPos = m_statPositiveCharge.load();
  int nNeg = m_statNegativeCharge.load();
  int nCombos = m_statTotalTripletCombos.load();
  int nValid = m_statValidTriplets.load();
  int nChi2Rej = m_statChi2Rejected.load();
  int nFitFail = m_statCircleFitFailed.load();
  int nIsoRej = m_statIsolationRejected.load();
  int nConsecPhiRej = m_statConsecPhiRejected.load();
  int nPhiSpreadRej = m_statPhiSpreadRejected.load();
  int nConsecDistRej = m_statConsecDistRejected.load();
  int nPairDistRej = m_statPairDistRejected.load();
  int nEdepRej = m_statEdepRejected.load();
  int nOutlier = m_statOutlierHitsRemoved.load();
  int nNeighbRej = m_statNeighbourRejected.load();
  int nPTRej = m_statComboPTRejected.load();
  int nMinPTRej = m_statComboMinPTRejected.load();
  int nTanLRej = m_statTanLambdaRejected.load();
  int nClusters = m_statClustersBuilt.load();
  int nH3 = m_statNhit3.load();
  int nH4 = m_statNhit4.load();
  int nH5 = m_statNhit5.load();
  int nH6 = m_statNhit6.load();
  int nH7p = m_statNhit7plus.load();
  int nSL = m_statStateAtLastHit.load();
  int nSIP = m_statStateAtIP.load();
  int nAnalyt = m_statAnalyticalPropSuccess.load();
  double trkPerEvt = (nEvt > 0) ? double(nTrk) / double(nEvt) : 0.0;
  double comboEff = (nCombos > 0) ? 100.0 * double(nValid) / double(nCombos) : 0.0;
  double analFrac = (nTrk > 0) ? 100.0 * double(nAnalyt) / double(nTrk) : 0.0;
  double avgOutlier = (nValid > 0) ? double(nOutlier) / double(nValid) : 0.0;
  // run summary
  info() << "Run summary" << endmsg;
  info() << "  Events processed              : " << nEvt << endmsg;
  info() << "  Tracks reconstructed          : " << nTrk << endmsg;
  info() << "  Tracks per event              : " << std::fixed << std::setprecision(2) << trkPerEvt << endmsg;
  info() << "  Positive tracks (omega > 0)   : " << nPos << endmsg;
  info() << "  Negative tracks (omega < 0)   : " << nNeg << endmsg;
  info() << "Exhaustive search" << endmsg;
  info() << "  Combinations evaluated        : " << nCombos << endmsg;
  info() << "  Tracks found (best combo)     : " << nValid << "  (" << std::setprecision(1) << comboEff
         << "% of combinations)" << endmsg;
  info() << "  Circle fit failed             : " << nFitFail << endmsg;
  info() << "  Rejected by chi2/ndf cut      : " << nChi2Rej << endmsg;
  info() << "  Rejected by consecutive dphi  : " << nConsecPhiRej << endmsg;
  info() << "  Rejected by phi spread        : " << nPhiSpreadRej << endmsg;
  info() << "  Rejected by consecutive dist  : " << nConsecDistRej << endmsg;
  info() << "  Rejected by pair distance     : " << nPairDistRej << endmsg;
  info() << "  Rejected by isolation cut     : " << nIsoRej << endmsg;
  info() << "  Hits removed by edep cut      : " << nEdepRej << endmsg;
  info() << "  Outlier hits removed          : " << nOutlier << "  (" << std::setprecision(2) << avgOutlier
         << " per track)" << endmsg;
  info() << "  Rejected by neighbour gate    : " << nNeighbRej << endmsg;
  info() << "  Rejected by MaxComboPT        : " << nPTRej << endmsg;
  info() << "  Rejected by MinComboPT        : " << nMinPTRej << endmsg;
  info() << "  Rejected by tanLambda cut     : " << nTanLRej << endmsg;
  info() << "  Hit clusters built            : " << nClusters << endmsg;
  info() << "Hits per track" << endmsg;
  info() << "  3 hits                        : " << nH3 << endmsg;
  info() << "  4 hits                        : " << nH4 << endmsg;
  info() << "  5 hits                        : " << nH5 << endmsg;
  info() << "  6 hits                        : " << nH6 << endmsg;
  info() << "  7 or more hits                : " << nH7p << endmsg;
  info() << "Track states stored" << endmsg;
  info() << "  AtLastHit                     : " << nSL << endmsg;
  info() << "  AtIP                          : " << nSIP << endmsg;
  info() << "SolenoidCrossing AtIP state" << endmsg;
  info() << "  Built                         : " << nAnalyt << "  (" << std::setprecision(1) << analFrac
         << "% of tracks)" << endmsg;
  info() << "  Barrel: no solenoid crossing  : " << m_statAnalBarrelNoIsect.load() << endmsg;
  info() << "  Barrel: both z solutions fail : " << m_statAnalBarrelBothZfail.load() << endmsg;
  info() << "  Barrel: z solutions swapped   : " << m_statAnalBarrelZswapped.load() << endmsg;
  info() << "  Endcap: |tanLambda| too small : " << m_statAnalEndcapSmallTanL.load() << endmsg;
  info() << "  Endcap: already past the cap  : " << m_statAnalEndcapPastCap.load() << endmsg;
  info() << "  Endcap: arc longer than 2 pi R: " << m_statAnalEndcapArcLong.load() << endmsg;
  info() << "Pair-coincidence filter" << endmsg;
  info() << "  Hits without a partner        : " << m_statHitsNonCoinc.load() << endmsg;
  info() << "  Combos with < 2 doublets      : " << m_statCombosRejByCoinc.load() << endmsg;
  info() << "  Tracks recovered by pass 2    : " << m_statCoincPass2Tracks.load() << endmsg;

  return Algorithm::finalize();
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

void CombinatorialTrackFinder::findTracks(const edm4hep::TrackerHitPlaneCollection* hits,
                                          edm4hep::TrackCollection& finalTracks,
                                          edm4hep::TrackerHitPlaneCollection& outputHits,
                                          const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks,
                                          edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinks) const {

  if (hits->size() < 3) {
    debug() << "Not enough hits for tracking. Need at least 3, found " << hits->size() << endmsg;
    return;
  }

  if (msgLevel(MSG::DEBUG) && hits->size() > 0) {
    const auto& dbgHit = (*hits)[0];
  }

  // Input hit -> all SimTrackerHits that contributed to it (one digi hit can have
  // several). Keyed on (collectionID, index): merged input collections all start
  // their indices at 0.
  auto hitKey = [](const auto& h) -> uint64_t {
    const auto id = h.id();
    return (static_cast<uint64_t>(id.collectionID) << 32) | static_cast<uint32_t>(id.index);
  };
  std::unordered_map<uint64_t, std::vector<edm4hep::SimTrackerHit>> recoToSimMap;
  for (const auto& link : recoSimLinks) {
    auto recoHit = link.getFrom();
    recoToSimMap[hitKey(recoHit)].push_back(link.getTo());
  }

  // DEBUG: the MC particles behind this event's hits
  if (msgLevel(MSG::DEBUG)) {
    debug() << "MC truth particles for this event (from RecoSim links):" << endmsg;
    std::unordered_set<int> seenMCIdx;
    int nLinked = 0;
    for (size_t gi = 0; gi < hits->size(); ++gi) {
      auto it = recoToSimMap.find(hitKey((*hits)[gi]));
      if (it == recoToSimMap.end())
        continue;
      nLinked++;
      for (const auto& sim : it->second) {
        auto mcPart = sim.getParticle();
        if (!mcPart.isAvailable())
          continue;
        int mcIdx = mcPart.id().index;
        if (!seenMCIdx.insert(mcIdx).second)
          continue; // already printed
        const auto& mom = mcPart.getMomentum();
        const auto& vtx = mcPart.getVertex();
        double mcPt = std::sqrt(mom.x * mom.x + mom.y * mom.y);
        double mcP = std::sqrt(mcPt * mcPt + mom.z * mom.z);
        double mcEta = (mcPt > 0.0) ? std::asinh(mom.z / mcPt) : 0.0;
        double mcPhi = std::atan2(mom.y, mom.x);
        double vtxR = std::sqrt(vtx.x * vtx.x + vtx.y * vtx.y) / 10.0; // mm->cm
        debug() << "  [MC truth] PDG=" << mcPart.getPDG() << "  q=" << mcPart.getCharge() << "  pT=" << std::fixed
                << std::setprecision(3) << mcPt << " GeV/c"
                << "  pz=" << mom.z << " GeV/c"
                << "  |p|=" << mcP << " GeV/c"
                << "  phi=" << mcPhi << " rad"
                << "  eta=" << mcEta << "  vtx_R=" << vtxR << " cm"
                << "  vtx_z=" << vtx.z / 10.0 << " cm" << endmsg;
      }
    }
    debug() << "  => " << nLinked << "/" << hits->size() << " hits linked to sim, " << seenMCIdx.size()
            << " unique MC particle(s)" << endmsg;
  }

  // Helper lambda: given a newly created output hit whose data was copied from
  // an input hit, propagate every RecoSim link of the input hit to the output
  // link collection. For cluster digi this preserves multi-SimHit links
  // (muon + delta-rays in the same readout cell) instead of collapsing them
  // to one; baseline digi yields exactly one link per hit as before.
  auto propagateLink = [&](const edm4hep::MutableTrackerHitPlane& outHit, const edm4hep::TrackerHitPlane& inHit) {
    auto it = recoToSimMap.find(hitKey(inHit));
    if (it == recoToSimMap.end())
      return;
    for (const auto& sim : it->second) {
      auto outLink = outputLinks.create();
      outLink.setFrom(outHit);
      outLink.setTo(sim);
      outLink.setWeight(1.0f);
    }
  };

  // hits already used by accepted tracks
  std::vector<bool> globalUsedHits(hits->size(), false);

  // SeedAndExtend: all candidates are built in the first iteration, then taken one per
  // iteration, best first, skipping those that share a hit with an accepted track.
  std::vector<k4rt::SECandidate> seCands;
  bool seBuilt = false;
  size_t seNext = 0;
  // hit positions (mm) of the accepted tracks, for the neighbour-track cuts
  std::vector<std::array<double, 3>> acceptedTrackHitPositions;
  int trackNumber = 0;

  // One track per iteration, until no candidate is left
  while (true) {
    trackNumber++;
    debug() << "\n=== TRACK FINDING ITERATION " << trackNumber << " ===" << endmsg;

    int remainingHits = 0;
    std::map<int, int> unusedHitsPerLayer;
    for (size_t i = 0; i < globalUsedHits.size(); ++i) {
      if (!globalUsedHits[i]) {
        remainingHits++;
        int layerID = getLayerID((*hits)[i].getCellID());
        unusedHitsPerLayer[layerID]++;
      }
    }

    debug() << "Remaining unused hits: " << remainingHits << endmsg;
    if (msgLevel(MSG::DEBUG)) {
      for (const auto& [layer, count] : unusedHitsPerLayer) {
        debug() << "  Layer " << layer << ": " << count << " unused hits" << endmsg;
      }
    }

    if (remainingHits < 3) {
      debug() << "Not enough unused hits remaining (" << remainingHits << ") to form another track" << endmsg;
      break;
    }

    // unused hits of this iteration
    using HitInfo = k4rt::HitInfo;

    std::vector<HitInfo> allHitInfo;
    std::map<int, std::vector<size_t>> hitIndicesByCompositeLayer;

    allHitInfo.reserve(remainingHits);

    // MaxHitEdepKeV: skip hits with a large energy deposit (EDM4hep edep is in GeV)
    const double edepThreshGeV = (m_maxHitEdepKeV.value() > 0.0) ? m_maxHitEdepKeV.value() * 1e-6 // keV -> GeV
                                                                 : -1.0;                          // disabled
    int evtEdepRej = 0;
    for (size_t i = 0; i < hits->size(); ++i) {
      if (globalUsedHits[i])
        continue;

      const auto& hit = (*hits)[i];

      if (edepThreshGeV > 0.0 && hit.getEDep() > edepThreshGeV) {
        debug() << "  Hit " << i << " excluded by edep cut: " << hit.getEDep() * 1e6 << " keV > "
                << m_maxHitEdepKeV.value() << " keV" << endmsg;
        ++evtEdepRej;
        continue;
      }

      const dd4hep::rec::Surface* surface = findSurface(hit);

      if (surface) {
        int layerID = getLayerID(hit.getCellID());
        int typeID = getTypeID(hit.getCellID());
        int compositeID = getCompositeID(hit.getCellID());

        allHitInfo.emplace_back(i, hit, layerID, typeID, compositeID);
        hitIndicesByCompositeLayer[compositeID].push_back(allHitInfo.size() - 1);

        debug() << "Available hit " << i << ": CellID=" << hit.getCellID() << ", Layer=" << layerID
                << ", Type=" << typeID << ", Composite=" << compositeID << endmsg;
      } else {
        debug() << "Could not find surface for unused hit with cellID " << hit.getCellID() << endmsg;
      }
    }
    m_statEdepRejected += evtEdepRej;
    if (evtEdepRej > 0)
      debug() << "  edep pre-filter removed " << evtEdepRej << " hits this iteration" << endmsg;

    if (allHitInfo.size() < 3) {
      debug() << "Not enough unused hits with valid surfaces (" << allHitInfo.size() << ") for iteration "
              << trackNumber << endmsg;
      break;
    }

    debug() << "Available hits for iteration " << trackNumber << ": " << allHitInfo.size() << " in "
            << hitIndicesByCompositeLayer.size() << " detector different layers/regions" << endmsg;

    std::vector<int> compositeLayerIDs;
    for (const auto& [compositeID, indices] : hitIndicesByCompositeLayer) {
      compositeLayerIDs.push_back(compositeID);
    }

    // layer keys sorted by layer number (innermost first), then by key
    auto layerNumber = [](int compositeID) { return k4rt::keyWithoutRegion(compositeID); };
    std::sort(compositeLayerIDs.begin(), compositeLayerIDs.end(), [&](int a, int b) {
      int la = layerNumber(a);
      int lb = layerNumber(b);
      if (la != lb)
        return la < lb; // smaller layer number => nearer to IP
      return a < b;     // tie-break by compositeID numeric order
    });

    if (compositeLayerIDs.size() < 3) {
      debug() << "Not enough detector regions with unused hits for triplet seeding" << endmsg;
      break;
    }

    // Field at the centroid of the hits, to turn MaxComboPT into a radius cap for the
    // k = 3 score: the radius (cm) of a track with pT = MaxComboPT
    // (there are at least 3 hits here: fewer layers stopped the search above)
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
      bRefForK3 = m_field.magneticField(centroid).z() / dd4hep::tesla;
    }
    const double absBRefForK3 = std::max(std::abs(bRefForK3), 1e-6);
    const double rCapK3 = k4rt::radiusFromPt(absBRefForK3, m_maxComboPT.value());

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
      debug() << "  No single-hit compositeIDs - using median edep as baseline: " << meanEdepSingle << " keV" << endmsg;
    }

    std::unordered_map<size_t, double> hitEdepScore;
    hitEdepScore.reserve(allHitInfo.size());
    if (m_edepScoreWeight.value() > 0.0) {
      for (size_t idx = 0; idx < allHitInfo.size(); ++idx) {
        const double edep_ = hitEdepKeV.at(idx);
        const double lowDev_ = std::max(0.0, meanEdepSingle - edep_);
        const double highDev_ = std::max(0.0, edep_ - 2.0 * meanEdepSingle);
        hitEdepScore[idx] = lowDev_ + highDev_;
        debug() << "  EdepScore hit[" << idx << "] compositeID=" << allHitInfo[idx].compositeID << " edep=" << edep_
                << " keV  baseline=" << meanEdepSingle << " keV"
                << " score=" << hitEdepScore[idx] << " keV" << endmsg;
      }

      const int maxHPC = m_maxHitsPerComposite.value();
      if (maxHPC > 0) {
        for (auto& [cid, cidIndices] : hitIndicesByCompositeLayer) {
          if (static_cast<int>(cidIndices.size()) <= maxHPC)
            continue;
          std::sort(cidIndices.begin(), cidIndices.end(),
                    [&](size_t a, size_t b) { return hitEdepScore[a] < hitEdepScore[b]; });
          cidIndices.resize(maxHPC);
          debug() << "  compositeID=" << cid << " trimmed to top " << maxHPC << " hits by edep similarity" << endmsg;
        }
      }
    }

    // hits per layer number (barrel and endcaps together)
    std::map<int, std::vector<size_t>> hitsByLayerNumber;
    for (size_t idx = 0; idx < allHitInfo.size(); ++idx) {
      int lyr = allHitInfo[idx].layerID;
      hitsByLayerNumber[lyr].push_back(idx);
    }

    std::vector<bool> usedHits(hits->size(), false);
    for (size_t i = 0; i < globalUsedHits.size(); ++i) {
      usedHits[i] = globalUsedHits[i];
    }

    struct TrackCandidate {
      edm4hep::Track track;
      double pT;
      std::vector<size_t> hitIndices;
      std::vector<int> usedCompositeIDs;

      TrackCandidate(edm4hep::Track t, double pt, const std::vector<size_t>& indices,
                     const std::vector<int>& composites)
          : track(t), pT(pt), hitIndices(indices), usedCompositeIDs(composites) {}
    };

    std::vector<TrackCandidate> trackCandidates;

    edm4hep::TrackCollection candidateTracks; // this iteration only

    // --- N-hit combinatorial search (Strategy = Exhaustive) ---
    // All C(N,k) combinations for k = MaxCombinatorialHits down to MinTrackHits,
    // each scored by the chi2/ndf of its circle fit, after road, isolation and
    // quality cuts; outliers are then removed from the best combination.

    // hit indices ordered by layer, inner -> outer
    std::vector<size_t> orderedHII;
    orderedHII.reserve(allHitInfo.size());
    for (const auto& cid : compositeLayerIDs) {
      for (size_t idx : hitIndicesByCompositeLayer[cid])
        orderedHII.push_back(idx);
    }

    const int nAvail = static_cast<int>(orderedHII.size());
    const int maxK = std::min(nAvail, m_maxCombinatorialHits.value());
    const int minK = std::max(3, m_minTrackHits.value());

    using CombResult = k4rt::CombResult;
    CombResult bestCombo;
    CombResult fallbackCombo; // best over-threshold combo; promoted if no valid combo exists
    bool foundCombination = false;

    // per-iteration counters (debug printout)
    int evtCombosEval = 0;    // combinations that reached the fit stage
    int evtIsoRej = 0;        // rejected by isolation cut
    int evtChi2Rej = 0;       // rejected by chi2/ndf threshold
    int evtFitFail = 0;       // fit itself returned false
    int evtConsecPhiRej = 0;  // rejected by MaxConsecDeltaPhi
    int evtPhiSpreadRej = 0;  // rejected by MaxComboPhiSpread
    int evtConsecDistRej = 0; // rejected by MaxConsecutiveHitDistance
    int evtPairDistRej = 0;   // rejected by MaxPairHitDistance
    int evtTanLambdaRej = 0;  // rejected by the tanLambda-deviation cut

    // --- Quality cuts ---

    // hit index -> layer key
    std::unordered_map<size_t, int> gIdxToComposite;
    gIdxToComposite.reserve(allHitInfo.size());
    for (const auto& hi : allHitInfo) {
      gIdxToComposite[hi.index] = hi.compositeID;
    }

    // Isolation: each hit of the combination must be HitIsolationCut (cm) away from
    // every other unused hit in a different layer. Hits in the same layer (secondaries
    // next to the track) are not counted.
    auto allIsolated = [&](const std::vector<size_t>& gIdxVec) -> bool {
      const double cut = m_hitIsolationCut.value();
      if (cut <= 0.0)
        return true;
      for (size_t gi : gIdxVec) {
        const auto refPos = (*hits)[gi].getPosition();
        Eigen::Vector3d ref(refPos[0] / 10.0, refPos[1] / 10.0, refPos[2] / 10.0);
        int giComposite = -999;
        auto cgit = gIdxToComposite.find(gi);
        if (cgit != gIdxToComposite.end())
          giComposite = cgit->second;
        for (size_t j = 0; j < hits->size(); ++j) {
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
          int jComposite = -999;
          auto cjit = gIdxToComposite.find(j);
          if (cjit == gIdxToComposite.end())
            continue; // not in this iteration's hit pool
          jComposite = cjit->second;
          if (jComposite == giComposite)
            continue;
          const auto p = (*hits)[j].getPosition();
          Eigen::Vector3d cand(p[0] / 10.0, p[1] / 10.0, p[2] / 10.0);
          double dist = (cand - ref).norm();
          if (dist < cut) {
            debug() << "    isolation FAIL: selected hit gIdx=" << gi << " compositeID=" << giComposite
                    << " too close to neighbour gIdx=" << j << " compositeID=" << jComposite << "  dist=" << dist
                    << " cm  cut=" << cut << " cm" << endmsg;
            return false;
          }
        }
      }
      return true;
    };

    // Road cuts: consecutive-layer dphi and total phi spread across the combo.
    // Both cuts reject combinations whose hits span an azimuthal window too wide
    // to belong to a single track.  Hits in gIdxVec are ordered inner->outer
    // (slot order matches compositeLayerIDs which is sorted by layer number).
    // Hit-set predicates (road cuts) live in ExhaustiveSearch.
    // (the debug stream is fetched per call: a cached one keeps the level of the last message)
    // the last value: largest deviation of a consecutive-pair tanLambda (dz / xy chord) from
    // the median tanLambda of the combination (0.5 = about 27 degrees of dip angle)
    const double maxComboTanLambdaDeviation = 0.5;
    const k4rt::RoadCutConfig roadCfg{m_maxConsecDeltaPhi.value(), m_maxComboPhiSpread.value(),
                                      m_maxConsecHitDist.value(), m_maxPairHitDist.value(), maxComboTanLambdaDeviation};
    auto roadCtx = [&]() {
      return k4rt::RoadCutContext{*hits, gIdxToComposite, roadCfg, msgLevel(MSG::DEBUG) ? &debug() : nullptr,
                                  &m_cidRegion};
    };
    auto passesRoadCuts = [&](const std::vector<size_t>& gIdxVec) -> bool {
      return k4rt::passesRoadCuts(gIdxVec, roadCtx());
    };
    auto passesDistanceCuts = [&](const std::vector<size_t>& gIdxVec) -> bool {
      return k4rt::passesDistanceCuts(gIdxVec, roadCtx());
    };
    auto passesTanLambdaCut = [&](const std::vector<size_t>& gIdxVec) -> bool {
      return k4rt::passesTanLambdaCut(gIdxVec, roadCtx());
    };

    // --- Pair coincidence (RequirePairCoincidence), for detectors built of layer ---
    // doublets (0,1), (2,3), (4,5): two hits of a doublet are coincident when
    // |dphi| < 0.15 rad, |deta| < 0.10 and their distance is below PairCoincidenceMaxDist.
    constexpr double kPairTolPhi = 0.15;
    constexpr double kPairTolEta = 0.10;
    const double kPairMaxDistMM = m_pairCoincMaxDist.value() * 10.0; // cm -> mm (hit positions are in mm)

    auto coincidentHits = [&](const edm4hep::TrackerHitPlane& a, const edm4hep::TrackerHitPlane& b) -> bool {
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
    };

    // partner layer of a doublet: 0 <-> 1, 2 <-> 3, 4 <-> 5
    auto pairPartnerCid = [](int composite) -> int {
      int absC = std::abs(composite);
      int layer = k4rt::layerOfKey(composite);
      int partnerLayer = layer ^ 1;
      int absPartner = absC - layer + partnerLayer; // same region and subdetector
      return (composite < 0) ? -absPartner : absPartner;
    };

    // Does hit at allHitInfo[hii] have any coincident partner in its paired layer?
    auto hasPartnerInPair = [&](size_t hii) -> bool {
      int partnerCid = pairPartnerCid(allHitInfo[hii].compositeID);
      auto it = hitIndicesByCompositeLayer.find(partnerCid);
      if (it == hitIndicesByCompositeLayer.end())
        return false;
      for (size_t partnerHii : it->second) {
        if (coincidentHits(allHitInfo[hii].hit, allHitInfo[partnerHii].hit))
          return true;
      }
      return false;
    };

    // Count complete doublets in a combo: a doublet is complete when the combo
    // contains hits in both layers of one of (0,1)/(2,3)/(4,5) AND they are
    // mutually coincident, both in the same region.
    auto countCompleteDoublets = [&](const std::vector<size_t>& gIdxVec) -> int {
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
            int regionB = k4rt::regionOfKey(cidB);
            if (regionA != regionB)
              continue;
            if (coincidentHits((*hits)[gi_a], (*hits)[gi_b])) {
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
    };

    // unused hits of each layer: hitsForCompAll[i] belongs to compositeLayerIDs[i]
    const int nDistinctAll = static_cast<int>(compositeLayerIDs.size());
    std::vector<std::vector<size_t>> hitsForCompAll(nDistinctAll);
    for (int ci = 0; ci < nDistinctAll; ++ci) {
      int cid = compositeLayerIDs[ci];
      for (size_t hii : hitIndicesByCompositeLayer.at(cid)) {
        if (!usedHits[allHitInfo[hii].index]) {
          hitsForCompAll[ci].push_back(hii);
        }
      }
    }

    // Pair coincidence: first pass only with hits that have a coincident partner
    // (hitsForCompSeed)
    const bool pcHard = m_requirePairCoinc.value();
    std::vector<std::vector<size_t>> hitsForCompSeed;
    if (pcHard) {
      hitsForCompSeed.assign(nDistinctAll, std::vector<size_t>{});
      int evtNonCoinc = 0;
      for (int ci = 0; ci < nDistinctAll; ++ci) {
        for (size_t hii : hitsForCompAll[ci]) {
          if (hasPartnerInPair(hii)) {
            hitsForCompSeed[ci].push_back(hii);
          } else {
            ++evtNonCoinc;
          }
        }
      }
      m_statHitsNonCoinc += evtNonCoinc;
      debug() << "  PairCoinc: " << evtNonCoinc << " hits demoted (no partner in paired layer)" << endmsg;
    }

    // --- Hit clustering (DoHitClustering) ---
    // Hits of neighbouring layers that pass the road cuts are joined (union-find);
    // each cluster is searched on its own. Without clustering: one cluster, all hits.
    struct ClusterDef {
      std::vector<int> cids;
      std::vector<std::vector<size_t>> hitsForComp;
    };
    std::vector<ClusterDef> activeClusters;

    if (m_doHitClustering.value() && nDistinctAll >= 2) {
      const size_t nhii = allHitInfo.size();
      std::vector<size_t> ufParent(nhii);
      for (size_t i = 0; i < nhii; ++i)
        ufParent[i] = i;

      auto ufFind = [&](size_t x) -> size_t {
        while (ufParent[x] != x) {
          ufParent[x] = ufParent[ufParent[x]]; // path halving
          x = ufParent[x];
        }
        return x;
      };
      auto ufUnite = [&](size_t a, size_t b) {
        a = ufFind(a);
        b = ufFind(b);
        if (a != b)
          ufParent[a] = b;
      };

      const double phiCut = m_maxConsecDeltaPhi.value();
      const double distCut = m_maxConsecHitDist.value(); // cm

      // join hits of neighbouring layers
      for (int ci = 0; ci < nDistinctAll - 1; ++ci) {
        for (size_t hii_a : hitsForCompAll[ci]) {
          const auto pa = (*hits)[allHitInfo[hii_a].index].getPosition();
          const double phiA = std::atan2(pa[1], pa[0]);
          for (size_t hii_b : hitsForCompAll[ci + 1]) {
            const auto pb = (*hits)[allHitInfo[hii_b].index].getPosition();
            if (phiCut > 0.0) {
              double dPhi = std::atan2(pb[1], pb[0]) - phiA;
              dPhi = k4rt::wrapAngle(dPhi);
              if (std::abs(dPhi) > phiCut)
                continue;
            }
            if (distCut > 0.0) {
              double dx = pb[0] - pa[0], dy = pb[1] - pa[1], dz = pb[2] - pa[2];
              if (std::sqrt(dx * dx + dy * dy + dz * dz) * 0.1 > distCut)
                continue;
            }
            ufUnite(hii_a, hii_b);
          }
        }
      }

      // also across one missing layer, with doubled cuts
      const double phiCut2 = (phiCut > 0.0) ? phiCut * 2.0 : -1.0;
      const double distCut2 = (distCut > 0.0) ? distCut * 2.0 : -1.0;
      for (int ci = 0; ci < nDistinctAll - 2; ++ci) {
        for (size_t hii_a : hitsForCompAll[ci]) {
          const auto pa = (*hits)[allHitInfo[hii_a].index].getPosition();
          const double phiA = std::atan2(pa[1], pa[0]);
          for (size_t hii_b : hitsForCompAll[ci + 2]) {
            const auto pb = (*hits)[allHitInfo[hii_b].index].getPosition();
            if (phiCut2 > 0.0) {
              double dPhi = std::atan2(pb[1], pb[0]) - phiA;
              dPhi = k4rt::wrapAngle(dPhi);
              if (std::abs(dPhi) > phiCut2)
                continue;
            }
            if (distCut2 > 0.0) {
              double dx = pb[0] - pa[0], dy = pb[1] - pa[1], dz = pb[2] - pa[2];
              if (std::sqrt(dx * dx + dy * dy + dz * dz) * 0.1 > distCut2)
                continue;
            }
            ufUnite(hii_a, hii_b);
          }
        }
      }

      std::unordered_map<size_t, std::vector<size_t>> rootToHIIs;
      for (int ci = 0; ci < nDistinctAll; ++ci) {
        for (size_t hii : hitsForCompAll[ci]) {
          rootToHIIs[ufFind(hii)].push_back(hii);
        }
      }

      // keep clusters with at least minK layers, layers in inner -> outer order
      for (auto& [root, hiis] : rootToHIIs) {
        std::unordered_map<int, std::vector<size_t>> cidMap;
        for (size_t hii : hiis)
          cidMap[allHitInfo[hii].compositeID].push_back(hii);
        if (static_cast<int>(cidMap.size()) < minK) {
          debug() << "  cluster root=" << root << " spans only " << cidMap.size()
                  << " layer(s) - skipped (minK=" << minK << ")" << endmsg;
          continue;
        }
        ClusterDef cd;
        for (int cid : compositeLayerIDs) {
          auto it = cidMap.find(cid);
          if (it == cidMap.end())
            continue;
          cd.cids.push_back(cid);
          cd.hitsForComp.push_back(std::move(it->second));
        }
        debug() << "  cluster root=" << root << " layers=" << cd.cids.size() << " hits=" << hiis.size() << endmsg;
        activeClusters.push_back(std::move(cd));
      }

      m_statClustersBuilt += static_cast<int>(activeClusters.size());
      debug() << "  Hit clustering [trk " << trackNumber << "]: " << nDistinctAll << " layer(s) -> "
              << activeClusters.size() << " cluster(s)" << endmsg;
    } else {
      ClusterDef cd;
      cd.cids = compositeLayerIDs;
      cd.hitsForComp = hitsForCompAll;
      activeClusters.push_back(std::move(cd));
    }

    // hard pair-coincidence mode: the first pass uses only the coincident hits
    std::vector<ClusterDef> activeClustersSeed;
    if (pcHard) {
      ClusterDef cd;
      cd.cids = compositeLayerIDs;
      cd.hitsForComp = hitsForCompSeed;
      activeClustersSeed.push_back(std::move(cd));
    }

    // passes: hard pair coincidence = coincident hits first, then all hits if the best
    // combination has < 3 complete doublets; otherwise one pass over all hits
    const int nPairPasses = pcHard ? 2 : 1;

    // --- SeedAndExtend: the next candidate whose hits are all free ---
    // Everything after the choice of hits is shared with the exhaustive search.
    if (m_useSeedExtend) {
      if (!seBuilt) {
        std::vector<std::pair<size_t, int>> pool;
        pool.reserve(allHitInfo.size());
        for (const auto& hi : allHitInfo)
          pool.emplace_back(hi.index, hi.compositeID);
        k4rt::SeedStats ss;
        const auto t0 = std::chrono::steady_clock::now();
        seCands = m_seedFinder->run(*hits, pool, ss);
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        seBuilt = true;
        if (m_useSeedExtend)
          debug() << "SeedAndExtend: " << ss.extended << " seeds extended, " << ss.skippedClaimed
                  << " skipped (hits claimed), " << ss.rejectedTracks << " rejected (" << ss.rejectedHoles
                  << " for holes, " << ss.rejectedShortPt << " short with low pT); confirmation +" << ss.confirmAdded
                  << " hits, " << ss.confirmSwapped << " swapped, " << ss.confirmRemoved << " removed; " << ss.merged
                  << " pieces merged, " << ss.loopers << " loopers, " << ss.shadowClaimed << " shadow hits claimed"
                  << endmsg;
        debug() << "SeedAndExtend: " << pool.size() << " hits, " << ss.rawTriplets << " raw combinations, " << ss.seeds
                << " seeds (" << ss.seedsAfterZ << " after z check, " << ss.extended << " extended) -> "
                << ss.candidates << " candidates in " << ms << " ms" << endmsg;
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
        bestCombo = CombResult{};
        bestCombo.gIdxVec = c.gIdx;
        for (size_t g : c.gIdx)
          bestCombo.hiiVec.push_back(gToHii[g]);
        bestCombo.chi2ndf = c.chi2ndf;
        bestCombo.avgProxCm = 0.0;
        bestCombo.x0 = c.x0;
        bestCombo.y0 = c.y0;
        bestCombo.radius = c.R;
        bestCombo.fitCov = c.cov;
        bestCombo.residuals = c.residuals;
        bestCombo.nHits = static_cast<int>(c.gIdx.size());
        bestCombo.kinkChi2 = c.kinkChi2;
        bestCombo.kinkNdf = c.kinkNdf;
        bestCombo.looper = c.looper;
        foundCombination = true;
        break;
      }
    }

    for (int pairPassIdx = 0; !m_useSeedExtend && pairPassIdx < nPairPasses; ++pairPassIdx) {
      // second pass only if the first one's best has fewer than 3 complete doublets
      if (pairPassIdx == 1 && foundCombination) {
        int nDoub = countCompleteDoublets(bestCombo.gIdxVec);
        if (nDoub >= 3)
          break;
        debug() << "  PairCoinc: pass-1 best has " << nDoub << " complete doublets - running pass-2 over all hits"
                << endmsg;
      }
      bool foundBeforeThisPass = foundCombination;
      int doubletsBeforeThisPass = foundCombination ? countCompleteDoublets(bestCombo.gIdxVec) : -1;
      const std::vector<ClusterDef>& clustersThisPass =
          (pcHard && pairPassIdx == 0) ? activeClustersSeed : activeClusters;

      // clusterLayerIDs and hitsForComp are this cluster's layers and hits
      for (const auto& cluster : clustersThisPass) {
        const std::vector<int>& clusterLayerIDs = cluster.cids;
        const std::vector<std::vector<size_t>>& hitsForComp = cluster.hitsForComp; // NOLINT(shadow)
        const int nDistinct = static_cast<int>(clusterLayerIDs.size());
        const int maxKeff = std::min(maxK, nDistinct);

        // For k = maxKeff ... minK: every choice of k layers, and for each every
        // combination of one hit per chosen layer.
        for (int k = maxKeff; k >= minK; --k) {
          // Enumerate C(nDistinct, k) - indices into clusterLayerIDs
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
              // every combination of one hit per layer
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
                  auto comboLabel = [&]() -> std::string {
                    std::string s = "[";
                    for (int ii = 0; ii < k; ++ii) {
                      if (ii)
                        s += ",";
                      s += std::to_string(gIdxVec[ii]);
                    }
                    return s + "]";
                  };

                  if (!passesRoadCuts(gIdxVec)) {
                    // which road cut failed (for the counters): repeat the
                    // consecutive-dphi cut, with the same exceptions as passesRoadCuts
                    bool consecFail = false;
                    {
                      const int nH = static_cast<int>(gIdxVec.size());
                      double sumR = 0.0, sumAbsZ = 0.0;
                      for (int i = 0; i < nH; ++i) {
                        const auto p = (*hits)[gIdxVec[i]].getPosition();
                        sumR += std::sqrt(p[0] * p[0] + p[1] * p[1]);
                        sumAbsZ += std::abs(p[2]);
                      }
                      const double thetaCombo = std::atan2(sumR / nH, sumAbsZ / nH);
                      constexpr double kNearBeamLo = 15.0 * M_PI / 180.0;
                      constexpr double kNearBeamHi = 160.0 * M_PI / 180.0;
                      const bool nearBeamPipe = (thetaCombo < kNearBeamLo || thetaCombo > kNearBeamHi);
                      const double maxConsec = m_maxConsecDeltaPhi.value();
                      if (!nearBeamPipe && maxConsec > 0.0) {
                        for (int ri = 0; ri < nH - 1 && !consecFail; ++ri) {
                          int cidA = 0, cidB = 0;
                          auto itA = gIdxToComposite.find(gIdxVec[ri]);
                          auto itB = gIdxToComposite.find(gIdxVec[ri + 1]);
                          if (itA != gIdxToComposite.end())
                            cidA = itA->second;
                          if (itB != gIdxToComposite.end())
                            cidB = itB->second;
                          const bool aBarrel = (k4rt::regionOfCid(cidA, &m_cidRegion) == 0);
                          const bool bBarrel = (k4rt::regionOfCid(cidB, &m_cidRegion) == 0);
                          if (aBarrel != bBarrel)
                            continue; // barrel-endcap step: no cut
                          const auto pa = (*hits)[gIdxVec[ri]].getPosition();
                          const auto pb = (*hits)[gIdxVec[ri + 1]].getPosition();
                          double dPhi = std::atan2(pb[1], pb[0]) - std::atan2(pa[1], pa[0]);
                          dPhi = k4rt::wrapAngle(dPhi);
                          if (std::abs(dPhi) > maxConsec)
                            consecFail = true;
                        }
                      }
                    }
                    if (consecFail) {
                      m_statConsecPhiRejected++;
                      evtConsecPhiRej++;
                      debug() << "  k=" << k << " " << comboLabel() << " -> REJECT consec dphi" << endmsg;
                    } else {
                      m_statPhiSpreadRejected++;
                      evtPhiSpreadRej++;
                      debug() << "  k=" << k << " " << comboLabel() << " -> REJECT phi spread" << endmsg;
                    }
                  } else if (!passesDistanceCuts(gIdxVec)) {
                    // which distance cut failed (for the counters)
                    bool consecFail = false;
                    {
                      const double mc = m_maxConsecHitDist.value();
                      if (mc > 0.0) {
                        const int nH = static_cast<int>(gIdxVec.size());
                        for (int ri = 0; ri < nH - 1 && !consecFail; ++ri) {
                          const auto pa = (*hits)[gIdxVec[ri]].getPosition();
                          const auto pb = (*hits)[gIdxVec[ri + 1]].getPosition();
                          double dx = pb[0] - pa[0], dy = pb[1] - pa[1], dz = pb[2] - pa[2];
                          if (std::sqrt(dx * dx + dy * dy + dz * dz) * 0.1 > mc)
                            consecFail = true;
                        }
                      }
                    }
                    if (consecFail) {
                      m_statConsecDistRejected++;
                      evtConsecDistRej++;
                      debug() << "  k=" << k << " " << comboLabel() << " -> REJECT consec dist" << endmsg;
                    } else {
                      m_statPairDistRejected++;
                      evtPairDistRej++;
                      debug() << "  k=" << k << " " << comboLabel() << " -> REJECT pair dist" << endmsg;
                    }
                  } else if (!allIsolated(gIdxVec)) {
                    m_statIsolationRejected++;
                    evtIsoRej++;
                    debug() << "  k=" << k << " " << comboLabel() << " -> REJECT isolation" << endmsg;
                  } else if (!passesTanLambdaCut(gIdxVec)) {
                    m_statTanLambdaRejected++;
                    evtTanLambdaRej++;
                    debug() << "  k=" << k << " " << comboLabel() << " -> REJECT tanLambda dev" << endmsg;
                  } else if (pcHard && countCompleteDoublets(gIdxVec) < 2) {
                    // hard pair coincidence: at least 2 complete doublets
                    m_statCombosRejByCoinc++;
                    debug() << "  k=" << k << " " << comboLabel() << " -> REJECT <2 complete doublets" << endmsg;
                  } else {
                    m_statTotalTripletCombos++;
                    evtCombosEval++;
                    std::vector<edm4hep::TrackerHitPlane> comboHits;
                    comboHits.reserve(k);
                    for (size_t gi : gIdxVec)
                      comboHits.push_back((*hits)[gi]);

                    // circle parameters start at 0, so the pT below is defined even if the fit fails
                    double cx0 = 0.0, cy0 = 0.0, cr = 0.0, cchi2ndf = 0.0;
                    Eigen::Matrix3d cCov = Eigen::Matrix3d::Zero();
                    std::vector<double> cRes;
                    bool fitOK = false;
                    try {
                      fitOK = fitCircleNHits(comboHits, cx0, cy0, cr, cchi2ndf, cCov, cRes);
                    } catch (...) {
                    }

                    // pT of the combination, from its circle radius
                    const double comboPT = k4rt::ptFromRadius(absBRefForK3, cr);
                    if (!fitOK) {
                      evtFitFail++;
                      m_statCircleFitFailed++;
                      debug() << "  k=" << k << " " << comboLabel() << " -> REJECT fit failed" << endmsg;
                    } else if (cchi2ndf >= m_nHitMaxChi2NDF.value()) {
                      evtChi2Rej++;
                      m_statChi2Rejected++;
                      debug() << "  k=" << k << " " << comboLabel() << " pT=" << std::fixed << std::setprecision(1)
                              << comboPT << " GeV"
                              << " chi2/ndf=" << cchi2ndf << " -> REJECT chi2 threshold=" << m_nHitMaxChi2NDF.value()
                              << endmsg;
                    } else if (comboPT > m_maxComboPT.value()) {
                      // pT above MaxComboPT: not a candidate, but kept as fallback in
                      // case every combination is above it (near-straight tracks)
                      m_statComboPTRejected++;
                      debug() << "  k=" << k << " " << comboLabel() << " -> REJECT pT=" << comboPT
                              << " > MaxComboPT=" << m_maxComboPT.value() << " (stored as fallback)" << endmsg;
                      // best fallback: more hits; then lower chi2/ndf (k >= 4) or
                      // lower pT (k = 3)
                      bool isBetterFallback = false;
                      if (k != fallbackCombo.nHits) {
                        isBetterFallback = (k > fallbackCombo.nHits);
                      } else if (k == 3) {
                        isBetterFallback = (cr < fallbackCombo.radius);
                      } else {
                        isBetterFallback = (cchi2ndf < fallbackCombo.chi2ndf);
                      }
                      if (isBetterFallback) {
                        fallbackCombo.hiiVec = hiiVec;
                        fallbackCombo.gIdxVec = gIdxVec;
                        fallbackCombo.chi2ndf = cchi2ndf;
                        fallbackCombo.avgProxCm = 0.0;
                        fallbackCombo.x0 = cx0;
                        fallbackCombo.y0 = cy0;
                        fallbackCombo.radius = cr;
                        fallbackCombo.fitCov = cCov;
                        fallbackCombo.residuals = cRes;
                        fallbackCombo.nHits = k;
                      }
                    } else if (m_minComboPT.value() > 0.0 && comboPT < m_minComboPT.value()) {
                      // pT below MinComboPT: rejected (no fallback)
                      m_statComboMinPTRejected++;
                      debug() << "  k=" << k << " " << comboLabel() << " -> REJECT pT=" << comboPT
                              << " < MinComboPT=" << m_minComboPT.value() << endmsg;
                    } else {
                      // Ranking of the combinations (lower score is better):
                      //  1. more hits always wins;
                      //  2. k >= 4: chi2/ndf + w * edep score / EdepNormKeV;
                      //  3. k = 3 (chi2 is always 0 for 3 points):
                      //     w * edep score / EdepNormKeV - min(R, R_cap) / R_cap,
                      //     i.e. straighter (higher pT) and MIP-like is better.
                      // The edep score is the mean hit score against the event
                      // baseline (see above).
                      double cEdepDev = 0.0;
                      const double w = m_edepScoreWeight.value();
                      if (w > 0.0 && !hitEdepScore.empty()) {
                        for (size_t hii_i : hiiVec)
                          cEdepDev += hitEdepScore.at(hii_i);
                        cEdepDev /= static_cast<double>(k);
                      }
                      const double edepNorm = std::max(m_edepNormKeV.value(), 1e-9);
                      bool isBetter;
                      if (k != bestCombo.nHits) {
                        isBetter = (k > bestCombo.nHits);
                      } else if (k == 3) {
                        // k = 3 (radius term saturates at R_cap; w = 0: radius only)
                        if (w > 0.0) {
                          const double cRadNorm = std::min(cr, rCapK3) / rCapK3;
                          const double bRadNorm = std::min(bestCombo.radius, rCapK3) / rCapK3;
                          const double cScore3 = w * cEdepDev / edepNorm - cRadNorm;
                          const double bScore3 = w * bestCombo.avgProxCm / edepNorm - bRadNorm;
                          isBetter = (cScore3 < bScore3);
                        } else {
                          isBetter = (cr > bestCombo.radius);
                        }
                      } else {
                        const double cScore = cchi2ndf + (w > 0.0 ? w * cEdepDev / edepNorm : 0.0);
                        const double bScore = bestCombo.chi2ndf + (w > 0.0 ? w * bestCombo.avgProxCm / edepNorm : 0.0);
                        isBetter = (cScore < bScore);
                      }
                      if (isBetter) {
                        bestCombo.hiiVec = hiiVec;
                        bestCombo.gIdxVec = gIdxVec;
                        bestCombo.chi2ndf = cchi2ndf;
                        bestCombo.avgProxCm = cEdepDev;
                        bestCombo.x0 = cx0;
                        bestCombo.y0 = cy0;
                        bestCombo.radius = cr;
                        bestCombo.fitCov = cCov;
                        bestCombo.residuals = cRes;
                        bestCombo.nHits = k;
                        foundCombination = true;
                        debug() << "  k=" << k << " " << comboLabel() << " pT=" << std::fixed << std::setprecision(1)
                                << comboPT << " GeV"
                                << " chi2/ndf=" << cchi2ndf << " -> NEW BEST" << endmsg;
                      } else {
                        debug() << "  k=" << k << " " << comboLabel() << " pT=" << std::fixed << std::setprecision(1)
                                << comboPT << " GeV"
                                << " chi2/ndf=" << cchi2ndf << " -> not better (best k=" << bestCombo.nHits
                                << " chi2/ndf=" << bestCombo.chi2ndf << ")" << endmsg;
                      }
                    }
                  }
                }

                // next combination of hits
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

            // next choice of k layers
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
      } // end cluster loop
      // run statistics: tracks found or improved by the second pass
      if (pairPassIdx == 1 && foundCombination) {
        int nDoubAfter = countCompleteDoublets(bestCombo.gIdxVec);
        if (!foundBeforeThisPass) {
          m_statCoincPass2Tracks++;
        } else if (nDoubAfter > doubletsBeforeThisPass) {
          m_statCoincPass2Tracks++;
        }
      }
    } // end pair-coincidence two-pass

    const int maxKeffAll = std::min(maxK, nDistinctAll);
    debug() << "  N-hit search [trk " << trackNumber << "]: " << nAvail << " avail hits in " << nDistinctAll
            << " regions, k=" << maxKeffAll << "->" << minK << " | evaluated=" << evtCombosEval
            << " consecPhiRej=" << evtConsecPhiRej << " phiSpreadRej=" << evtPhiSpreadRej
            << " consecDistRej=" << evtConsecDistRej << " pairDistRej=" << evtPairDistRej << " isoRej=" << evtIsoRej
            << " tanLRej=" << evtTanLambdaRej << " chi2Rej=" << evtChi2Rej << " fitFail=" << evtFitFail
            << " edepRej=" << evtEdepRej
            << (foundCombination ? (" | BEST k=" + std::to_string(bestCombo.nHits) +
                                    " chi2/ndf=" + std::to_string(bestCombo.chi2ndf) +
                                    " edepDev=" + std::to_string(bestCombo.avgProxCm) + " keV")
                                 : " | NO COMBINATION FOUND")
            << endmsg;

    if (!foundCombination) {
      if (fallbackCombo.nHits > 0) {
        // every combination was above MaxComboPT: use the fallback
        bestCombo = fallbackCombo;
        foundCombination = true;
        debug() << "  Promoted fallback combo: k=" << bestCombo.nHits << " chi2/ndf=" << bestCombo.chi2ndf
                << " pT=" << std::fixed << std::setprecision(1) << k4rt::ptFromRadius(absBRefForK3, bestCombo.radius)
                << " GeV (all combos exceeded MaxComboPT=" << m_maxComboPT.value() << " GeV)" << endmsg;
      } else {
        debug() << "No valid N-hit combination found for track " << trackNumber << " - stopping" << endmsg;
        break;
      }
    }

    // --- Guided hit choice in crowded layers (doCrowdedLayerHitSelection) ---
    // Layers with one hit are "clean". In each crowded layer the hit is chosen by
    //   A (>= 3 clean layers): distance to the circle through the clean hits,
    //      score 0.75 * dr / SigmaHitDefault + 0.25 * edep score / EdepNormKeV;
    //   B (2 clean layers): kink angle to the line through the two clean hits,
    //      score 0.75 * kink / MaxConsecDeltaPhi + 0.25 * edep score / EdepNormKeV.
    // The result replaces the best combination if its refit has a lower chi2/ndf.
    if (m_doCrowdedLayerHitSelection && !m_useSeedExtend) {
      int nCleanCIDs = 0;
      int nCrowdedCIDs = 0;
      for (const auto& [cid, cidHIIs] : hitIndicesByCompositeLayer) {
        if (static_cast<int>(cidHIIs.size()) == 1)
          ++nCleanCIDs;
        else if (static_cast<int>(cidHIIs.size()) > 1)
          ++nCrowdedCIDs;
      }
      const int totalCIDs = nCleanCIDs + nCrowdedCIDs;

      // clean hits, ordered by layer key
      auto collectCleanHits = [&](std::vector<edm4hep::TrackerHitPlane>& cleanHits, std::vector<size_t>& cleanHIIs,
                                  std::vector<size_t>& cleanGIdxs) {
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
      };

      // refit the guided choice; keep it if it is different and better
      auto applyGuidedResult = [&](std::vector<edm4hep::TrackerHitPlane>& guidedHits, std::vector<size_t>& guidedHIIs,
                                   std::vector<size_t>& guidedGIdxs, const char* branchName) {
        double gx0 = 0, gy0 = 0, gr = 0, gchi2 = 0;
        Eigen::Matrix3d gCov;
        std::vector<double> gRes;
        bool guidedOK = false;
        try {
          guidedOK = fitCircleNHits(guidedHits, gx0, gy0, gr, gchi2, gCov, gRes);
        } catch (...) {
        }

        if (!guidedOK) {
          debug() << branchName << " fit failed - keeping general combo" << endmsg;
          return;
        }
        if (gchi2 > m_nHitMaxChi2NDF.value()) {
          debug() << branchName << " fit chi2/ndf=" << gchi2 << " > threshold - keeping general combo" << endmsg;
          return;
        }
        bool differentChoice = false;
        for (size_t gHII : guidedHIIs) {
          bool inGeneral = false;
          for (size_t bHII : bestCombo.hiiVec)
            if (bHII == gHII) {
              inGeneral = true;
              break;
            }
          if (!inGeneral) {
            differentChoice = true;
            break;
          }
        }
        if (differentChoice && gchi2 < bestCombo.chi2ndf) {
          debug() << branchName << ": " << (differentChoice ? "OVERRIDING" : "CONFIRMING") << " (" << nCrowdedCIDs
                  << " crowded CID(s))"
                  << " | guided chi2/ndf=" << gchi2 << " | general chi2/ndf=" << bestCombo.chi2ndf << endmsg;
          bestCombo.hiiVec = guidedHIIs;
          bestCombo.gIdxVec = guidedGIdxs;
          bestCombo.chi2ndf = gchi2;
          bestCombo.x0 = gx0;
          bestCombo.y0 = gy0;
          bestCombo.radius = gr;
          bestCombo.fitCov = gCov;
          bestCombo.residuals = gRes;
          bestCombo.nHits = static_cast<int>(guidedHIIs.size());
        } else {
          debug() << branchName << ": confirms general combo"
                  << " (all same hits, guided chi2/ndf=" << gchi2 << ")" << endmsg;
        }
      };

      // --- A: circle through the clean hits ---
      if (nCleanCIDs >= 3 && nCrowdedCIDs >= 1) {
        debug() << "Guided crowded-layer Branch A: " << nCleanCIDs << " clean CIDs + " << nCrowdedCIDs
                << " crowded CID(s)" << endmsg;

        std::vector<edm4hep::TrackerHitPlane> cleanHits;
        std::vector<size_t> cleanHIIs, cleanGIdxs;
        collectCleanHits(cleanHits, cleanHIIs, cleanGIdxs);

        if (static_cast<int>(cleanHits.size()) >= 3) {
          double cx0 = 0, cy0 = 0, cr = 0, cchi2 = 0;
          Eigen::Matrix3d cCov;
          std::vector<double> cRes;
          bool cleanOK = false;
          try {
            cleanOK = fitCircleNHits(cleanHits, cx0, cy0, cr, cchi2, cCov, cRes);
          } catch (...) {
          }

          if (cleanOK && cr > 0) {
            const double drNorm = std::max(m_sigmaHitDefault.value(), 1e-9); // cm
            const double edepNorm = std::max(m_edepNormKeV.value(), 1e-9);   // keV

            std::vector<edm4hep::TrackerHitPlane> guidedHits = cleanHits;
            std::vector<size_t> guidedHIIs = cleanHIIs;
            std::vector<size_t> guidedGIdxs = cleanGIdxs;
            bool allCrowdedResolved = true;

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
                debug() << "    BranchA CID=" << cid << " hii=" << hii << " dr=" << dr << " cm edepDev=" << edepDev
                        << " keV score=" << score << endmsg;
                if (score < bestScore) {
                  bestScore = score;
                  bestHII = hii;
                }
              }
              if (bestHII == SIZE_MAX) {
                allCrowdedResolved = false;
                break;
              }
              guidedHIIs.push_back(bestHII);
              guidedGIdxs.push_back(allHitInfo[bestHII].index);
              guidedHits.push_back(allHitInfo[bestHII].hit);
            }

            if (allCrowdedResolved)
              applyGuidedResult(guidedHits, guidedHIIs, guidedGIdxs, "Guided crowded-layer BranchA");
          }
        }
      }
      // --- B: line through the two clean hits ---
      else if (nCleanCIDs == 2 && nCrowdedCIDs >= 1 && totalCIDs >= 3) {
        debug() << "Guided crowded-layer Branch B: 2 clean CIDs + " << nCrowdedCIDs << " crowded CID(s)" << endmsg;

        std::vector<edm4hep::TrackerHitPlane> cleanHits;
        std::vector<size_t> cleanHIIs, cleanGIdxs;
        collectCleanHits(cleanHits, cleanHIIs, cleanGIdxs);

        if (static_cast<int>(cleanHits.size()) == 2) {
          // direction from the inner (A) to the outer (B) clean hit, cm
          const auto& posA_mm = cleanHits[0].getPosition();
          const auto& posB_mm = cleanHits[1].getPosition();
          double ax = posA_mm.x * 0.1, ay = posA_mm.y * 0.1, az = posA_mm.z * 0.1;
          double bx = posB_mm.x * 0.1, by = posB_mm.y * 0.1, bz = posB_mm.z * 0.1;
          double dabx = bx - ax, daby = by - ay, dabz = bz - az;
          double lenAB = std::sqrt(dabx * dabx + daby * daby + dabz * dabz);

          const double sigmaKappa = std::max(m_maxConsecDeltaPhi.value(), 1e-4); // rad
          const double edepNorm = std::max(m_edepNormKeV.value(), 1e-9);         // keV

          std::vector<edm4hep::TrackerHitPlane> guidedHits = cleanHits;
          std::vector<size_t> guidedHIIs = cleanHIIs;
          std::vector<size_t> guidedGIdxs = cleanGIdxs;
          bool allCrowdedResolved = (lenAB > 1e-9); // degenerate if A==B

          if (allCrowdedResolved) {
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
                debug() << "    BranchB CID=" << cid << " hii=" << hii << " kappa=" << kappa
                        << " rad edepDev=" << edepDev << " keV score=" << score << endmsg;
                if (score < bestScore) {
                  bestScore = score;
                  bestHII = hii;
                }
              }
              if (bestHII == SIZE_MAX) {
                allCrowdedResolved = false;
                break;
              }
              guidedHIIs.push_back(bestHII);
              guidedGIdxs.push_back(allHitInfo[bestHII].index);
              guidedHits.push_back(allHitInfo[bestHII].hit);
            }

            if (allCrowdedResolved)
              applyGuidedResult(guidedHits, guidedHIIs, guidedGIdxs, "Guided crowded-layer BranchB");
          }
        }
      }
    }

    // --- Neighbour-track cuts ---
    // A combination with a hit within NeighbourTrackMaxDist (cm) of an accepted track
    // may be a delta-ray shadow. It is kept only with >= NeighbourTrackMinHits hits
    // and chi2/ndf <= NeighbourTrackMaxChi2NDF; otherwise its hits near the accepted
    // track are marked as used and the search goes on.
    // Not applied in SeedAndExtend, whose tracks already own their hits.
    if (!m_useSeedExtend && m_neighbourTrackMaxDist.value() > 0.0 && !acceptedTrackHitPositions.empty()) {
      const double distThreshMm = m_neighbourTrackMaxDist.value() * 10.0; // cm -> mm
      const double distThreshSq = distThreshMm * distThreshMm;
      bool isNeighbour = false;
      for (size_t hii : bestCombo.hiiVec) {
        const auto& hp = allHitInfo[hii].hit.getPosition(); // mm
        for (const auto& ap : acceptedTrackHitPositions) {
          double dx = hp.x - ap[0], dy = hp.y - ap[1], dz = hp.z - ap[2];
          if (dx * dx + dy * dy + dz * dz < distThreshSq) {
            isNeighbour = true;
            break;
          }
        }
        if (isNeighbour)
          break;
      }

      if (isNeighbour) {
        const int nHits = static_cast<int>(bestCombo.hiiVec.size());
        const bool passHits = (nHits >= m_neighbourTrackMinHits.value());
        const bool passChi2 = (bestCombo.chi2ndf <= m_neighbourTrackMaxChi2NDF.value());
        if (!passHits || !passChi2) {
          debug() << "Neighbour-gate: iter " << trackNumber << " REJECTED shadow combo"
                  << " hits=" << nHits << " chi2/ndf=" << bestCombo.chi2ndf
                  << "  (required >=" << m_neighbourTrackMinHits.value()
                  << " hits AND <=" << m_neighbourTrackMaxChi2NDF.value() << " chi2/ndf)" << endmsg;
          ++m_statNeighbourRejected;
          // only the hits near the accepted track are used up
          int nConsumed = 0;
          for (size_t hii : bestCombo.hiiVec) {
            const auto& hp = allHitInfo[hii].hit.getPosition(); // mm
            bool hitIsNeighbour = false;
            for (const auto& ap : acceptedTrackHitPositions) {
              double dx = hp.x - ap[0], dy = hp.y - ap[1], dz = hp.z - ap[2];
              if (dx * dx + dy * dy + dz * dz < distThreshSq) {
                hitIsNeighbour = true;
                break;
              }
            }
            if (hitIsNeighbour) {
              globalUsedHits[allHitInfo[hii].index] = true;
              ++nConsumed;
            }
          }
          debug() << "  Neighbour-gate: consumed " << nConsumed << " of " << bestCombo.hiiVec.size()
                  << " combo hits (within " << m_neighbourTrackMaxDist.value() * 10.0 << " mm of accepted track)"
                  << endmsg;
          continue;
        } else {
          debug() << "Neighbour-gate: iter " << trackNumber << " passes stricter cuts"
                  << " (hits=" << nHits << ", chi2/ndf=" << bestCombo.chi2ndf << ")" << endmsg;
        }
      }
    }

    {
      debug() << "  Best combo: " << bestCombo.nHits << " hits, chi2/ndf=" << bestCombo.chi2ndf << endmsg;
      for (int bi = 0; bi < bestCombo.nHits; ++bi) {
        size_t gi = bestCombo.gIdxVec[bi];
        size_t hii = bestCombo.hiiVec[bi];
        const auto pos = (*hits)[gi].getPosition();
        double res = (bi < static_cast<int>(bestCombo.residuals.size())) ? bestCombo.residuals[bi] : -1.0;
        debug() << "    hit[" << bi << "] gIdx=" << gi << " compositeID=" << allHitInfo[hii].compositeID << " pos=("
                << std::fixed << std::setprecision(1) << pos[0] / 10.0 << "," << pos[1] / 10.0 << "," << pos[2] / 10.0
                << ") cm  residual=" << std::setprecision(4) << res << " cm" << endmsg;
      }
    }

    // --- Outlier removal: drop the worst hit while its residual is above ---
    //    OutlierSigma * SigmaHitDefault, and refit (not in SeedAndExtend)
    {
      const double outlierCutCm = m_outlierSigma.value() * m_sigmaHitDefault.value();
      for (int iter = 0; !m_useSeedExtend && iter < m_maxOutlierIterations.value(); ++iter) {
        if (static_cast<int>(bestCombo.gIdxVec.size()) <= minK)
          break;

        size_t worstIdx = 0;
        double worstRes = 0.0;
        for (size_t ri = 0; ri < bestCombo.residuals.size(); ++ri) {
          if (bestCombo.residuals[ri] > worstRes) {
            worstRes = bestCombo.residuals[ri];
            worstIdx = ri;
          }
        }
        if (worstRes <= outlierCutCm)
          break;

        {
          size_t wGi = bestCombo.gIdxVec[worstIdx];
          size_t wHii = bestCombo.hiiVec[worstIdx];
          const auto wPos = (*hits)[wGi].getPosition();
          debug() << "  Outlier iter " << iter + 1 << ": removing hit[" << worstIdx << "] gIdx=" << wGi
                  << " compositeID=" << allHitInfo[wHii].compositeID << " pos=(" << std::fixed << std::setprecision(1)
                  << wPos[0] / 10.0 << "," << wPos[1] / 10.0 << "," << wPos[2] / 10.0
                  << ") cm  residual=" << std::setprecision(4) << worstRes << " cm (cut=" << outlierCutCm << " cm)"
                  << endmsg;
        }
        m_statOutlierHitsRemoved++;

        bestCombo.gIdxVec.erase(bestCombo.gIdxVec.begin() + worstIdx);
        bestCombo.hiiVec.erase(bestCombo.hiiVec.begin() + worstIdx);
        bestCombo.residuals.erase(bestCombo.residuals.begin() + worstIdx);
        bestCombo.nHits = static_cast<int>(bestCombo.gIdxVec.size());

        if (bestCombo.nHits < minK)
          break;

        std::vector<edm4hep::TrackerHitPlane> refitHits;
        for (size_t gi : bestCombo.gIdxVec)
          refitHits.push_back((*hits)[gi]);
        double rx0, ry0, rr, rchi2;
        Eigen::Matrix3d rCov;
        std::vector<double> rRes;
        bool refitOK = false;
        try {
          refitOK = fitCircleNHits(refitHits, rx0, ry0, rr, rchi2, rCov, rRes);
        } catch (...) {
        }
        if (refitOK) {
          bestCombo.x0 = rx0;
          bestCombo.y0 = ry0;
          bestCombo.radius = rr;
          bestCombo.chi2ndf = rchi2;
          bestCombo.fitCov = rCov;
          bestCombo.residuals = rRes;
        }
      }
    }

    for (size_t gi : bestCombo.gIdxVec)
      usedHits[gi] = true;

    // --- pT of the chosen combination ---
    // from its circle radius and the field at the centroid of its hits (inside or outside
    // the solenoid); it orders the candidates and is used by the SolenoidCrossing AtIP state
    const int Nh = bestCombo.nHits;
    {
      double sumXf = 0, sumYf = 0, sumZf = 0;
      for (size_t gi : bestCombo.gIdxVec) {
        const auto p = (*hits)[gi].getPosition();
        sumXf += p[0] / 10.0;
        sumYf += p[1] / 10.0;
        sumZf += p[2] / 10.0;
      }
      dd4hep::Position avgPosN(sumXf / Nh, sumYf / Nh, sumZf / Nh);
      double bFieldN = m_field.magneticField(avgPosN).z() / dd4hep::tesla;
      double pTN = k4rt::ptFromRadius(bFieldN, bestCombo.radius); // GeV/c

      debug() << "N-hit circle fit: " << Nh << " hits  pT=" << pTN << " GeV/c"
              << "  chi2/ndf=" << bestCombo.chi2ndf << endmsg;
      bestCombo.pT = pTN;

      auto nhSeedTrack = candidateTracks.create(); // empty track (no states)
      m_statTotalTripletCombos++;
      m_statValidTriplets++;

      std::vector<int> usedComps;
      for (size_t hii : bestCombo.hiiVec)
        usedComps.push_back(allHitInfo[hii].compositeID);
      trackCandidates.emplace_back(nhSeedTrack, pTN, bestCombo.gIdxVec, usedComps);
    }

    // one candidate per iteration: the best combination
    if (trackCandidates.empty()) {
      debug() << "No track candidates found in iteration " << trackNumber << " - stopping track search" << endmsg;
      break;
    }

    std::sort(trackCandidates.begin(), trackCandidates.end(),
              [](const TrackCandidate& a, const TrackCandidate& b) { return a.pT > b.pT; });

    const auto& candidate = trackCandidates.front();

    for (size_t idx : candidate.hitIndices) {
      globalUsedHits[idx] = true;
    }

    bool validIndices = true;
    for (size_t idx : candidate.hitIndices) {
      if (idx >= usedHits.size()) {
        error() << "Final processing: hit index " << idx << " out of bounds (usedHits size: " << usedHits.size() << ")"
                << endmsg;
        validIndices = false;
        break;
      }
    }

    if (!validIndices) {
      continue;
    }

    debug() << "Processing track candidate " << trackNumber << " with pT=" << candidate.pT
            << " GeV/c  nHits=" << candidate.hitIndices.size() << endmsg;

    // the track's hits, by index (cellIDs need not be unique)
    std::vector<edm4hep::TrackerHitPlane> trackHits;
    trackHits.reserve(3);

    for (size_t idx : candidate.hitIndices) {
      if (idx < hits->size()) {
        trackHits.push_back((*hits)[idx]);
      } else {
        warning() << "Hit index " << idx << " out of bounds - skipping hit" << endmsg;
      }
    }

    if (trackHits.size() != candidate.hitIndices.size()) {
      warning() << "Could not retrieve all hits for track candidate " << trackNumber << " - expected "
                << candidate.hitIndices.size() << ", got " << trackHits.size() << endmsg;
      continue;
    }

    debug() << "Successfully built trackHits vector with " << trackHits.size() << " hits for track " << trackNumber
            << endmsg;

    debug() << "N-hit track " << trackNumber << " has " << trackHits.size()
            << " hits after combinatorial selection and outlier rejection" << endmsg;

    // for the neighbour-track cuts of the next iterations
    if (m_neighbourTrackMaxDist.value() > 0.0) {
      for (const auto& th : trackHits) {
        const auto& p = th.getPosition(); // mm
        acceptedTrackHitPositions.push_back({p.x, p.y, p.z});
      }
    }

    auto finalTrack = finalTracks.create();

    debug() << "Created new final track object for track " << trackNumber << endmsg;

    // The fitted helix as EDM4hep track states (makeEdm4hepStateFromCircle): at the first
    // hit, at the last hit and, with IPStateMode = Helix, at the closest approach to the origin.
    if (trackHits.size() >= 3) {
      const double x0 = bestCombo.x0; // circle centre and radius, cm
      const double y0 = bestCombo.y0;
      const double radius = bestCombo.radius;
      if (radius > 0) {
        debug() << "N-hit circle fit for track " << trackNumber << " (" << trackHits.size()
                << " hits): chi2/ndf=" << bestCombo.chi2ndf << "  centre (" << x0 << ", " << y0 << ") cm  radius "
                << radius << " cm  pT=" << bestCombo.pT << " GeV/c" << endmsg;
        struct StoredState {
          StatePoint where;
          int location;
          bool wanted;
        };
        const StoredState stored[] = {
            {StatePoint::FirstHit, edm4hep::TrackState::AtFirstHit, true},
            {StatePoint::LastHit, edm4hep::TrackState::AtLastHit, true},
            {StatePoint::Origin, edm4hep::TrackState::AtIP, m_ipStateMode.value() == "Helix"},
        };
        for (const auto& entry : stored) {
          if (!entry.wanted)
            continue;
          edm4hep::TrackState state;
          if (!makeEdm4hepStateFromCircle(x0, y0, radius, bestCombo.fitCov, trackHits, entry.where, state)) {
            warning() << "Could not build track state " << entry.location << " of track " << trackNumber
                      << " (degenerate hits)" << endmsg;
            continue;
          }
          state.location = entry.location;
          finalTrack.addToTrackStates(state);
        }
        // the circle fit has one radial residual per hit and 3 parameters
        const int fitNdf = std::max(1, static_cast<int>(trackHits.size()) - 3);
        finalTrack.setChi2(static_cast<float>(bestCombo.chi2ndf * fitNdf)); // raw chi2
        finalTrack.setNdf(fitNdf);
      } else {
        warning() << "Circle fit failed for track " << trackNumber << endmsg;
      }
    }

    // circle-fit ndf (n_hits - 3)
    finalTrack.setNdf(std::max(1, static_cast<int>(trackHits.size()) - 3));
    // SeedAndExtend: the stored chi2 is the kink chi2 (with multiple scattering);
    // loopers get track type 1
    if (bestCombo.kinkNdf > 0) {
      finalTrack.setChi2(static_cast<float>(bestCombo.kinkChi2));
      finalTrack.setNdf(bestCombo.kinkNdf);
    }
    if (bestCombo.looper)
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
    if (m_ipStateMode.value() == "SolenoidCrossing" && allHitsOutsideSolenoid(trackHits) && foundCombination &&
        haveFirstHitState && firstHitState.omega != 0.0f && finalTrack.trackerHits_size() > 0) {
      // outer circle (cm)
      double R_outer_anal = bestCombo.radius; // cm
      double cx_o_anal = bestCombo.x0;        // cm
      double cy_o_anal = bestCombo.y0;        // cm

      // the innermost hit anchors the crossing
      auto refHit = finalTrack.getTrackerHits(0);
      double xr = refHit.getPosition()[0] / 10.0; // mm -> cm
      double yr = refHit.getPosition()[1] / 10.0;
      double zr = refHit.getPosition()[2] / 10.0;

      // barrel hit: crossing of the solenoid barrel; endcap hit: of its end face
      bool isEndcapHit = (getTypeID(refHit.getCellID()) != 0);
      debug() << "  [Analytical] Ref hit typeID=" << getTypeID(refHit.getCellID()) << "  -> "
              << (isEndcapHit ? "endcap" : "barrel") << " crossing" << endmsg;

      // EDM4hep: omega > 0 means clockwise seen from +z
      const bool clockwise = (firstHitState.omega > 0.0f);
      bool atIPok = analyticalInnerPropagation(cx_o_anal, cy_o_anal, R_outer_anal, clockwise, bestCombo.pT,
                                               firstHitState.tanLambda, xr, yr, zr, isEndcapHit, finalTrack);

      if (atIPok) {
        m_statAnalyticalPropSuccess++;

        // truth printouts (DEBUG only)
        if (msgLevel(MSG::DEBUG)) {
          // printout: d0 of the AtIP state vs the true production radius
          try {
            edm4hep::MCParticle bestMCpart;
            bool hasMCpart = false;
            for (const auto& th : trackHits) {
              auto it = recoToSimMap.find(hitKey(th));
              if (it == recoToSimMap.end())
                continue;
              // Prefer a muon contributor when the digi-hit has several.
              for (const auto& sim : it->second) {
                auto mcp = sim.getParticle();
                if (!mcp.isAvailable())
                  continue;
                if (!hasMCpart) {
                  bestMCpart = mcp;
                  hasMCpart = true;
                }
                if (std::abs(mcp.getPDG()) == 13) {
                  bestMCpart = mcp;
                  break;
                }
              }
            }
            if (hasMCpart) {
              const auto& vtx = bestMCpart.getVertex();                    // mm (EDM4hep)
              double rxy_truth = std::sqrt(vtx.x * vtx.x + vtx.y * vtx.y); // mm
              double d0_reco_mm = 0.0;
              for (int si = 0; si < static_cast<int>(finalTrack.trackStates_size()); ++si) {
                auto ts = finalTrack.getTrackStates(si);
                if (ts.location == edm4hep::TrackState::AtIP) {
                  d0_reco_mm = ts.D0;
                  break;
                }
              }
              bool suspicious = (std::abs(d0_reco_mm) > rxy_truth + 20.0);
              debug() << "  [TruthRxy] rxy_vtx=" << std::fixed << std::setprecision(2) << rxy_truth
                      << " mm  d0_reco=" << d0_reco_mm << " mm"
                      << "  PDG=" << bestMCpart.getPDG() << (suspicious ? "  ** |d0|>rxy+20mm **" : "") << endmsg;
            }
          } catch (const std::exception& e) {
            debug() << "[TruthRxy] Exception: " << e.what() << endmsg;
          }

          // printout: all states vs the matched MC particle
          try {
            if (m_mcParticles.exist()) {
              const auto* mcParts = m_mcParticles.get();
              if (mcParts && !mcParts->empty()) {
                edm4hep::Track tv = finalTrack;
                compareStatesWithTruth(tv, *mcParts, trackNumber);
              }
            }
          } catch (const std::exception& e) {
            debug() << "[TruthCmp] Exception: " << e.what() << endmsg;
          }
        }
      } else {
        debug() << "  [Analytical] AtIP propagation failed for track " << trackNumber << endmsg;
      }
    }

    debug() << "Successfully created track " << trackNumber << " with " << finalTrack.trackerHits_size() << " hits"
            << " and " << finalTrack.trackStates_size() << " track states" << endmsg;

    // ---- per-track run statistics ----
    m_statTracksReconstructed++;
    bool hasLastHit = false, hasIP = false;
    for (int sIdx = 0; sIdx < static_cast<int>(finalTrack.trackStates_size()); ++sIdx) {
      auto ts = finalTrack.getTrackStates(sIdx);
      switch (ts.location) {
      case edm4hep::TrackState::AtLastHit:
        hasLastHit = true;
        break;
      case edm4hep::TrackState::AtIP:
        hasIP = true;
        break;
      default:
        break;
      }
    }
    if (hasLastHit)
      m_statStateAtLastHit++;
    if (hasIP)
      m_statStateAtIP++;

    {
      int nFinalHits = static_cast<int>(trackHits.size());
      if (nFinalHits <= 3)
        m_statNhit3++;
      else if (nFinalHits == 4)
        m_statNhit4++;
      else if (nFinalHits == 5)
        m_statNhit5++;
      else if (nFinalHits == 6)
        m_statNhit6++;
      else
        m_statNhit7plus++;
    }
    m_statValidTriplets++;

    // DEBUG: AtFirstHit state vs the MC particle with most hits on the track
    if (msgLevel(MSG::DEBUG)) {
      std::unordered_map<int, std::pair<edm4hep::MCParticle, int>> mcCount;
      for (const auto& th : trackHits) {
        auto it = recoToSimMap.find(hitKey(th));
        if (it == recoToSimMap.end())
          continue;
        // each MC particle counts once per hit
        std::unordered_set<int> seenInThisHit;
        for (const auto& sim : it->second) {
          auto mcPart = sim.getParticle();
          if (!mcPart.isAvailable())
            continue;
          int mcIdx = mcPart.id().index;
          if (!seenInThisHit.insert(mcIdx).second)
            continue;
          mcCount[mcIdx].first = mcPart;
          mcCount[mcIdx].second++;
        }
      }
      int bestCount = 0;
      edm4hep::MCParticle bestMC;
      bool hasMC = false;
      for (const auto& [idx, pair] : mcCount) {
        if (pair.second > bestCount) {
          bestCount = pair.second;
          bestMC = pair.first;
          hasMC = true;
        }
      }
      if (hasMC) {
        const auto& mom = bestMC.getMomentum();
        double mcPt = std::sqrt(mom.x * mom.x + mom.y * mom.y);
        double mcP = std::sqrt(mcPt * mcPt + mom.z * mom.z);
        double mcEta = (mcPt > 0.0) ? std::asinh(mom.z / mcPt) : 0.0;
        double mcPhi = std::atan2(mom.y, mom.x);

        double recoPt = 0.0, recoPhi = 0.0, recoD0mm = 0.0, recoZ0mm = 0.0, recoEta = 0.0;
        for (int sIdx = 0; sIdx < static_cast<int>(finalTrack.trackStates_size()); ++sIdx) {
          auto ts = finalTrack.getTrackStates(sIdx);
          if (ts.location == edm4hep::TrackState::AtFirstHit) {
            double omegaMm = ts.omega; // 1/mm
            if (std::abs(omegaMm) > 1e-12) {
              recoPt = getPT(ts); // field at the first hit
            }
            recoPhi = ts.phi;
            recoD0mm = ts.D0;
            recoZ0mm = ts.Z0;
            recoEta = (std::abs(ts.tanLambda) > 0.0) ? std::asinh(ts.tanLambda) : 0.0;
            break;
          }
        }

        double dPt = (mcPt > 0) ? (recoPt - mcPt) / mcPt * 100.0 : 0.0;
        double dPhi = recoPhi - mcPhi;
        dPhi = k4rt::wrapAngle(dPhi);

        debug() << "[RECO vs TRUTH trk " << trackNumber << "]"
                << "  PDG=" << bestMC.getPDG() << "  q_truth=" << bestMC.getCharge() << "  hits_matched=" << bestCount
                << "/" << trackHits.size() << endmsg;
        debug() << "  Truth: pT=" << std::fixed << std::setprecision(3) << mcPt << " GeV/c  phi=" << mcPhi
                << " rad  eta=" << mcEta << "  |p|=" << mcP << " GeV/c" << endmsg;
        debug() << "  Reco : pT=" << recoPt << " GeV/c  phi=" << recoPhi << " rad  eta=" << recoEta
                << "  d0=" << recoD0mm << " mm  z0=" << recoZ0mm << " mm" << endmsg;
        debug() << "  Delta: dpT/pT=" << std::setprecision(1) << dPt << "%"
                << "  dphi=" << std::setprecision(4) << dPhi << " rad"
                << "  deta=" << recoEta - mcEta << endmsg;
      } else {
        debug() << "[RECO vs TRUTH trk " << trackNumber << "] no MC link found for track hits" << endmsg;
      }
    }

    // run statistics: sign of omega of the AtFirstHit state
    {
      double omegaLastHit = 0.0;
      bool hasLastHitOmega = false;
      for (int sIdx = 0; sIdx < static_cast<int>(finalTrack.trackStates_size()); ++sIdx) {
        auto ts = finalTrack.getTrackStates(sIdx);
        if (ts.location == edm4hep::TrackState::AtFirstHit) {
          omegaLastHit = ts.omega;
          hasLastHitOmega = true;
          break;
        }
      }
      if (hasLastHitOmega) {
        // omega > 0: clockwise (positive charge in a field along +z)
        if (omegaLastHit > 0.0)
          m_statPositiveCharge++;
        else
          m_statNegativeCharge++;
        debug() << "Charge classification: omega=" << omegaLastHit << " 1/mm (AtFirstHit)  ->  "
                << (omegaLastHit > 0 ? "positive" : "negative") << endmsg;
      }
    }
  } // end of the track loop

  debug() << "Track reconstruction complete: created " << finalTracks.size() << " final tracks across "
          << (trackNumber - 1) << " iterations" << endmsg;
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
  debug() << "  [Anal.Entry] " << (isEndcap ? "endcap" : "barrel") << "  R_o=" << std::fixed << std::setprecision(2)
          << R_outer << " d_oc=" << d_oc << " outerPOCA=" << (d_oc - R_outer)
          << " cm" // outer circle closest approach to IP
          << "  pT=" << std::setprecision(3) << pT << " GeV"
          << "  z_ref=" << std::setprecision(2) << z_ref << " cm" << endmsg;

  // Backward direction on the outer circle (+1: increasing angle, for a clockwise track);
  // the same sign gives the reversed curvature of the inner helix.
  const double inner_sign = clockwise ? 1.0 : -1.0;

  // --- Solenoid boundary crossing ---
  double xc_sol, yc_sol, z_sol;

  if (!isEndcap) {
    // --- barrel: outer circle x solenoid cylinder ---
    if (d_oc < 1e-6) {
      debug() << "  [Analytical] Barrel: outer circle center at IP - degenerate" << endmsg;
      m_statAnalBarrelNoIsect++;
      return false;
    }
    if (R_outer <= std::abs(d_oc - R_sol) || R_outer >= d_oc + R_sol) {
      debug() << "  [Analytical] Barrel: outer circle does not cross solenoid"
              << "  d_oc=" << d_oc << " R_outer=" << R_outer << " R_sol=" << R_sol << endmsg;
      m_statAnalBarrelNoIsect++;
      return false;
    }

    double h = (R_sol * R_sol + d_oc * d_oc - R_outer * R_outer) / (2.0 * d_oc);
    double l = std::sqrt(std::max(0.0, R_sol * R_sol - h * h));
    double ex = cx_o / d_oc, ey = cy_o / d_oc;
    double px = -ey, py = ex;

    double P1x = h * ex + l * px, P1y = h * ey + l * py;
    double P2x = h * ex - l * px, P2y = h * ey - l * py;

    // the crossing reached first going backward from the hit
    double theta_hit = std::atan2(y_ref - cy_o, x_ref - cx_o);
    auto bwdArc = [&](double Px, double Py) -> double {
      double tc = std::atan2(Py - cy_o, Px - cx_o);
      return clockwise ? k4rt::wrapAngleZeroToTwoPi(tc - theta_hit) : k4rt::wrapAngleZeroToTwoPi(theta_hit - tc);
    };
    double arc1 = bwdArc(P1x, P1y);
    double arc2 = bwdArc(P2x, P2y);
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
        debug() << "  [Analytical] Barrel: both crossings have |z_sol|>|z_ref|"
                << " z_ref=" << z_ref << " z_sol=" << z_sol << " z_sol_alt=" << z_sol_alt << " - skip" << endmsg;
        m_statAnalBarrelBothZfail++;
        return false;
      }
      debug() << "  [Analytical] Barrel: swapped to alt crossing"
              << " (|z_sol|=" << std::abs(z_sol) << " > |z_ref|=" << std::abs(z_ref) << ") -> z_sol_alt=" << z_sol_alt
              << " cm" << endmsg;
      m_statAnalBarrelZswapped++;
      xc_sol = xc_alt;
      yc_sol = yc_alt;
      s_bwd = s_bwd_alt;
      z_sol = z_sol_alt;
    }

    debug() << "  [Anal.Brl] crossing=(" << std::fixed << std::setprecision(2) << xc_sol << ", " << yc_sol
            << ") cm  bwd_arc=" << s_bwd << " cm  z_sol=" << z_sol << " cm" << endmsg;

  } else {
    // --- end disk: arc s backward from the hit with z(s) = +-Z_sol ---
    if (std::abs(tanLambda) < 1e-6) {
      debug() << "  [Analytical] Endcap: |tanLambda| too small to extrapolate z" << endmsg;
      m_statAnalEndcapSmallTanL++;
      return false;
    }

    double z_cap = (z_ref < 0.0) ? -Z_sol : +Z_sol;

    // z_ref - tanLambda * s = z_cap
    double s_endcap = (z_ref - z_cap) / tanLambda;
    if (s_endcap < 0.0) {
      // hit already at |z| < Z_sol: use the barrel crossing instead
      debug() << "  [Analytical] Endcap: z_ref=" << z_ref << " cm already inside z_cap=" << z_cap
              << " cm - falling back to barrel formula" << endmsg;
      m_statAnalEndcapPastCap++;
      if (d_oc < 1e-6 || R_outer <= std::abs(d_oc - R_sol) || R_outer >= d_oc + R_sol) {
        debug() << "  [Analytical] Endcap->Barrel fallback: no barrel crossing" << endmsg;
        return false;
      }
      double h = (R_sol * R_sol + d_oc * d_oc - R_outer * R_outer) / (2.0 * d_oc);
      double l = std::sqrt(std::max(0.0, R_sol * R_sol - h * h));
      double ex = cx_o / d_oc, ey = cy_o / d_oc;
      double px = -ey, py = ex;
      double P1x = h * ex + l * px, P1y = h * ey + l * py;
      double P2x = h * ex - l * px, P2y = h * ey - l * py;
      double theta_hit = std::atan2(y_ref - cy_o, x_ref - cx_o);
      auto bwdArc2 = [&](double Px, double Py) -> double {
        double tc = std::atan2(Py - cy_o, Px - cx_o);
        return clockwise ? k4rt::wrapAngleZeroToTwoPi(tc - theta_hit) : k4rt::wrapAngleZeroToTwoPi(theta_hit - tc);
      };
      double arc1 = bwdArc2(P1x, P1y), arc2 = bwdArc2(P2x, P2y);
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
        debug() << "  [Analytical] Endcap->Barrel fallback: z_sol=" << z_sol << " cm out of range" << endmsg;
        return false;
      }
      debug() << "  [Anal.Brl*] fallback crossing=(" << std::fixed << std::setprecision(2) << xc_sol << ", " << yc_sol
              << ") cm  bwd_arc=" << s_bwd_fb << " cm  z_sol=" << z_sol << " cm" << endmsg;
    } else {
      if (s_endcap > 2.0 * M_PI * R_outer) {
        debug() << "  [Analytical] Endcap: backward arc would exceed one loop - skip" << endmsg;
        m_statAnalEndcapArcLong++;
        return false;
      }

      double theta_ref = std::atan2(y_ref - cy_o, x_ref - cx_o);
      double d_theta = s_endcap / R_outer;
      double theta_cap = theta_ref + inner_sign * d_theta;

      xc_sol = cx_o + R_outer * std::cos(theta_cap);
      yc_sol = cy_o + R_outer * std::sin(theta_cap);
      z_sol = z_cap;

      debug() << "  [Anal.End] crossing=(" << std::fixed << std::setprecision(2) << xc_sol << ", " << yc_sol << ", "
              << z_sol << ") cm  bwd_arc=" << s_endcap << " cm" << endmsg;
    } // end if (s_endcap < 0.0)
  } // end if (!isEndcap)

  // --- tangent and outward radial direction at the crossing ---
  double rx = (xc_sol - cx_o) / R_outer;
  double ry = (yc_sol - cy_o) / R_outer;
  double tx = clockwise ? -ry : ry;
  double ty = clockwise ? rx : -rx;
  debug() << "  [Analytical] Tangent: (" << tx << ", " << ty << ")  outward radial: (" << rx << ", " << ry << ")"
          << endmsg;

  // --- inner field, at the origin ---
  double B_inner = m_field.magneticField(dd4hep::Position(0.0, 0.0, 0.0)).z() / dd4hep::tesla;
  if (std::abs(B_inner) < 1e-6) {
    debug() << "  [Analytical] Inner B field vanishes (" << B_inner << " T)" << endmsg;
    return false;
  }

  // --- inner circle: same pT, field reversed, so the centre is on the other side ---
  double R_inner = k4rt::radiusFromPt(B_inner, pT); // cm
  double cx_i = xc_sol + R_inner * rx;
  double cy_i = yc_sol + R_inner * ry;

  debug() << "  [Anal.Inner] R_i=" << std::fixed << std::setprecision(2) << R_inner << " cx_i=" << cx_i
          << " cy_i=" << cy_i << " B_i=" << std::setprecision(3) << B_inner << " T"
          << "  d0=" << std::setprecision(2) << (std::sqrt(cx_i * cx_i + cy_i * cy_i) - R_inner) << " cm" << endmsg;

  // --- inner helix at its perigee ---
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

  debug() << "  [Analytical] Inner perigee:"
          << "  d0=" << d0_inner_cm << " cm"
          << "  phi=" << phi0_inner << " rad"
          << "  omega=" << omega_inner_mm << " /mm"
          << "  z0=" << z0_inner_cm << " cm"
          << "  (bwd arc=" << s_sol_to_ip << " cm)" << endmsg;

  // --- stored as AtIP ---
  edm4hep::TrackState atIPState = createTrackState(d0_inner_cm * 10.0, // cm -> mm
                                                   phi0_inner, omega_inner_mm,
                                                   z0_inner_cm * 10.0, // cm -> mm
                                                   tanLambda, edm4hep::TrackState::AtIP);
  finalTrack.addToTrackStates(atIPState);

  debug() << "  [Analytical] AtIP (" << (isEndcap ? "endcap" : "barrel") << ") saved:"
          << "  d0=" << d0_inner_cm << " cm"
          << "  phi=" << phi0_inner << " rad"
          << "  omega=" << omega_inner_mm << " /mm"
          << "  z0=" << z0_inner_cm << " cm" << endmsg;

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

  // --- Step 1: Algebraic WLS circle fit ---
  // Circle equation: x^2 + y^2 + D*x + E*y + F = 0
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
    debug() << "fitCircleNHits: degenerate algebraic fit (disc=" << disc << ")" << endmsg;
    return false;
  }
  radius = std::sqrt(disc) / 2.0;

  // --- Step 2: Gauss-Newton geometric refinement ---
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
    debug() << "fitCircleNHits: non-positive radius after refinement" << endmsg;
    return false;
  }

  // --- Step 3: Covariance on [x0, y0, R] from final Jacobian ---
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

  // --- Step 4: Per-hit residuals and chi2/NDF ---
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

  debug() << "fitCircleNHits (" << N << " hits): center=(" << x0 << "," << y0 << ") cm  R=" << radius
          << " cm  chi2/ndf=" << chi2ndf << endmsg;
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
  const double scanLimitMm = 10000.0;                         // 10 m: beyond any FCC-ee solenoid
  const double radiusForZScanMm = 0.5 * m_solenoidBoundaryMm; // the radial scan runs first

  // Bz in tesla at a distance (mm) along the scan direction
  auto getBz = [&](double distanceMm) -> double {
    const double distanceCm = distanceMm / 10.0; // DD4hep positions are in cm
    dd4hep::Position position(distanceCm, 0.0, 0.0);
    if (alongZ)
      position = dd4hep::Position(radiusForZScanMm / 10.0, 0.0, distanceCm);
    return m_field.magneticField(position).z() / dd4hep::tesla;
  };

  const double bzCentral = getBz(startMm);
  debug() << "  " << tag << " scan: B_z = " << bzCentral << " T near the IP" << endmsg;
  if (std::abs(bzCentral) < 1e-6) {
    warning() << "  " << tag << ": no field near the IP; the solenoid boundary is set to " << scanLimitMm << " mm"
              << endmsg;
    return scanLimitMm;
  }
  // true when the field at this distance is no longer the central solenoid field
  auto outsideCentralField = [&](double distanceMm) -> bool {
    return std::abs(getBz(distanceMm) - bzCentral) > 0.5 * std::abs(bzCentral);
  };

  // 1. Coarse scan: the 1 mm interval [lowMm, highMm] where the central field ends
  double lowMm = -1.0;
  double highMm = -1.0;
  for (double distanceMm = startMm + stepMm; distanceMm <= scanLimitMm; distanceMm += stepMm) {
    if (outsideCentralField(distanceMm)) {
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
  debug() << "  " << tag << " coarse interval: " << lowMm << " mm (B_z = " << getBz(lowMm) << " T) to " << highMm
          << " mm (B_z = " << getBz(highMm) << " T)" << endmsg;

  // 2. Binary search inside the interval, down to 0.01 mm
  for (int iteration = 0; iteration < 40; ++iteration) {
    const double middleMm = 0.5 * (lowMm + highMm);
    if (outsideCentralField(middleMm))
      highMm = middleMm;
    else
      lowMm = middleMm;
    if (highMm - lowMm < 0.01)
      break;
  }

  const double boundaryMm = 0.5 * (lowMm + highMm);
  debug() << "  " << tag << " boundary = " << std::fixed << std::setprecision(2) << boundaryMm << " mm" << endmsg;
  return boundaryMm;
}

// ============================================================================
// TRUTH COMPARISON
// ============================================================================
void CombinatorialTrackFinder::compareStatesWithTruth(const edm4hep::Track& track,
                                                      const edm4hep::MCParticleCollection& mcParticles,
                                                      int trackNumber) const {

  // MC particle with the highest pT (single-particle samples)
  const edm4hep::MCParticle* mu = nullptr;
  double bestPT2 = -1.0;
  for (const auto& p : mcParticles) {
    const auto& m = p.getMomentum();
    double pT2 = m.x * m.x + m.y * m.y;
    if (pT2 > bestPT2) {
      bestPT2 = pT2;
      mu = &p;
    }
  }
  if (!mu) {
    debug() << "[TruthCmp] No MCParticle for track " << trackNumber << endmsg;
    return;
  }

  const auto& vtx = mu->getVertex();   // mm (EDM4hep convention)
  const auto& mom = mu->getMomentum(); // GeV/c
  double pT_true = std::sqrt(mom.x * mom.x + mom.y * mom.y);
  double phi0_true = std::atan2(mom.y, mom.x);
  double tanL_true = (pT_true > 1e-9) ? mom.z / pT_true : 0.0;
  // straight-line d0 and z0 of the production vertex
  double d0_true_mm = -vtx.x * std::sin(phi0_true) + vtx.y * std::cos(phi0_true);
  double z0_true_mm = vtx.z;

  // fields inside (origin) and outside (50 cm beyond the boundary) the solenoid
  double B_ip = m_field.magneticField(dd4hep::Position(0.0, 0.0, 0.0)).z() / dd4hep::tesla;
  double R_outer_cm = (m_solenoidBoundaryMm / 10.0) + 50.0;
  double B_outer = m_field.magneticField(dd4hep::Position(R_outer_cm, 0.0, 0.0)).z() / dd4hep::tesla;

  // --- Print header ---
  debug() << "\n"
          << "+==========================================================================+\n"
          << "|  TRUTH COMPARISON  Track " << std::setw(3) << trackNumber << "   PDG=" << std::setw(5) << mu->getPDG()
          << "   pT=" << std::fixed << std::setprecision(2) << std::setw(7) << pT_true << " GeV/c"
          << "   tanLambda=" << std::setw(6) << std::setprecision(3) << tanL_true << "            |\n"
          << "|  MC vertex: (" << std::setw(7) << std::setprecision(2) << vtx.x << ", " << std::setw(7) << vtx.y << ", "
          << std::setw(7) << vtx.z << ") mm"
          << "   B_ip=" << std::setw(6) << std::setprecision(3) << B_ip << " T   B_outer=" << std::setw(6) << B_outer
          << " T           |\n"
          << "+=================+=============+=============+============+=============+\n"
          << "|  State          |  D0 (mm)    |  Z0 (mm)    |  phi (rad) |  pT (GeV/c) |\n"
          << "+=================+=============+=============+============+=============+\n";

  // Truth row
  debug() << "|  Truth (MC IP)  |" << std::setw(11) << std::setprecision(3) << d0_true_mm << "  |" << std::setw(11)
          << z0_true_mm << "  |" << std::setw(10) << std::setprecision(4) << phi0_true << "  |" << std::setw(11)
          << std::setprecision(2) << pT_true << "  |" << endmsg;
  debug() << "+=================+=============+=============+============+=============+\n";

  // Helper: print one state row + residual row
  struct StateRow {
    const char* label;
    int loc;
    double B_for_pT;
  };
  StateRow rows[] = {
      {"AtFirstHit", edm4hep::TrackState::AtFirstHit, B_outer},
      {"AtVertex", edm4hep::TrackState::AtVertex, B_ip},
      {"AtIP(analyt)", edm4hep::TrackState::AtIP, B_ip},
  };
  for (auto& row : rows) {
    bool found = false;
    for (int i = 0; i < static_cast<int>(track.trackStates_size()); ++i) {
      auto ts = track.getTrackStates(i);
      if (ts.location != row.loc)
        continue;
      found = true;
      double pT_reco = (std::abs(ts.omega) > 1e-12 && std::abs(row.B_for_pT) > 1e-6)
                           ? k4rt::ptFromOmega(row.B_for_pT, ts.omega)
                           : 0.0;
      double resD0 = ts.D0 - d0_true_mm;
      double resZ0 = ts.Z0 - z0_true_mm;
      double resPhi = ts.phi - phi0_true;
      double resPT = pT_reco - pT_true;
      // phi residual: wrap to (-pi, pi]
      resPhi = k4rt::wrapAngle(resPhi);
      debug() << "|  " << std::left << std::setw(15) << row.label << "|" << std::right << std::setw(11)
              << std::setprecision(3) << ts.D0 << "  |" << std::setw(11) << ts.Z0 << "  |" << std::setw(10)
              << std::setprecision(4) << ts.phi << "  |" << std::setw(11) << std::setprecision(2) << pT_reco << "  |"
              << endmsg;
      debug() << "|  " << std::left << std::setw(15) << "  d (reco-truth)" << "|" << std::right << std::setw(9)
              << std::setprecision(3) << (resD0 >= 0 ? "+" : "") << resD0 << "   |" << std::setw(9)
              << (resZ0 >= 0 ? "+" : "") << resZ0 << "   |" << std::setw(8) << std::setprecision(4)
              << (resPhi >= 0 ? "+" : "") << resPhi << "    |" << std::setw(9) << std::setprecision(2)
              << (resPT >= 0 ? "+" : "") << resPT << "   |" << endmsg;
      break;
    }
    if (!found)
      debug() << "|  " << std::left << std::setw(15) << row.label
              << "|  (not saved)                                              |" << endmsg;
  }
  debug() << "+=================+=============+=============+============+=============+\n" << endmsg;
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
