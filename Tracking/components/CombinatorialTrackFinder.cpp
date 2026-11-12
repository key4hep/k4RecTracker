/**
 * CombinatorialTrackFinder: track finding on the planar hits of one tracking system
 * (overview in CombinatorialTrackFinder.h).
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

#include "CombinatorialTrackFinder.h"
#include "ExhaustiveSearch.h"
#include "TrackCandidate.h"
#include "TrackingMath.h"

#include "DD4hep/DD4hepUnits.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <iomanip>
#include <limits>

DECLARE_COMPONENT(CombinatorialTrackFinder)

using k4rt::CircleFit;
using k4rt::kCmPerMm;
using k4rt::kMmPerCm;
using k4rt::TrackCandidate;

namespace {

// ---- fixed constants ----
constexpr int kMaterialSamplesPerLayer = 24;        // layer material: sampled on up to 24 sensors,
constexpr double kLayerMaterialHalfDepthCm = 150.0; // each side of a layer at most 150 cm deep
constexpr double kPairTolerancePhi = 0.15;          // two hits of a doublet are coincident within
constexpr double kPairToleranceEta = 0.10;          // these |dphi| (rad) and |deta|
constexpr int kMinCompleteDoublets = 2;             // pair coincidence: complete doublets needed
constexpr double kGuidedGeometricWeight = 0.75;     // crowded layers: score = 0.75 * geometric term
constexpr double kGuidedEdepWeight = 0.25;          //                       + 0.25 * edep term
constexpr double kCrossingZToleranceCm = 1e-3;      // SolenoidCrossing: crossing |z| <= |z_hit| + 10 um

// key of a hit: (collectionID, index); merged input collections all start their indices at 0
template <typename Hit>
uint64_t hitKey(const Hit& hit) {
  const auto id = hit.id();
  return (static_cast<uint64_t>(id.collectionID) << 32) | static_cast<uint32_t>(id.index);
}

enum class StatePoint { FirstHit, LastHit, Origin }; ///< where a stored track state is placed

// formula helpers, defined at the end of the file
bool makeEdm4hepStateFromCircle(const CircleFit& circle, const std::vector<edm4hep::TrackerHitPlane>& trackHits,
                                StatePoint where, edm4hep::TrackState& state);
edm4hep::TrackState stateWithDefaultCovariance(double d0, double phi, double omega, double z0, double tanLambda,
                                               int location);
void circleCylinderCrossings(double centreX, double centreY, double radius, double centreDistance,
                             double solenoidRadius, double crossingX[2], double crossingY[2]);
double backwardArcAngle(double centreX, double centreY, double angleOfHit, bool clockwise, double pointX,
                        double pointY);
bool nextLayerChoice(std::vector<int>& chosen, int numberOfLayers);
bool nextHitChoice(std::vector<int>& position, const std::vector<const std::vector<size_t>*>& hitsOfLayer);

} // namespace

// ============================================================================================
// set-up
// ============================================================================================

CombinatorialTrackFinder::CombinatorialTrackFinder(const std::string& name, ISvcLocator* pSvcLocator)
    : MultiTransformer(
          name, pSvcLocator,
          {KeyValues{"InputHitCollection", {"TrackerHits"}}, KeyValues{"HeaderCollectionName", {"EventHeader"}},
           KeyValues{"InputRecoSimLinkCollection", {"TrackerHitSimTrackerHitLink"}}},
          {KeyValues{"OutputTrackCollection", {"DisplacedTracks"}}, KeyValues{"OutputHitCollection", {"TrackHits"}},
           KeyValues{"OutputRecoSimLinkCollection", {"DisplacedTrackHitSimLinks"}}}) {}

StatusCode CombinatorialTrackFinder::initialize() {
  StatusCode status = Algorithm::initialize();
  if (!status.isSuccess())
    return status;
  if (!m_geoSvc.retrieve()) {
    error() << "Failed to retrieve GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  m_detector = m_geoSvc->getDetector();
  if (!m_detector) {
    error() << "Failed to access detector" << endmsg;
    return StatusCode::FAILURE;
  }
  if (setupCellIDDecoding().isFailure() || setupSurfaces().isFailure() || readSearchOptions().isFailure())
    return StatusCode::FAILURE;
  if (m_strategyChoice == Strategy::SeedAndExtend && setupSeedFinder().isFailure())
    return StatusCode::FAILURE;
  info() << "Found " << m_surfaces.size() << " detector surfaces in " << m_surfacesByLayer.size() << " layers"
         << endmsg;
  findSolenoid();
  return StatusCode::SUCCESS;
}

// the cellID decoder and the readout fields of the layer key
StatusCode CombinatorialTrackFinder::setupCellIDDecoding() {
  const std::string encoding = m_geoSvc->constantAsString(m_encodingStringParameter);
  m_bitFieldCoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(encoding);
  const auto& description = m_bitFieldCoder->fieldDescription();
  if (!readoutHasField(m_layerFieldName.value())) {
    error() << "Readout '" << m_encodingStringParameter.value() << "' (" << description << ") has no field '"
            << m_layerFieldName.value() << "'. Set LayerFieldName to match this readout." << endmsg;
    return StatusCode::FAILURE;
  }
  m_hasRegionField = readoutHasField(m_regionFieldName.value());
  if (!m_hasRegionField) {
    warning() << "Readout '" << m_encodingStringParameter.value() << "' (" << description << ") has no field '"
              << m_regionFieldName.value() << "'; treating every hit as barrel. Set RegionFieldName "
              << "(\"type\" for MuonSystemReadoutID, \"side\" for GlobalTrackerReadoutID)." << endmsg;
  }
  // several subdetectors in one readout: "system" goes into the layer key, so equally numbered
  // layers of different subdetectors stay apart
  const bool severalDetectors = m_detectorNames.value().size() > 1;
  m_useSystemInLayerKey = severalDetectors && readoutHasField("system");
  if (severalDetectors && !m_useSystemInLayerKey) {
    warning() << m_detectorNames.value().size() << " detectors share this readout but "
              << "'system' is NOT folded into the layer key -- layers with the same "
              << "index in different subdetectors will be merged." << endmsg;
  }
  info() << "CellID decoding: layer field '" << m_layerFieldName.value() << "', region field '"
         << m_regionFieldName.value() << (m_hasRegionField ? "'" : "' (absent, barrel-only)")
         << ", system in layer key: " << (m_useSystemInLayerKey ? "yes" : "no") << endmsg;
  return StatusCode::SUCCESS;
}

// surface map (one detector or merged), field, surfaces by layer, region of every layer key, and
// the layer geometry
StatusCode CombinatorialTrackFinder::setupSurfaces() {
  auto* surfaceManager = m_detector->extension<dd4hep::rec::SurfaceManager>();
  const std::vector<std::string> detectorNames = m_detectorNames.value();
  if (detectorNames.empty()) {
    error() << "DetectorNames is empty: give the name(s) of the detector(s) to track in" << endmsg;
    return StatusCode::FAILURE;
  }
  if (detectorNames.size() == 1) {
    m_surfaceMap = surfaceManager->map(detectorNames.front());
    if (!m_surfaceMap) {
      error() << "Could not find surface map for detector: " << detectorNames.front() << endmsg;
      return StatusCode::FAILURE;
    }
  } else {
    m_mergedSurfaceMap.clear();
    for (const auto& detectorName : detectorNames) {
      const auto* detectorMap = surfaceManager->map(detectorName);
      if (!detectorMap) {
        error() << "Could not find surface map for detector: " << detectorName << endmsg;
        return StatusCode::FAILURE;
      }
      m_mergedSurfaceMap.insert(detectorMap->begin(), detectorMap->end());
      info() << "  merged " << detectorMap->size() << " surfaces from '" << detectorName << "'" << endmsg;
    }
    m_surfaceMap = &m_mergedSurfaceMap;
    info() << "Merged surface map over " << detectorNames.size() << " detectors: " << m_surfaceMap->size()
           << " surfaces" << endmsg;
  }

  m_field = m_detector->field();
  const dd4hep::Direction fieldAtCentre = m_field.magneticField(dd4hep::Position(0, 0, 0));
  info() << "Magnetic field at the detector centre: (" << fieldAtCentre.x() / dd4hep::tesla << ", "
         << fieldAtCentre.y() / dd4hep::tesla << ", " << fieldAtCentre.z() / dd4hep::tesla << ") T" << endmsg;

  m_commonConfig.minTrackHits = m_minTrackHits.value();
  m_commonConfig.minPtGeV = m_minPt.value();
  m_commonConfig.maxPairHitDistanceCm = m_maxPairHitDistance.value();
  m_commonConfig.bFieldTesla = fieldAtCentre.z() / dd4hep::tesla;
  m_commonConfig.sigmaHitDefaultCm = m_sigmaHitDefault.value();

  m_surfaces.clear();
  for (const auto& [cellID, surface] : *m_surfaceMap) {
    const auto* concreteSurface = dynamic_cast<const dd4hep::rec::Surface*>(surface);
    if (concreteSurface)
      m_surfaces.push_back(concreteSurface);
  }
  m_surfacesByLayer = groupSurfacesByLayer();

  // layer key -> region from the readout's region field; two layers must not share a key
  m_regionOfLayerKey.clear();
  int keysWithOtherRegion = 0;
  for (const auto* surface : m_surfaces) {
    const dd4hep::DetElement detElement = surface->detElement();
    if (!detElement.isValid())
      continue;
    const uint64_t volumeID = detElement.volumeID();
    const int layerKey = layerKeyOf(volumeID);
    const int region = regionOf(volumeID);
    const auto [entry, inserted] = m_regionOfLayerKey.emplace(layerKey, region);
    if (!inserted && entry->second != region) {
      error() << "layer key " << layerKey << " is produced by surfaces in different regions (" << entry->second
              << " and " << region << "): the 100*system+layer encoding collides for this detector" << endmsg;
      return StatusCode::FAILURE;
    }
    if (inserted && region != k4rt::regionOfLayer(layerKey, nullptr))
      ++keysWithOtherRegion;
  }
  if (keysWithOtherRegion > 0)
    info() << keysWithOtherRegion << " layer key(s) whose region differs from the key-range rule; "
           << "using the readout region field" << endmsg;

  m_layerGeometry.build(m_surfacesByLayer);
  info() << "Layer geometry (inside-out by closest approach to the IP):" << endmsg;
  for (const auto& layer : m_layerGeometry.layers()) {
    info() << "  layer key " << std::setw(6) << layer.layerKey << "  " << std::setw(7) << k4rt::kindName(layer.kind)
           << "  r=[" << std::fixed << std::setprecision(2) << layer.rMin_cm << "," << layer.rMax_cm << "]"
           << "  z=[" << layer.zMin_cm << "," << layer.zMax_cm << "] cm"
           << "  dMin=" << layer.closestApproachCm << "  nSurf=" << layer.numberOfSurfaces << std::defaultfloat
           << std::setprecision(6) << endmsg;
  }
  for (const auto& problem : m_layerGeometry.validate())
    warning() << "LayerGeometry: " << problem << endmsg;
  return StatusCode::SUCCESS;
}

// text options -> enums; the doublets of the pair coincidence
StatusCode CombinatorialTrackFinder::readSearchOptions() {
  const std::string& ipStateMode = m_ipStateMode.value();
  if (ipStateMode == "None") {
    m_ipState = IPState::None;
  } else if (ipStateMode == "Helix") {
    m_ipState = IPState::Helix;
  } else if (ipStateMode == "SolenoidCrossing") {
    m_ipState = IPState::SolenoidCrossing;
  } else {
    error() << "Unknown IPStateMode '" << ipStateMode << "' (expected None, Helix or SolenoidCrossing)" << endmsg;
    return StatusCode::FAILURE;
  }

  if (m_strategy.value() == "Exhaustive") {
    m_strategyChoice = Strategy::Exhaustive;
  } else if (m_strategy.value() == "SeedAndExtend") {
    m_strategyChoice = Strategy::SeedAndExtend;
  } else {
    error() << "Unknown Strategy '" << m_strategy.value() << "' (expected Exhaustive or SeedAndExtend)" << endmsg;
    return StatusCode::FAILURE;
  }

  const std::vector<int>& doubletLayers = m_doubletLayers.value(); // two by two
  if (doubletLayers.size() % 2 != 0) {
    error() << "DoubletLayers needs an even number of layer numbers (two per doublet)" << endmsg;
    return StatusCode::FAILURE;
  }
  m_doubletLayerPairs.clear();
  m_partnerLayer.clear();
  for (size_t i = 0; i + 1 < doubletLayers.size(); i += 2) {
    m_doubletLayerPairs.emplace_back(doubletLayers[i], doubletLayers[i + 1]);
    m_partnerLayer[doubletLayers[i]] = doubletLayers[i + 1];
    m_partnerLayer[doubletLayers[i + 1]] = doubletLayers[i];
  }
  return StatusCode::SUCCESS;
}

// SeedAndExtend: the material of every layer and the seed finder
StatusCode CombinatorialTrackFinder::setupSeedFinder() {
  k4rt::SeedConfig seedConfig;
  const std::string& policy = m_seedLayerPolicy.value();
  if (policy == "innermost") {
    seedConfig.policy = k4rt::SeedLayerPolicy::Innermost;
  } else if (policy == "outermost") {
    seedConfig.policy = k4rt::SeedLayerPolicy::Outermost;
  } else if (policy == "consecutive") {
    seedConfig.policy = k4rt::SeedLayerPolicy::Consecutive;
  } else {
    error() << "Unknown SeedLayerPolicy '" << policy << "' (expected innermost, outermost or consecutive)" << endmsg;
    return StatusCode::FAILURE;
  }
  seedConfig.maxCandidatesPerLayer = m_maxCandidatesPerLayer.value();
  seedConfig.shortTrackMaxHits = m_shortTrackMaxHits.value();
  seedConfig.shadowHitWindow_cm = m_shadowHitWindow.value();
  seedConfig.kinkMaxChi2Triplet = m_kinkMaxChi2Triplet.value();
  seedConfig.seedKinkMaxChi2NDF = m_seedKinkMaxChi2NDF.value();
  seedConfig.trackKinkMaxChi2NDF = m_trackKinkMaxChi2NDF.value();
  seedConfig.common = m_commonConfig;

  // layer material, for the multiple scattering in the kink chi2
  m_layerGeometry.measureMaterial(*m_detector, m_surfacesByLayer, kMaterialSamplesPerLayer, kLayerMaterialHalfDepthCm);
  auto& materialLine = info();
  materialLine << "SeedAndExtend: B = " << seedConfig.common.bFieldTesla << " T; layer x/X0 [%]:";
  for (const auto& layer : m_layerGeometry.layers())
    materialLine << " " << layer.layerKey << ":" << std::fixed << std::setprecision(2) << 100.0 * layer.xOverX0;
  materialLine << std::defaultfloat << std::setprecision(6) << endmsg;

  // adapter: the sensor plane of a hit, for the seed finder
  auto sensorPlane = [this](uint64_t cellID, Eigen::Vector3d& origin, Eigen::Vector3d& normal) {
    const dd4hep::rec::Surface* surface = surfaceOf(cellID);
    if (!surface)
      return false;
    const auto surfaceOrigin = surface->origin();
    const auto surfaceNormal = surface->normal(surfaceOrigin);
    origin = Eigen::Vector3d(surfaceOrigin.x(), surfaceOrigin.y(), surfaceOrigin.z());
    normal = Eigen::Vector3d(surfaceNormal.x(), surfaceNormal.y(), surfaceNormal.z()).normalized();
    return true;
  };
  m_seedFinder = std::make_unique<k4rt::SeedFinder>(m_layerGeometry, seedConfig, sensorPlane);
  for (const auto& problem : m_seedFinder->prepare())
    warning() << "SeedFinder: " << problem << endmsg;
  auto& passesLine = info();
  passesLine << "Strategy SeedAndExtend: passes";
  for (int seedSize : m_seedFinder->passes())
    passesLine << " [" << seedSize << " hits: " << m_seedFinder->numSets(seedSize) << " layer sets]";
  passesLine << " (policy " << policy << ")" << endmsg;
  if (m_seedFinder->numSets(m_seedFinder->passes().front()) == 0) {
    error() << "SeedAndExtend: no seed layer sets - check SeedLayerPolicy" << endmsg;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

// the solenoid boundary (radius first: the z scan runs at half of it)
void CombinatorialTrackFinder::findSolenoid() {
  m_solenoidRadiusMm = findSolenoidBoundary(false);
  info() << "Solenoid Bz boundary (radial): R = " << std::fixed << std::setprecision(4) << m_solenoidRadiusMm
         << " mm  (" << m_solenoidRadiusMm / kMmPerCm << " cm)" << endmsg;
  m_solenoidHalfLengthMm = findSolenoidBoundary(true);
  info() << "Solenoid Bz boundary (endcap): |Z| = " << std::fixed << std::setprecision(4) << m_solenoidHalfLengthMm
         << " mm  (" << m_solenoidHalfLengthMm / kMmPerCm << " cm)" << endmsg;
}

// ============================================================================================
// the search in one event
// ============================================================================================

namespace {

/// A group of hits searched on its own: its layers (inner to outer) and the hit indices of each.
struct HitCluster {
  std::vector<int> layerKeys;
  std::vector<std::vector<size_t>> hitsPerLayer;
};

/// Union-find for the hit clustering: unite() joins two groups, find() gives a group's representative.
struct UnionFind {
  std::vector<size_t> parent;
  explicit UnionFind(size_t size) : parent(size) {
    for (size_t i = 0; i < size; ++i)
      parent[i] = i;
  }
  size_t find(size_t element) {
    while (parent[element] != element) {
      parent[element] = parent[parent[element]]; // path halving
      element = parent[element];
    }
    return element;
  }
  void unite(size_t first, size_t second) {
    first = find(first);
    second = find(second);
    if (first != second)
      parent[first] = second;
  }
};

/// The best combination so far, and the best above MaxComboPT (used only if nothing is below it).
struct SearchResult {
  TrackCandidate best;
  TrackCandidate fallback;
  bool found = false;
};

} // namespace

// One track per iteration, until no candidate is left.
//   Exhaustive: collect the usable hits, find the best combination, refine it (crowded layers,
//   neighbour tracks, outliers), store it.
//   SeedAndExtend: the SeedFinder builds all candidates once; they are stored best first, skipping
//   those sharing a hit with an accepted track.
class CombinatorialTrackFinder::EventSearch {
public:
  EventSearch(const CombinatorialTrackFinder& finder, const edm4hep::TrackerHitPlaneCollection& hits,
              edm4hep::TrackCollection& outputTracks, edm4hep::TrackerHitPlaneCollection& outputHits,
              const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks,
              edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinks);

  void run();

private:
  // ---- the two strategies ----
  void runExhaustive();
  void runSeedAndExtend();

  // ---- the usable hits of one iteration ----
  bool collectPool();
  void prepareExhaustiveIteration();
  void computeEdepScores();
  void findIsolationNeighbours();
  double fieldAtCentroid(const std::vector<size_t>& hitIndices) const;

  // ---- Exhaustive: the best combination ----
  bool searchCombinations(TrackCandidate& best);
  HitCluster pairCoincidenceCluster(const HitCluster& allHits) const;
  std::vector<HitCluster> buildClusters(const HitCluster& allHits) const;
  void joinLayers(UnionFind& groups, const HitCluster& allHits, const std::vector<size_t>& poolPosition, int step,
                  double phiLimit, double distanceLimit) const;
  void searchCluster(const HitCluster& cluster, const k4rt::RoadCuts& roadCuts, SearchResult& result) const;
  bool passesPreFitCuts(const std::vector<size_t>& combination, const k4rt::RoadCuts& roadCuts) const;
  void offerCombination(const std::vector<size_t>& combination, SearchResult& result) const;
  double meanEdepScore(const std::vector<size_t>& combination) const;
  bool isBetterCandidate(int numberOfHits, const CircleFit& fit, double edepScoreKeV, const TrackCandidate& best) const;
  static bool isBetterFallback(int numberOfHits, const CircleFit& fit, const TrackCandidate& fallback);

  // ---- Exhaustive: refinement of the chosen combination ----
  void selectCrowdedLayerHits(TrackCandidate& best) const;
  std::vector<size_t> cleanLayerHits() const;
  bool addBestHitOfCrowdedLayers(const std::function<double(size_t)>& geometricTerm, std::vector<size_t>& guided) const;
  void applyGuidedChoice(const std::vector<size_t>& guided, TrackCandidate& best) const;
  bool rejectedAsNeighbour(const TrackCandidate& best);
  bool nearAcceptedTrack(size_t hitIndex, double distanceSquaredMm2) const;
  void removeOutliers(TrackCandidate& best) const;

  // ---- hit tests of the exhaustive search ----
  bool isolated(const std::vector<size_t>& combination) const;
  bool coincident(size_t first, size_t second) const;
  int completeDoublets(const std::vector<size_t>& combination) const;

  // ---- output ----
  void buildTrack(const TrackCandidate& candidate);
  std::optional<edm4hep::TrackState> addHelixStates(const TrackCandidate& candidate,
                                                    const std::vector<edm4hep::TrackerHitPlane>& trackHits,
                                                    edm4hep::MutableTrack& track) const;
  void addSolenoidCrossingState(const TrackCandidate& candidate, const std::vector<edm4hep::TrackerHitPlane>& trackHits,
                                const std::optional<edm4hep::TrackState>& firstHitState,
                                edm4hep::MutableTrack& track) const;
  void copyHitToOutput(const edm4hep::TrackerHitPlane& hit, edm4hep::MutableTrack& track);

  // ---- the event ----
  const CombinatorialTrackFinder& finder;
  const edm4hep::TrackerHitPlaneCollection& hits;
  edm4hep::TrackCollection& outputTracks;
  edm4hep::TrackerHitPlaneCollection& outputHits;
  edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinks;
  std::unordered_map<uint64_t, std::vector<edm4hep::SimTrackerHit>> recoToSimHits; ///< hitKey -> its SimTrackerHits
  // per input hit (index = hit index)
  std::vector<bool> hasSurface;                                 ///< on a known detector surface
  std::vector<k4rt::HitGeometry> hitGeometry;                   ///< position, azimuth, layer key, region
  std::vector<bool> used;                                       ///< by an accepted track (or given up as a neighbour)
  std::vector<std::array<double, 3>> acceptedTrackHitPositions; ///< mm, for the neighbour-track cuts
  int trackNumber = 0;                                          ///< current iteration (messages)
  bool pairCoincidence = false;                                 ///< RequirePairCoincidence

  // ---- the current iteration ----
  std::vector<size_t> poolHits;                         ///< the unused usable hits, in collection order
  std::map<int, std::vector<size_t>> hitsByLayer;       ///< layer key -> pool hits (after MaxHitsPerCompositeID)
  std::vector<int> layerOrder;                          ///< layer keys of the pool, inner -> outer
  std::vector<double> edepScore;                        ///< keV per hit (0 when not used)
  std::vector<std::vector<size_t>> isolationNeighbours; ///< pool hits closer than HitIsolationCut, other layer
  double absoluteFieldTesla = 0.0;                      ///< |Bz| at the pool centroid, for the pT
  double radiusAtMaxPtCm = 0.0;                         ///< radius of a track with pT = MaxComboPT
  int maxHitsPerCombination = 0;
  int minHitsPerCombination = 0;
};

std::tuple<edm4hep::TrackCollection, edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection>
CombinatorialTrackFinder::operator()(const edm4hep::TrackerHitPlaneCollection& hits,
                                     const edm4hep::EventHeaderCollection& headers,
                                     const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks) const {
  edm4hep::TrackCollection trackCollection;
  edm4hep::TrackerHitPlaneCollection outputHits;
  edm4hep::TrackerHitSimTrackerHitLinkCollection outputLinks;

  unsigned int eventNumber = 0;
  if (!headers.empty())
    eventNumber = headers[0].getEventNumber();
  debug() << "Event " << eventNumber << ": " << hits.size() << " hits" << endmsg;

  if (hits.size() < 3) { // a track needs at least 3 hits
    m_tracksPerEvent += 0;
    return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks));
  }
  try {
    EventSearch search(*this, hits, trackCollection, outputHits, recoSimLinks, outputLinks);
    search.run();
  } catch (const std::exception& exception) {
    error() << "Exception during track finding: " << exception.what() << endmsg;
    m_tracksPerEvent += 0;
    return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks));
  }
  m_tracksPerEvent += trackCollection.size();
  if (msgLevel(MSG::DEBUG))
    printTrackSummary(trackCollection);
  return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks));
}

// DEBUG after each event: one line per track, from its state at the first hit (the truth
// comparison belongs in the analysis)
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
        const dd4hep::Position firstHit(position[0] / kMmPerCm, position[1] / kMmPerCm, position[2] / kMmPerCm);
        bField = m_field.magneticField(firstHit).z() / dd4hep::tesla;
      }
      const double pT = (std::abs(state.omega) > 1e-12) ? k4rt::ptFromOmega(bField, state.omega) : 0.0;
      debug() << "  track " << i + 1 << ": " << track.trackerHits_size() << " hits, chi2/ndf "
              << track.getChi2() / std::max(1, track.getNdf()) << ", pT " << pT << " GeV, tanLambda " << state.tanLambda
              << endmsg;
      break;
    }
  }
}

CombinatorialTrackFinder::EventSearch::EventSearch(const CombinatorialTrackFinder& finderIn,
                                                   const edm4hep::TrackerHitPlaneCollection& hitsIn,
                                                   edm4hep::TrackCollection& outputTracksIn,
                                                   edm4hep::TrackerHitPlaneCollection& outputHitsIn,
                                                   const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks,
                                                   edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinksIn)
    : finder(finderIn), hits(hitsIn), outputTracks(outputTracksIn), outputHits(outputHitsIn),
      outputLinks(outputLinksIn), hasSurface(hitsIn.size(), false), hitGeometry(hitsIn.size()),
      used(hitsIn.size(), false), pairCoincidence(finderIn.m_requirePairCoincidence.value()) {
  for (const auto& link : recoSimLinks)
    recoToSimHits[hitKey(link.getFrom())].push_back(link.getTo());

  // what the search needs of every hit, computed once
  for (size_t i = 0; i < hits.size(); ++i) {
    const auto hit = hits[i];
    if (!finder.surfaceOf(hit.getCellID()))
      continue;
    hasSurface[i] = true;
    const auto& position = hit.getPosition(); // mm
    k4rt::HitGeometry& geometry = hitGeometry[i];
    geometry.x = position[0];
    geometry.y = position[1];
    geometry.z = position[2];
    geometry.transverseRadius = std::sqrt(position[0] * position[0] + position[1] * position[1]);
    geometry.absoluteZ = std::abs(position[2]);
    geometry.azimuth = std::atan2(position[1], position[0]);
    geometry.layerKey = finder.layerKeyOf(hit.getCellID());
    geometry.region = k4rt::regionOfLayer(geometry.layerKey, &finder.m_regionOfLayerKey);
  }
}

void CombinatorialTrackFinder::EventSearch::run() {
  if (finder.m_strategyChoice == Strategy::SeedAndExtend)
    runSeedAndExtend();
  else
    runExhaustive();
}

void CombinatorialTrackFinder::EventSearch::runExhaustive() {
  while (true) {
    ++trackNumber;
    if (!collectPool())
      break;
    prepareExhaustiveIteration();
    TrackCandidate best;
    if (!searchCombinations(best))
      break;
    if (finder.m_doCrowdedLayerHitSelection)
      selectCrowdedLayerHits(best);
    if (rejectedAsNeighbour(best))
      continue;
    removeOutliers(best);
    buildTrack(best);
  }
}

void CombinatorialTrackFinder::EventSearch::runSeedAndExtend() {
  ++trackNumber;
  if (!collectPool())
    return;
  // all candidates at once, from the hits usable at the start
  std::vector<std::pair<size_t, int>> pool;
  pool.reserve(poolHits.size());
  for (size_t hitIndex : poolHits)
    pool.emplace_back(hitIndex, hitGeometry[hitIndex].layerKey);
  const std::vector<TrackCandidate> candidates = finder.m_seedFinder->run(hits, pool);

  // stored best first, each only if all its hits are still free
  size_t next = 0;
  while (true) {
    while (next < candidates.size()) {
      bool anyHitUsed = false;
      for (size_t hitIndex : candidates[next].hits) {
        if (used[hitIndex]) {
          anyHitUsed = true;
          break;
        }
      }
      if (!anyHitUsed)
        break;
      ++next;
    }
    if (next == candidates.size())
      break;
    buildTrack(candidates[next]);
    ++next;
    ++trackNumber;
    if (!collectPool()) // the stop rule of the exhaustive search: too few hits or layers left
      break;
  }
}

// the unused hits on a known surface and below MaxHitEdepKeV, and their layers inner -> outer;
// false (stop) with fewer than 3 hits or layers
bool CombinatorialTrackFinder::EventSearch::collectPool() {
  int remainingHits = 0;
  for (size_t i = 0; i < used.size(); ++i) {
    if (!used[i])
      remainingHits++;
  }
  if (remainingHits < 3)
    return false;

  double edepThresholdGeV = -1.0; // MaxHitEdepKeV (EDM4hep edep in GeV); -1 = off
  if (finder.m_maxHitEdepKeV.value() > 0.0)
    edepThresholdGeV = finder.m_maxHitEdepKeV.value() * 1e-6;
  poolHits.clear();
  hitsByLayer.clear();
  for (size_t i = 0; i < hits.size(); ++i) {
    if (used[i] || !hasSurface[i])
      continue;
    if (edepThresholdGeV > 0.0 && hits[i].getEDep() > edepThresholdGeV)
      continue;
    poolHits.push_back(i);
    hitsByLayer[hitGeometry[i].layerKey].push_back(i);
  }
  if (poolHits.size() < 3)
    return false;

  // layer keys by layer number (innermost first), then by key
  layerOrder.clear();
  for (const auto& [layerKey, layerHits] : hitsByLayer)
    layerOrder.push_back(layerKey);
  std::sort(layerOrder.begin(), layerOrder.end(), [](int left, int right) {
    const int leftLayer = k4rt::keyWithoutRegion(left);
    const int rightLayer = k4rt::keyWithoutRegion(right);
    if (leftLayer != rightLayer)
      return leftLayer < rightLayer;
    return left < right;
  });
  return layerOrder.size() >= 3;
}

// Exhaustive: field and pT cap, energy-deposit scores, size range of a combination, isolation
void CombinatorialTrackFinder::EventSearch::prepareExhaustiveIteration() {
  // field at the pool centroid: MaxComboPT as a radius cap for the k = 3 ranking
  absoluteFieldTesla = std::max(std::abs(fieldAtCentroid(poolHits)), 1e-6);
  radiusAtMaxPtCm = k4rt::radiusFromPt(absoluteFieldTesla, finder.m_maxPt.value());

  computeEdepScores();

  int availableHits = 0; // after MaxHitsPerCompositeID
  for (int layerKey : layerOrder)
    availableHits += static_cast<int>(hitsByLayer[layerKey].size());
  maxHitsPerCombination = std::min(availableHits, finder.m_maxCombinatorialHits.value());
  minHitsPerCombination = std::max(3, finder.m_commonConfig.minTrackHits);

  if (finder.m_hitIsolationCut.value() > 0.0)
    findIsolationNeighbours();
}

// Energy-deposit score: the right hit of a crowded layer usually deposits like the hits of the
// single-hit layers (the baseline), delta rays and secondaries much more or less.
//   score = distance of edep from [baseline, 2 x baseline] (keV, lower = better)
// MaxHitsPerCompositeID keeps only the best-scoring hits of crowded layers.
void CombinatorialTrackFinder::EventSearch::computeEdepScores() {
  std::vector<double> edepKeV(hits.size(), 0.0);
  for (size_t hitIndex : poolHits)
    edepKeV[hitIndex] = hits[hitIndex].getEDep() * 1e6; // GeV -> keV

  // baseline: mean edep of the single-hit layers, else the median of all hits
  double baselineKeV = 0.0;
  int singleHitLayers = 0;
  for (const auto& [layerKey, layerHits] : hitsByLayer) {
    if (layerHits.size() == 1) {
      baselineKeV += edepKeV[layerHits[0]];
      ++singleHitLayers;
    }
  }
  if (singleHitLayers > 0) {
    baselineKeV /= static_cast<double>(singleHitLayers);
  } else {
    std::vector<double> sortedEdeps;
    sortedEdeps.reserve(poolHits.size());
    for (size_t hitIndex : poolHits)
      sortedEdeps.push_back(edepKeV[hitIndex]);
    std::sort(sortedEdeps.begin(), sortedEdeps.end());
    const size_t middle = sortedEdeps.size() / 2;
    if (sortedEdeps.size() % 2 == 0)
      baselineKeV = 0.5 * (sortedEdeps[middle - 1] + sortedEdeps[middle]);
    else
      baselineKeV = sortedEdeps[middle];
  }

  edepScore.assign(hits.size(), 0.0);
  if (!(finder.m_edepScoreWeight.value() > 0.0))
    return;
  for (size_t hitIndex : poolHits) {
    const double below = std::max(0.0, baselineKeV - edepKeV[hitIndex]);
    const double above = std::max(0.0, edepKeV[hitIndex] - 2.0 * baselineKeV);
    edepScore[hitIndex] = below + above;
  }
  const int maxHitsPerLayer = finder.m_maxHitsPerLayer.value();
  if (maxHitsPerLayer <= 0)
    return;
  for (auto& [layerKey, layerHits] : hitsByLayer) {
    if (static_cast<int>(layerHits.size()) <= maxHitsPerLayer)
      continue;
    std::sort(layerHits.begin(), layerHits.end(),
              [this](size_t left, size_t right) { return edepScore[left] < edepScore[right]; });
    layerHits.resize(maxHitsPerLayer);
  }
}

// HitIsolationCut: for every pool hit, the pool hits of other layers closer than the cut (cm)
// (same-layer hits, e.g. secondaries next to the track, do not count)
void CombinatorialTrackFinder::EventSearch::findIsolationNeighbours() {
  const double cut = finder.m_hitIsolationCut.value();
  isolationNeighbours.assign(hits.size(), {});
  for (size_t hitIndex : poolHits) {
    const auto position = hits[hitIndex].getPosition();
    const Eigen::Vector3d reference(position[0] / kMmPerCm, position[1] / kMmPerCm, position[2] / kMmPerCm);
    for (size_t otherIndex : poolHits) {
      if (otherIndex == hitIndex || hitGeometry[otherIndex].layerKey == hitGeometry[hitIndex].layerKey)
        continue;
      const auto otherPosition = hits[otherIndex].getPosition();
      const Eigen::Vector3d other(otherPosition[0] / kMmPerCm, otherPosition[1] / kMmPerCm,
                                  otherPosition[2] / kMmPerCm);
      if ((other - reference).norm() < cut)
        isolationNeighbours[hitIndex].push_back(otherIndex);
    }
  }
}

// Bz (T) at the centroid of the hits
double CombinatorialTrackFinder::EventSearch::fieldAtCentroid(const std::vector<size_t>& hitIndices) const {
  const int numberOfHits = static_cast<int>(hitIndices.size());
  double sumX = 0.0, sumY = 0.0, sumZ = 0.0; // cm
  for (size_t hitIndex : hitIndices) {
    const auto position = hits[hitIndex].getPosition();
    sumX += position[0] / kMmPerCm;
    sumY += position[1] / kMmPerCm;
    sumZ += position[2] / kMmPerCm;
  }
  const dd4hep::Position centroid(sumX / numberOfHits, sumY / numberOfHits, sumZ / numberOfHits);
  return finder.m_field.magneticField(centroid).z() / dd4hep::tesla;
}

// ============================================================================================
// Exhaustive: the best combination
// ============================================================================================

// all C(N,k) combinations for k = MaxCombinatorialHits down to MinTrackHits, after road, isolation
// and quality cuts, ranked by their circle fit
bool CombinatorialTrackFinder::EventSearch::searchCombinations(TrackCandidate& best) {
  k4rt::RoadCutConfig roadCutConfig;
  roadCutConfig.maxConsecDeltaPhi = finder.m_maxConsecutiveDeltaPhi.value();
  roadCutConfig.maxComboPhiSpread = finder.m_maxCombinationPhiSpread.value();
  roadCutConfig.maxConsecHitDist = finder.m_maxConsecutiveHitDistance.value();
  roadCutConfig.maxPairHitDist = finder.m_commonConfig.maxPairHitDistanceCm;
  roadCutConfig.maxComboTanLambdaDev = finder.m_maxCombinationTanLambdaDeviation.value();
  const k4rt::RoadCuts roadCuts(roadCutConfig, hitGeometry);

  HitCluster allHits; // all pool hits, layer by layer from the inside out
  allHits.layerKeys = layerOrder;
  for (int layerKey : layerOrder)
    allHits.hitsPerLayer.push_back(hitsByLayer.at(layerKey));
  const std::vector<HitCluster> clusters = buildClusters(allHits);

  // pair coincidence: hits with a coincident partner first, then all hits if a doublet stays
  // incomplete; otherwise one pass over all hits
  std::vector<HitCluster> firstPassClusters;
  if (pairCoincidence)
    firstPassClusters.push_back(pairCoincidenceCluster(allHits));
  const int numberOfDoublets = static_cast<int>(finder.m_doubletLayerPairs.size());
  SearchResult result;
  const int numberOfPasses = pairCoincidence ? 2 : 1;
  for (int pass = 0; pass < numberOfPasses; ++pass) {
    if (pass == 1 && result.found && completeDoublets(result.best.hits) >= numberOfDoublets)
      break;
    const bool coincidentHitsOnly = (pairCoincidence && pass == 0);
    for (const HitCluster& cluster : (coincidentHitsOnly ? firstPassClusters : clusters))
      searchCluster(cluster, roadCuts, result);
  }

  if (result.found) {
    best = result.best;
    return true;
  }
  if (result.fallback.numberOfHits() > 0) {
    best = result.fallback; // every combination was above MaxComboPT
    return true;
  }
  return false;
}

// pair coincidence, first pass: only hits with a coincident hit in the other layer of their doublet
// (same region and subdetector)
HitCluster CombinatorialTrackFinder::EventSearch::pairCoincidenceCluster(const HitCluster& allHits) const {
  HitCluster cluster;
  cluster.layerKeys = allHits.layerKeys;
  cluster.hitsPerLayer.assign(allHits.layerKeys.size(), std::vector<size_t>{});
  for (size_t layerIndex = 0; layerIndex < allHits.layerKeys.size(); ++layerIndex) {
    for (size_t hitIndex : allHits.hitsPerLayer[layerIndex]) {
      // partner layer key: the same key with the layer number swapped
      const int layerKey = hitGeometry[hitIndex].layerKey;
      const int layerNumber = k4rt::layerNumberOfKey(layerKey);
      const auto partner = finder.m_partnerLayer.find(layerNumber);
      if (partner == finder.m_partnerLayer.end())
        continue; // not in a doublet
      const int absolutePartnerKey = std::abs(layerKey) - layerNumber + partner->second;
      const int partnerKey = (layerKey < 0) ? -absolutePartnerKey : absolutePartnerKey;
      const auto partnerHits = hitsByLayer.find(partnerKey);
      if (partnerHits == hitsByLayer.end())
        continue;
      for (size_t partnerIndex : partnerHits->second) {
        if (coincident(hitIndex, partnerIndex)) {
          cluster.hitsPerLayer[layerIndex].push_back(hitIndex);
          break;
        }
      }
    }
  }
  return cluster;
}

// DoHitClustering: hits of neighbouring layers within the road cuts are joined (union-find), each
// cluster searched on its own; without clustering one cluster with all hits
std::vector<HitCluster> CombinatorialTrackFinder::EventSearch::buildClusters(const HitCluster& allHits) const {
  const int numberOfLayers = static_cast<int>(allHits.layerKeys.size());
  if (!(finder.m_doHitClustering.value() && numberOfLayers >= 2))
    return {allHits};

  // union-find on the pool position of a hit; clusters come in the (hash) order of their
  // representative position, the order the results always had
  std::vector<size_t> poolPosition(hits.size(), 0);
  for (size_t position = 0; position < poolHits.size(); ++position)
    poolPosition[poolHits[position]] = position;
  UnionFind groups(poolHits.size());
  const double phiCut = finder.m_maxConsecutiveDeltaPhi.value();
  const double distanceCut = finder.m_maxConsecutiveHitDistance.value(); // cm
  joinLayers(groups, allHits, poolPosition, 1, phiCut, distanceCut);
  joinLayers(groups, allHits, poolPosition, 2, (phiCut > 0.0) ? phiCut * 2.0 : -1.0,
             (distanceCut > 0.0) ? distanceCut * 2.0 : -1.0); // across one missing layer, cuts doubled

  std::unordered_map<size_t, std::vector<size_t>> hitsOfGroup; // representative -> hits
  for (int layerIndex = 0; layerIndex < numberOfLayers; ++layerIndex) {
    for (size_t hitIndex : allHits.hitsPerLayer[layerIndex])
      hitsOfGroup[groups.find(poolPosition[hitIndex])].push_back(hitIndex);
  }

  std::vector<HitCluster> clusters; // with at least minHitsPerCombination layers, inner -> outer
  for (auto& [representative, groupHits] : hitsOfGroup) {
    std::unordered_map<int, std::vector<size_t>> hitsOfLayer;
    for (size_t hitIndex : groupHits)
      hitsOfLayer[hitGeometry[hitIndex].layerKey].push_back(hitIndex);
    if (static_cast<int>(hitsOfLayer.size()) < minHitsPerCombination)
      continue;
    HitCluster cluster;
    for (int layerKey : layerOrder) {
      auto layerHits = hitsOfLayer.find(layerKey);
      if (layerHits == hitsOfLayer.end())
        continue;
      cluster.layerKeys.push_back(layerKey);
      cluster.hitsPerLayer.push_back(std::move(layerHits->second));
    }
    clusters.push_back(std::move(cluster));
  }
  return clusters;
}

// joins hits of layers `step` apart (1: neighbours, 2: across one missing layer) within the phi
// and distance limits
void CombinatorialTrackFinder::EventSearch::joinLayers(UnionFind& groups, const HitCluster& allHits,
                                                       const std::vector<size_t>& poolPosition, int step,
                                                       double phiLimit, double distanceLimit) const {
  const int numberOfLayers = static_cast<int>(allHits.hitsPerLayer.size());
  for (int layerIndex = 0; layerIndex + step < numberOfLayers; ++layerIndex) {
    for (size_t innerHit : allHits.hitsPerLayer[layerIndex]) {
      const auto innerPosition = hits[innerHit].getPosition();
      const double innerAzimuth = std::atan2(innerPosition[1], innerPosition[0]);
      for (size_t outerHit : allHits.hitsPerLayer[layerIndex + step]) {
        const auto outerPosition = hits[outerHit].getPosition();
        if (phiLimit > 0.0) {
          const double deltaPhi = k4rt::wrapAngle(std::atan2(outerPosition[1], outerPosition[0]) - innerAzimuth);
          if (std::abs(deltaPhi) > phiLimit)
            continue;
        }
        if (distanceLimit > 0.0) {
          const double dx = outerPosition[0] - innerPosition[0];
          const double dy = outerPosition[1] - innerPosition[1];
          const double dz = outerPosition[2] - innerPosition[2];
          if (std::sqrt(dx * dx + dy * dy + dz * dz) * kCmPerMm > distanceLimit)
            continue;
        }
        groups.unite(poolPosition[innerHit], poolPosition[outerHit]);
      }
    }
  }
}

// all combinations of one cluster: k = most hits down to the fewest, every choice of k layers and
// every choice of one hit per chosen layer
void CombinatorialTrackFinder::EventSearch::searchCluster(const HitCluster& cluster, const k4rt::RoadCuts& roadCuts,
                                                          SearchResult& result) const {
  const int numberOfLayers = static_cast<int>(cluster.layerKeys.size());
  const int largestSize = std::min(maxHitsPerCombination, numberOfLayers);
  std::vector<size_t> combination;
  for (int size = largestSize; size >= minHitsPerCombination; --size) {
    std::vector<int> chosenLayers(size); // first choice 0, 1, ..., size-1
    for (int i = 0; i < size; ++i)
      chosenLayers[i] = i;
    do {
      // the hits of each chosen layer (a layer without hits gives no combination)
      std::vector<const std::vector<size_t>*> hitsOfChosenLayer(size);
      bool anyLayerEmpty = false;
      for (int i = 0; i < size; ++i) {
        hitsOfChosenLayer[i] = &cluster.hitsPerLayer[chosenLayers[i]];
        if (hitsOfChosenLayer[i]->empty()) {
          anyLayerEmpty = true;
          break;
        }
      }
      if (anyLayerEmpty)
        continue;
      std::vector<int> position(size, 0); // the hit taken in each chosen layer
      combination.assign(size, 0);
      do {
        for (int i = 0; i < size; ++i)
          combination[i] = (*hitsOfChosenLayer[i])[position[i]];
        if (passesPreFitCuts(combination, roadCuts))
          offerCombination(combination, result);
      } while (nextHitChoice(position, hitsOfChosenLayer));
    } while (nextLayerChoice(chosenLayers, numberOfLayers));
  }
}

// the cuts before the circle fit, cheapest first (the counters record the rejections)
bool CombinatorialTrackFinder::EventSearch::passesPreFitCuts(const std::vector<size_t>& combination,
                                                             const k4rt::RoadCuts& roadCuts) const {
  if (!roadCuts.passesRoadCuts(combination)) {
    ++finder.m_rejectedRoadCuts;
    return false;
  }
  if (!roadCuts.passesDistanceCuts(combination)) {
    ++finder.m_rejectedDistanceCuts;
    return false;
  }
  if (!isolated(combination))
    return false; // HitIsolationCut
  if (!roadCuts.passesTanLambdaCut(combination)) {
    ++finder.m_rejectedTanLambdaCut;
    return false;
  }
  if (pairCoincidence && completeDoublets(combination) < kMinCompleteDoublets)
    return false;
  return true;
}

// no hit of the combination has an isolation neighbour outside the combination
bool CombinatorialTrackFinder::EventSearch::isolated(const std::vector<size_t>& combination) const {
  if (finder.m_hitIsolationCut.value() <= 0.0)
    return true;
  for (size_t hitIndex : combination) {
    for (size_t neighbour : isolationNeighbours[hitIndex]) {
      if (std::find(combination.begin(), combination.end(), neighbour) == combination.end())
        return false;
    }
  }
  return true;
}

// two hits of a doublet are coincident within |dphi| < 0.15, |deta| < 0.10 and PairCoincidenceMaxDist
bool CombinatorialTrackFinder::EventSearch::coincident(size_t first, size_t second) const {
  const double maxDistanceMm = finder.m_pairCoincidenceMaxDistance.value() * kMmPerCm;
  const auto firstPosition = hits[first].getPosition();
  const auto secondPosition = hits[second].getPosition();
  const double dx = secondPosition[0] - firstPosition[0];
  const double dy = secondPosition[1] - firstPosition[1];
  const double dz = secondPosition[2] - firstPosition[2];
  if ((dx * dx + dy * dy + dz * dz) > (maxDistanceMm * maxDistanceMm))
    return false;
  const double deltaPhi = k4rt::wrapAngle(std::atan2(secondPosition[1], secondPosition[0]) -
                                          std::atan2(firstPosition[1], firstPosition[0]));
  if (std::abs(deltaPhi) > kPairTolerancePhi)
    return false;
  const double firstRadius = std::sqrt(firstPosition[0] * firstPosition[0] + firstPosition[1] * firstPosition[1]);
  const double secondRadius = std::sqrt(secondPosition[0] * secondPosition[0] + secondPosition[1] * secondPosition[1]);
  if (firstRadius < 1e-6 || secondRadius < 1e-6)
    return false;
  const double firstEta = std::asinh(firstPosition[2] / firstRadius);
  const double secondEta = std::asinh(secondPosition[2] / secondRadius);
  return std::abs(secondEta - firstEta) <= kPairToleranceEta;
}

// complete doublets: hits in both layers of a doublet, in the same region, coincident
int CombinatorialTrackFinder::EventSearch::completeDoublets(const std::vector<size_t>& combination) const {
  int completeCount = 0;
  for (const auto& [firstLayer, secondLayer] : finder.m_doubletLayerPairs) {
    bool complete = false;
    for (size_t firstHit : combination) {
      const int firstKey = hitGeometry[firstHit].layerKey;
      if (k4rt::layerNumberOfKey(firstKey) != firstLayer)
        continue;
      for (size_t secondHit : combination) {
        const int secondKey = hitGeometry[secondHit].layerKey;
        if (secondHit == firstHit || k4rt::layerNumberOfKey(secondKey) != secondLayer)
          continue;
        if (k4rt::regionFromKey(secondKey) == k4rt::regionFromKey(firstKey) && coincident(firstHit, secondHit)) {
          complete = true;
          break;
        }
      }
      if (complete)
        break;
    }
    if (complete)
      ++completeCount;
  }
  return completeCount;
}

// fit a combination; past the chi2 and pT cuts it replaces the best one if better; above
// MaxComboPT it can only become the fallback (near-straight tracks)
void CombinatorialTrackFinder::EventSearch::offerCombination(const std::vector<size_t>& combination,
                                                             SearchResult& result) const {
  ++finder.m_combinationsTried;
  const std::optional<CircleFit> fit = k4rt::fitCircle(hits, combination, finder.m_commonConfig.sigmaHitDefaultCm);
  if (!fit)
    return;
  if (fit->chi2ndf >= finder.m_maxCircleChi2NDF.value()) {
    ++finder.m_rejectedChi2Cut;
    return;
  }
  const int numberOfHits = static_cast<int>(combination.size());
  const double pT = k4rt::ptFromRadius(absoluteFieldTesla, fit->radius); // GeV
  if (pT > finder.m_maxPt.value()) {
    ++finder.m_rejectedPtCuts;
    if (isBetterFallback(numberOfHits, *fit, result.fallback)) {
      result.fallback.hits = combination;
      result.fallback.circle = *fit;
      result.fallback.edepScoreKeV = 0.0;
    }
    return;
  }
  if (finder.m_commonConfig.minPtGeV > 0.0 && pT < finder.m_commonConfig.minPtGeV) {
    ++finder.m_rejectedPtCuts; // below MinComboPT: no fallback
    return;
  }
  const double score = meanEdepScore(combination);
  if (isBetterCandidate(numberOfHits, *fit, score, result.best)) {
    result.best.hits = combination;
    result.best.circle = *fit;
    result.best.edepScoreKeV = score;
    result.found = true;
  }
}

// mean energy-deposit score (keV; 0 when EdepScoreWeight is off)
double CombinatorialTrackFinder::EventSearch::meanEdepScore(const std::vector<size_t>& combination) const {
  if (!(finder.m_edepScoreWeight.value() > 0.0))
    return 0.0;
  double sum = 0.0;
  for (size_t hitIndex : combination)
    sum += edepScore[hitIndex];
  return sum / static_cast<double>(combination.size());
}

// ranking (lower score is better): more hits always win; k >= 4: chi2/ndf + w * edep / EdepNormKeV;
// k = 3 (chi2 is 0): w * edep / EdepNormKeV - min(R, R_cap) / R_cap (straighter and MIP-like
// wins), with w = 0 the larger radius
bool CombinatorialTrackFinder::EventSearch::isBetterCandidate(int numberOfHits, const CircleFit& fit,
                                                              double edepScoreKeV, const TrackCandidate& best) const {
  if (numberOfHits != best.numberOfHits())
    return numberOfHits > best.numberOfHits();
  const double weight = finder.m_edepScoreWeight.value();
  const double edepNorm = std::max(finder.m_edepNormKeV.value(), 1e-9);
  if (numberOfHits == 3) {
    if (!(weight > 0.0))
      return fit.radius > best.circle.radius;
    const double radiusTerm = std::min(fit.radius, radiusAtMaxPtCm) / radiusAtMaxPtCm; // saturates at R_cap
    const double bestRadiusTerm = std::min(best.circle.radius, radiusAtMaxPtCm) / radiusAtMaxPtCm;
    const double score = weight * edepScoreKeV / edepNorm - radiusTerm;
    const double bestScore = weight * best.edepScoreKeV / edepNorm - bestRadiusTerm;
    return score < bestScore;
  }
  const double score = fit.chi2ndf + (weight > 0.0 ? weight * edepScoreKeV / edepNorm : 0.0);
  const double bestScore = best.circle.chi2ndf + (weight > 0.0 ? weight * best.edepScoreKeV / edepNorm : 0.0);
  return score < bestScore;
}

// the best fallback: more hits, then lower chi2/ndf (k >= 4) or lower pT (k = 3)
bool CombinatorialTrackFinder::EventSearch::isBetterFallback(int numberOfHits, const CircleFit& fit,
                                                             const TrackCandidate& fallback) {
  if (numberOfHits != fallback.numberOfHits())
    return numberOfHits > fallback.numberOfHits();
  if (numberOfHits == 3)
    return fit.radius < fallback.circle.radius;
  return fit.chi2ndf < fallback.circle.chi2ndf;
}

// ============================================================================================
// Exhaustive: refinement of the chosen combination
// ============================================================================================

// DoCrowdedLayerHitSelection: layers with one hit are "clean"; in each crowded layer the hit with
// the lowest 0.75 * geometric term + 0.25 * edep score / EdepNormKeV is chosen, the geometric term
// being the distance to the circle through >= 3 clean hits (/ SigmaHitDefault), or with 2 clean
// hits the kink to their line (/ MaxConsecDeltaPhi). Kept if the refit has a lower chi2/ndf.
void CombinatorialTrackFinder::EventSearch::selectCrowdedLayerHits(TrackCandidate& best) const {
  int cleanLayers = 0;
  int crowdedLayers = 0;
  for (const auto& [layerKey, layerHits] : hitsByLayer) {
    if (layerHits.size() == 1)
      ++cleanLayers;
    else if (layerHits.size() > 1)
      ++crowdedLayers;
  }
  if (crowdedLayers < 1)
    return;

  if (cleanLayers >= 3) { // circle through the clean hits
    const std::vector<size_t> cleanHits = cleanLayerHits();
    if (cleanHits.size() < 3)
      return;
    const std::optional<CircleFit> cleanCircle =
        k4rt::fitCircle(hits, cleanHits, finder.m_commonConfig.sigmaHitDefaultCm);
    if (!cleanCircle || !(cleanCircle->radius > 0))
      return;
    const double distanceNorm = std::max(finder.m_sigmaHitDefault.value(), 1e-9); // cm
    auto distanceToCircle = [&](size_t hitIndex) {
      const auto position = hits[hitIndex].getPosition(); // mm
      const double dx = position.x * kCmPerMm - cleanCircle->centreX;
      const double dy = position.y * kCmPerMm - cleanCircle->centreY;
      return std::abs(std::sqrt(dx * dx + dy * dy) - cleanCircle->radius) / distanceNorm;
    };
    std::vector<size_t> guided = cleanHits;
    if (addBestHitOfCrowdedLayers(distanceToCircle, guided))
      applyGuidedChoice(guided, best);
  } else if (cleanLayers == 2 && cleanLayers + crowdedLayers >= 3) { // line through the two clean hits
    const std::vector<size_t> cleanHits = cleanLayerHits();
    if (cleanHits.size() != 2)
      return;
    const auto innerPosition = hits[cleanHits[0]].getPosition(); // mm
    const auto outerPosition = hits[cleanHits[1]].getPosition();
    const double innerX = innerPosition.x * kCmPerMm, innerY = innerPosition.y * kCmPerMm,
                 innerZ = innerPosition.z * kCmPerMm;
    const double outerX = outerPosition.x * kCmPerMm, outerY = outerPosition.y * kCmPerMm,
                 outerZ = outerPosition.z * kCmPerMm;
    const double lineX = outerX - innerX, lineY = outerY - innerY, lineZ = outerZ - innerZ;
    const double lineLength = std::sqrt(lineX * lineX + lineY * lineY + lineZ * lineZ);
    if (!(lineLength > 1e-9))
      return;                                                                        // the two clean hits coincide
    const double kinkNorm = std::max(finder.m_maxConsecutiveDeltaPhi.value(), 1e-4); // rad
    auto kinkAngle = [&](size_t hitIndex) {               // angle at the inner clean hit between the line and the hit
      const auto position = hits[hitIndex].getPosition(); // mm
      const double toHitX = position.x * kCmPerMm - innerX;
      const double toHitY = position.y * kCmPerMm - innerY;
      const double toHitZ = position.z * kCmPerMm - innerZ;
      const double toHitLength = std::sqrt(toHitX * toHitX + toHitY * toHitY + toHitZ * toHitZ);
      double kink = 0.0;
      if (toHitLength > 1e-9) {
        double cosine = (lineX * toHitX + lineY * toHitY + lineZ * toHitZ) / (lineLength * toHitLength);
        cosine = std::max(-1.0, std::min(1.0, cosine));
        kink = std::acos(cosine);
      }
      return kink / kinkNorm;
    };
    std::vector<size_t> guided = cleanHits;
    if (addBestHitOfCrowdedLayers(kinkAngle, guided))
      applyGuidedChoice(guided, best);
  }
}

// the hit of every single-hit layer, by layer key
std::vector<size_t> CombinatorialTrackFinder::EventSearch::cleanLayerHits() const {
  std::vector<size_t> cleanHits;
  for (const auto& [layerKey, layerHits] : hitsByLayer) {
    if (layerHits.size() == 1)
      cleanHits.push_back(layerHits[0]);
  }
  return cleanHits;
}

// adds the lowest-score hit of every crowded layer (by layer key) to `guided`; false if a crowded
// layer has no candidate
bool CombinatorialTrackFinder::EventSearch::addBestHitOfCrowdedLayers(
    const std::function<double(size_t)>& geometricTerm, std::vector<size_t>& guided) const {
  const double edepNorm = std::max(finder.m_edepNormKeV.value(), 1e-9); // keV
  for (const auto& [layerKey, layerHits] : hitsByLayer) {
    if (layerHits.size() <= 1)
      continue;
    size_t bestHit = SIZE_MAX;
    double bestScore = std::numeric_limits<double>::max();
    for (size_t hitIndex : layerHits) {
      const double score =
          kGuidedGeometricWeight * geometricTerm(hitIndex) + kGuidedEdepWeight * (edepScore[hitIndex] / edepNorm);
      if (score < bestScore) {
        bestScore = score;
        bestHit = hitIndex;
      }
    }
    if (bestHit == SIZE_MAX)
      return false;
    guided.push_back(bestHit);
  }
  return true;
}

// the guided choice replaces the best one if it differs and its refit has a lower chi2/ndf
void CombinatorialTrackFinder::EventSearch::applyGuidedChoice(const std::vector<size_t>& guided,
                                                              TrackCandidate& best) const {
  const std::optional<CircleFit> fit = k4rt::fitCircle(hits, guided, finder.m_commonConfig.sigmaHitDefaultCm);
  if (!fit || fit->chi2ndf > finder.m_maxCircleChi2NDF.value())
    return;
  bool differentChoice = false;
  for (size_t hitIndex : guided) {
    if (std::find(best.hits.begin(), best.hits.end(), hitIndex) == best.hits.end()) {
      differentChoice = true;
      break;
    }
  }
  if (differentChoice && fit->chi2ndf < best.circle.chi2ndf) {
    best.hits = guided;
    best.circle = *fit;
  }
}

// Neighbour-track cuts (delta-ray shadows): a combination with a hit within NeighbourTrackMaxDist
// of an accepted track needs NeighbourTrackMinHits and NeighbourTrackMaxChi2NDF; otherwise its hits
// near the accepted track are used up and true is returned (next iteration)
bool CombinatorialTrackFinder::EventSearch::rejectedAsNeighbour(const TrackCandidate& best) {
  if (!(finder.m_neighbourTrackMaxDist.value() > 0.0) || acceptedTrackHitPositions.empty())
    return false;
  const double distanceMm = finder.m_neighbourTrackMaxDist.value() * kMmPerCm;
  const double distanceSquared = distanceMm * distanceMm;
  bool isNeighbour = false;
  for (size_t hitIndex : best.hits) {
    if (nearAcceptedTrack(hitIndex, distanceSquared)) {
      isNeighbour = true;
      break;
    }
  }
  if (!isNeighbour)
    return false;
  const bool enoughHits = (best.numberOfHits() >= finder.m_neighbourTrackMinHits.value());
  const bool goodFit = (best.circle.chi2ndf <= finder.m_neighbourTrackMaxChi2NDF.value());
  if (enoughHits && goodFit)
    return false;
  for (size_t hitIndex : best.hits) {
    if (nearAcceptedTrack(hitIndex, distanceSquared))
      used[hitIndex] = true;
  }
  return true;
}

bool CombinatorialTrackFinder::EventSearch::nearAcceptedTrack(size_t hitIndex, double distanceSquaredMm2) const {
  const auto position = hits[hitIndex].getPosition(); // mm
  for (const auto& accepted : acceptedTrackHitPositions) {
    const double dx = position.x - accepted[0];
    const double dy = position.y - accepted[1];
    const double dz = position.z - accepted[2];
    if (dx * dx + dy * dy + dz * dz < distanceSquaredMm2)
      return true;
  }
  return false;
}

// drop the worst hit while its residual is above OutlierSigma * SigmaHitDefault, and refit
void CombinatorialTrackFinder::EventSearch::removeOutliers(TrackCandidate& best) const {
  const double outlierCutCm = finder.m_outlierSigma.value() * finder.m_sigmaHitDefault.value();
  for (int iteration = 0; iteration < finder.m_maxOutlierIterations.value(); ++iteration) {
    if (best.numberOfHits() <= minHitsPerCombination)
      break;
    size_t worst = 0;
    double worstResidual = 0.0;
    for (size_t i = 0; i < best.circle.residuals.size(); ++i) {
      if (best.circle.residuals[i] > worstResidual) {
        worstResidual = best.circle.residuals[i];
        worst = i;
      }
    }
    if (worstResidual <= outlierCutCm)
      break;
    best.hits.erase(best.hits.begin() + worst);
    best.circle.residuals.erase(best.circle.residuals.begin() + worst);
    if (best.numberOfHits() < minHitsPerCombination)
      break;
    const std::optional<CircleFit> refit = k4rt::fitCircle(hits, best.hits, finder.m_commonConfig.sigmaHitDefaultCm);
    if (refit)
      best.circle = *refit;
  }
}

// ============================================================================================
// output: the track and its states
// ============================================================================================

// the candidate as an EDM4hep track: hits marked used, helix states, hits and reco-sim links copied
void CombinatorialTrackFinder::EventSearch::buildTrack(const TrackCandidate& candidate) {
  std::vector<edm4hep::TrackerHitPlane> trackHits;
  trackHits.reserve(candidate.hits.size());
  for (size_t hitIndex : candidate.hits) {
    used[hitIndex] = true;
    trackHits.push_back(hits[hitIndex]);
  }
  if (finder.m_neighbourTrackMaxDist.value() > 0.0) { // for the neighbour-track cuts
    for (const auto& hit : trackHits) {
      const auto& position = hit.getPosition(); // mm
      acceptedTrackHitPositions.push_back({position.x, position.y, position.z});
    }
  }

  auto track = outputTracks.create();
  const std::optional<edm4hep::TrackState> firstHitState = addHelixStates(candidate, trackHits, track);
  // SeedAndExtend: the kink chi2 (multiple scattering); else the circle-fit chi2 (ndf = hits - 3)
  if (candidate.kinkNdf > 0) {
    track.setChi2(static_cast<float>(candidate.kinkChi2));
    track.setNdf(candidate.kinkNdf);
  } else {
    const int fitNdf = std::max(1, static_cast<int>(trackHits.size()) - 3);
    track.setChi2(static_cast<float>(candidate.circle.chi2ndf * fitNdf));
    track.setNdf(fitNdf);
  }
  if (candidate.looper)
    track.setType(1);
  for (const auto& hit : trackHits)
    copyHitToOutput(hit, track);
  if (finder.m_ipState == IPState::SolenoidCrossing)
    addSolenoidCrossingState(candidate, trackHits, firstHitState, track);

  finder.m_hitsPerTrack += trackHits.size();
}

// the helix as AtFirstHit, AtLastHit and (IPStateMode = Helix) AtIP; returns the AtFirstHit state
std::optional<edm4hep::TrackState>
CombinatorialTrackFinder::EventSearch::addHelixStates(const TrackCandidate& candidate,
                                                      const std::vector<edm4hep::TrackerHitPlane>& trackHits,
                                                      edm4hep::MutableTrack& track) const {
  struct StoredState {
    StatePoint where;
    int location;
    bool wanted;
  };
  const StoredState storedStates[] = {
      {StatePoint::FirstHit, edm4hep::TrackState::AtFirstHit, true},
      {StatePoint::LastHit, edm4hep::TrackState::AtLastHit, true},
      {StatePoint::Origin, edm4hep::TrackState::AtIP, finder.m_ipState == IPState::Helix},
  };
  std::optional<edm4hep::TrackState> firstHitState;
  for (const auto& stored : storedStates) {
    if (!stored.wanted)
      continue;
    edm4hep::TrackState state;
    if (!makeEdm4hepStateFromCircle(candidate.circle, trackHits, stored.where, state)) {
      finder.warning() << "Could not build track state " << stored.location << " of track " << trackNumber
                       << " (degenerate hits)" << endmsg;
      continue;
    }
    state.location = stored.location;
    track.addToTrackStates(state);
    if (stored.location == edm4hep::TrackState::AtFirstHit)
      firstHitState = state;
  }
  return firstHitState;
}

// IPStateMode = SolenoidCrossing (tracks with all hits outside the solenoid): rotation sense and
// tanLambda from the AtFirstHit state, pT from the radius and the field at the hits' centroid
void CombinatorialTrackFinder::EventSearch::addSolenoidCrossingState(
    const TrackCandidate& candidate, const std::vector<edm4hep::TrackerHitPlane>& trackHits,
    const std::optional<edm4hep::TrackState>& firstHitState, edm4hep::MutableTrack& track) const {
  if (!finder.allHitsOutsideSolenoid(trackHits) || !firstHitState || firstHitState->omega == 0.0f)
    return;
  OuterHelix outer;
  outer.centreXCm = candidate.circle.centreX;
  outer.centreYCm = candidate.circle.centreY;
  outer.radiusCm = candidate.circle.radius;
  outer.clockwise = (firstHitState->omega > 0.0f); // EDM4hep: omega > 0 = clockwise from +z
  outer.ptGeV = k4rt::ptFromRadius(fieldAtCentroid(candidate.hits), candidate.circle.radius);
  outer.tanLambda = firstHitState->tanLambda;
  const auto& position = trackHits[0].getPosition(); // the innermost hit anchors the crossing
  outer.innermostHitCm = Eigen::Vector3d(position[0] / kMmPerCm, position[1] / kMmPerCm, position[2] / kMmPerCm);
  outer.innermostHitInEndcap = (finder.regionOf(trackHits[0].getCellID()) != 0);
  const std::optional<edm4hep::TrackState> atIP = finder.solenoidCrossingState(outer);
  if (atIP)
    track.addToTrackStates(*atIP);
  finder.m_solenoidCrossingBuilt += atIP.has_value();
}

// a copy of the input hit in the output, attached to the track, with all its reco-sim links (a hit
// can have several SimHits, e.g. a muon and a delta ray in one cell)
void CombinatorialTrackFinder::EventSearch::copyHitToOutput(const edm4hep::TrackerHitPlane& hit,
                                                            edm4hep::MutableTrack& track) {
  auto outputHit = outputHits.create();
  outputHit.setCellID(hit.getCellID());
  outputHit.setTime(hit.getTime());
  outputHit.setEDep(hit.getEDep());
  outputHit.setEDepError(hit.getEDepError());
  outputHit.setPosition(hit.getPosition());
  outputHit.setCovMatrix(hit.getCovMatrix());
  outputHit.setDu(hit.getDu());
  outputHit.setDv(hit.getDv());
  outputHit.setU(hit.getU());
  outputHit.setV(hit.getV());
  const auto simHits = recoToSimHits.find(hitKey(hit));
  if (simHits != recoToSimHits.end()) {
    for (const auto& simHit : simHits->second) {
      auto link = outputLinks.create();
      link.setFrom(outputHit);
      link.setTo(simHit);
      link.setWeight(1.0f);
    }
  }
  track.addToTrackerHits(outputHit);
}

// ============================================================================================
// detector surfaces and cellID decoding
// ============================================================================================

const dd4hep::rec::Surface* CombinatorialTrackFinder::surfaceOf(uint64_t cellID) const {
  const auto entry = m_surfaceMap->find(cellID);
  if (entry == m_surfaceMap->end())
    return nullptr;
  return dynamic_cast<const dd4hep::rec::Surface*>(entry->second);
}

int CombinatorialTrackFinder::regionOf(uint64_t cellID) const {
  if (!m_hasRegionField)
    return 0; // single-region readout: all barrel
  return m_bitFieldCoder->get(cellID, m_regionFieldName.value());
}

int CombinatorialTrackFinder::layerNumberOf(uint64_t cellID) const {
  return m_bitFieldCoder->get(cellID, m_layerFieldName.value());
}

// barrel 100*system + layer, +endcap 10000 + that, -endcap -(10000 + that); the system term only
// when several subdetectors share one readout
int CombinatorialTrackFinder::layerKeyOf(uint64_t cellID) const {
  int key = layerNumberOf(cellID);
  if (m_useSystemInLayerKey)
    key += 100 * m_bitFieldCoder->get(cellID, "system");
  const int region = regionOf(cellID);
  if (region == 0)
    return key;
  if (region == 1)
    return k4rt::kEndcapKeyOffset + key;
  return -k4rt::kEndcapKeyOffset - key;
}

bool CombinatorialTrackFinder::readoutHasField(const std::string& name) const {
  for (const auto& field : m_bitFieldCoder->fields()) {
    if (field.name() == name)
      return true;
  }
  return false;
}

std::map<int, std::vector<const dd4hep::rec::Surface*>> CombinatorialTrackFinder::groupSurfacesByLayer() const {
  std::map<int, std::vector<const dd4hep::rec::Surface*>> surfacesByLayer;
  for (const auto& surface : m_surfaces) {
    const dd4hep::DetElement detElement = surface->detElement();
    if (!detElement.isValid())
      continue;
    surfacesByLayer[layerKeyOf(detElement.volumeID())].push_back(surface); // same key as the hits
  }
  if (!msgLevel(MSG::INFO))
    return surfacesByLayer;
  info() << "Layer distribution:" << endmsg;
  for (const auto& [layerKey, surfaces] : surfacesByLayer) {
    const int region = k4rt::regionFromKey(layerKey);
    const char* regionName = (region == 0) ? "Barrel" : (region > 0 ? "Positive Endcap" : "Negative Endcap");
    info() << "  " << regionName << " Layer " << k4rt::keyWithoutRegion(layerKey) << " (ID=" << layerKey
           << "): " << surfaces.size() << " surfaces" << endmsg;
  }
  return surfacesByLayer;
}

// ============================================================================================
// the solenoid
// ============================================================================================

// The field is followed outward from near the IP in 1 mm steps until it is no longer the central
// field (Bz moved more than half-way to zero or beyond: a return field or a 0 T gap), then found
// to 0.01 mm by bisection. alongZ: along z at half the solenoid radius (on the axis compensating
// solenoids can change the field near the IP).
double CombinatorialTrackFinder::findSolenoidBoundary(bool alongZ) const {
  const std::string tag = alongZ ? "[SolScanZ] |Z|" : "[SolScan] R";
  const double stepMm = 1.0;
  const double startMm = 10.0;
  const double scanLimitMm = 10000.0; // beyond any FCC-ee solenoid

  const double bzCentral = fieldAlongScan(startMm, alongZ);
  if (std::abs(bzCentral) < 1e-6) {
    warning() << "  " << tag << ": no field near the IP; the solenoid boundary is set to " << scanLimitMm << " mm"
              << endmsg;
    return scanLimitMm;
  }
  double lowMm = -1.0; // coarse scan: the 1 mm interval where the central field ends
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
  for (int iteration = 0; iteration < 40; ++iteration) { // bisection down to 0.01 mm
    const double middleMm = 0.5 * (lowMm + highMm);
    if (outsideCentralField(middleMm, alongZ, bzCentral))
      highMm = middleMm;
    else
      lowMm = middleMm;
    if (highMm - lowMm < 0.01)
      break;
  }
  return 0.5 * (lowMm + highMm);
}

double CombinatorialTrackFinder::fieldAlongScan(double distanceMm, bool alongZ) const {
  const double distanceCm = distanceMm / kMmPerCm; // DD4hep positions are in cm
  dd4hep::Position position(distanceCm, 0.0, 0.0);
  if (alongZ) // at half the solenoid radius (the radial scan runs first)
    position = dd4hep::Position(0.5 * m_solenoidRadiusMm / kMmPerCm, 0.0, distanceCm);
  return m_field.magneticField(position).z() / dd4hep::tesla;
}

bool CombinatorialTrackFinder::outsideCentralField(double distanceMm, bool alongZ, double bzCentral) const {
  return std::abs(fieldAlongScan(distanceMm, alongZ) - bzCentral) > 0.5 * std::abs(bzCentral);
}

bool CombinatorialTrackFinder::isInsideSolenoid(const Eigen::Vector3d& positionCm) const {
  const double transverseRadius = std::sqrt(positionCm.x() * positionCm.x() + positionCm.y() * positionCm.y());
  const double absoluteZ = std::abs(positionCm.z());
  return (transverseRadius < m_solenoidRadiusMm / kMmPerCm) && (absoluteZ < m_solenoidHalfLengthMm / kMmPerCm);
}

bool CombinatorialTrackFinder::allHitsOutsideSolenoid(const std::vector<edm4hep::TrackerHitPlane>& trackHits) const {
  for (const auto& hit : trackHits) {
    const auto& position = hit.getPosition();
    const Eigen::Vector3d positionCm(position[0] / kMmPerCm, position[1] / kMmPerCm, position[2] / kMmPerCm);
    if (isInsideSolenoid(positionCm))
      return false;
  }
  return true;
}

// The solenoid boundary = a cylinder at R_sol closed by disks at |z| = Z_sol; the outer circle is
// followed back from the innermost hit to it, then the inner helix (inner field, same pT and
// tanLambda) is attached and its perigee returned.
std::optional<edm4hep::TrackState> CombinatorialTrackFinder::solenoidCrossingState(const OuterHelix& outer) const {
  const double solenoidHalfLengthCm = m_solenoidHalfLengthMm / kMmPerCm;
  const double centreX = outer.centreXCm;
  const double centreY = outer.centreYCm;
  const double radius = outer.radiusCm;
  // backward on the outer circle (+1: angle grows, clockwise track); also the reversed curvature of
  // the inner helix
  const double innerSign = outer.clockwise ? 1.0 : -1.0;

  // 1. the crossing of the boundary (cm)
  Eigen::Vector3d crossing;
  if (!outer.innermostHitInEndcap) {
    if (!solenoidBarrelCrossing(outer, true, crossing))
      return std::nullopt;
  } else {
    // endcap hit: the end disk, at the arc s back from the hit with z(s) = +-Z_sol
    if (std::abs(outer.tanLambda) < 1e-6)
      return std::nullopt;
    const double hitZ = outer.innermostHitCm.z();
    const double diskZ = (hitZ < 0.0) ? -solenoidHalfLengthCm : +solenoidHalfLengthCm;
    const double arcToDisk = (hitZ - diskZ) / outer.tanLambda; // hitZ - tanLambda * s = diskZ
    if (arcToDisk < 0.0) {
      // the hit is already at |z| < Z_sol: the barrel crossing instead
      if (!solenoidBarrelCrossing(outer, false, crossing))
        return std::nullopt;
    } else {
      if (arcToDisk > 2.0 * M_PI * radius)
        return std::nullopt;
      const double angleOfHit = std::atan2(outer.innermostHitCm.y() - centreY, outer.innermostHitCm.x() - centreX);
      const double angleAtDisk = angleOfHit + innerSign * (arcToDisk / radius);
      crossing =
          Eigen::Vector3d(centreX + radius * std::cos(angleAtDisk), centreY + radius * std::sin(angleAtDisk), diskZ);
    }
  }

  // 2. the inner field, at the origin
  const double innerField = m_field.magneticField(dd4hep::Position(0.0, 0.0, 0.0)).z() / dd4hep::tesla;
  if (std::abs(innerField) < 1e-6)
    return std::nullopt;

  // 3. inner circle: same pT and tangent at the crossing, bent the other way, so its centre lies
  //    on the outward radial line of the outer circle, R_inner = pT / (0.3 B_inner) beyond it
  const double radialX = (crossing.x() - centreX) / radius;
  const double radialY = (crossing.y() - centreY) / radius;
  const double innerRadius = k4rt::radiusFromPt(innerField, outer.ptGeV); // cm
  const double innerCentreX = crossing.x() + innerRadius * radialX;
  const double innerCentreY = crossing.y() + innerRadius * radialY;

  // 4. perigee: d0 = |inner centre| - R_inner, direction perpendicular to the line to the centre
  const double d0Cm = std::sqrt(innerCentreX * innerCentreX + innerCentreY * innerCentreY) - innerRadius;
  const double phi0 = k4rt::wrapAngle(std::atan2(innerCentreY, innerCentreX) - innerSign * M_PI / 2.0);
  const double omega = innerSign / (innerRadius * kMmPerCm); // 1/mm, signed

  // 5. z at the perigee: back from the crossing, the shorter way (< pi)
  const double angleOfCrossing = std::atan2(crossing.y() - innerCentreY, crossing.x() - innerCentreX);
  const double angleOfPerigee = std::atan2(-innerCentreY, -innerCentreX);
  double backwardAngle = (omega > 0) ? k4rt::wrapAngleZeroToTwoPi(angleOfCrossing - angleOfPerigee)
                                     : k4rt::wrapAngleZeroToTwoPi(angleOfPerigee - angleOfCrossing);
  if (backwardAngle > M_PI)
    backwardAngle = 2.0 * M_PI - backwardAngle;
  const double z0Cm = crossing.z() - outer.tanLambda * (innerRadius * backwardAngle);

  return stateWithDefaultCovariance(d0Cm * kMmPerCm, phi0, omega, z0Cm * kMmPerCm, outer.tanLambda,
                                    edm4hep::TrackState::AtIP);
}

bool CombinatorialTrackFinder::solenoidBarrelCrossing(const OuterHelix& outer, bool tryOtherCrossing,
                                                      Eigen::Vector3d& crossingCm) const {
  const double solenoidRadiusCm = m_solenoidRadiusMm / kMmPerCm;
  const double centreX = outer.centreXCm;
  const double centreY = outer.centreYCm;
  const double radius = outer.radiusCm;
  const double centreDistance = std::sqrt(centreX * centreX + centreY * centreY); // from the IP
  if (centreDistance < 1e-6)
    return false;
  if (radius <= std::abs(centreDistance - solenoidRadiusCm) || radius >= centreDistance + solenoidRadiusCm)
    return false; // the circles miss each other
  double crossingX[2];
  double crossingY[2];
  circleCylinderCrossings(centreX, centreY, radius, centreDistance, solenoidRadiusCm, crossingX, crossingY);

  // the two crossings in the order reached going back from the hit
  const double hitX = outer.innermostHitCm.x();
  const double hitY = outer.innermostHitCm.y();
  const double hitZ = outer.innermostHitCm.z();
  const double angleOfHit = std::atan2(hitY - centreY, hitX - centreX);
  double backwardAngle[2];
  for (int i = 0; i < 2; ++i)
    backwardAngle[i] = backwardArcAngle(centreX, centreY, angleOfHit, outer.clockwise, crossingX[i], crossingY[i]);
  const int first = (backwardAngle[0] <= backwardAngle[1]) ? 0 : 1;
  const int order[2] = {first, 1 - first};

  // the first (with tryOtherCrossing also the other) crossing on the IP side of the hit in z
  const int crossingsToTry = tryOtherCrossing ? 2 : 1;
  for (int attempt = 0; attempt < crossingsToTry; ++attempt) {
    const int crossingIndex = order[attempt];
    const double z = hitZ - outer.tanLambda * (radius * backwardAngle[crossingIndex]);
    if (!(std::abs(z) > std::abs(hitZ) + kCrossingZToleranceCm)) {
      crossingCm = Eigen::Vector3d(crossingX[crossingIndex], crossingY[crossingIndex], z);
      return true;
    }
  }
  return false;
}

// ============================================================================================
// formula helpers (declared at the top)
// ============================================================================================

namespace {

// The fitted helix (circle + straight line of z versus transverse arc length) as an EDM4hep track
// state. Geometric, so valid for either field sign: the particle moves outward (from the hit
// nearest the origin), omega > 0 = clockwise seen from +z, and d0, phi, z0 are given at the circle
// point closest to the reference point (innermost hit, outermost hit or origin), EDM4hep
// convention: pca = ref + (-d0 sin(phi), d0 cos(phi), z0), centre = pca + (1/omega)(sin(phi), -cos(phi)).
// Circle in cm, state in mm. false (state unchanged) for degenerate input.
bool makeEdm4hepStateFromCircle(const CircleFit& circle, const std::vector<edm4hep::TrackerHitPlane>& trackHits,
                                StatePoint where, edm4hep::TrackState& state) {
  const double centreX = circle.centreX * kMmPerCm; // mm
  const double centreY = circle.centreY * kMmPerCm;
  const double radius = circle.radius * kMmPerCm;
  if (radius <= 0.0 || trackHits.size() < 2)
    return false;

  // 1. hits sorted by distance from the origin (the order an outward-moving particle crosses them)
  struct HitPoint {
    double x, y, z;   // mm
    double distance;  // from the origin, mm
    double varianceZ; // mm^2
  };
  std::vector<HitPoint> points;
  for (const auto& hit : trackHits) {
    const auto& position = hit.getPosition();
    HitPoint point;
    point.x = position[0];
    point.y = position[1];
    point.z = position[2];
    point.distance = std::sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
    point.varianceZ = hit.getCovMatrix(edm4hep::Cartesian::z, edm4hep::Cartesian::z);
    points.push_back(point);
  }
  std::sort(points.begin(), points.end(),
            [](const HitPoint& left, const HitPoint& right) { return left.distance < right.distance; });

  // 2. sense of rotation seen from +z (+1 = counter-clockwise): sum of the cross products of
  //    consecutive hits about the centre
  double crossSum = 0.0;
  for (size_t i = 0; i + 1 < points.size(); ++i) {
    const double fromCentreX = points[i].x - centreX;
    const double fromCentreY = points[i].y - centreY;
    const double nextFromCentreX = points[i + 1].x - centreX;
    const double nextFromCentreY = points[i + 1].y - centreY;
    crossSum += fromCentreX * nextFromCentreY - fromCentreY * nextFromCentreX;
  }
  if (crossSum == 0.0)
    return false;
  const double sense = (crossSum > 0.0) ? 1.0 : -1.0;

  // 3. reference point and the closest circle point (pca)
  double referenceX = 0.0, referenceY = 0.0, referenceZ = 0.0;
  size_t nearestHit = 0; // the hit next to the reference point, to count the arc
  if (where == StatePoint::FirstHit) {
    referenceX = points.front().x;
    referenceY = points.front().y;
    referenceZ = points.front().z;
  } else if (where == StatePoint::LastHit) {
    nearestHit = points.size() - 1;
    referenceX = points.back().x;
    referenceY = points.back().y;
    referenceZ = points.back().z;
  }
  const double centreToReferenceX = referenceX - centreX;
  const double centreToReferenceY = referenceY - centreY;
  const double centreToReference =
      std::sqrt(centreToReferenceX * centreToReferenceX + centreToReferenceY * centreToReferenceY);
  if (centreToReference < 1e-6)
    return false; // reference at the centre: no closest point
  const double anglePca = std::atan2(centreToReferenceY, centreToReferenceX);
  const double pcaX = centreX + radius * centreToReferenceX / centreToReference;
  const double pcaY = centreY + radius * centreToReferenceY / centreToReference;

  // 4. phi = direction of travel at the pca (radius vector turned by +-90 degrees), omega > 0 for
  //    clockwise motion, d0 from pca - ref = (-d0 sin(phi), d0 cos(phi))
  const double phi = k4rt::wrapAngle(anglePca + sense * M_PI / 2.0);
  const double omega = -sense / radius; // 1/mm
  const double d0 = -(pcaX - referenceX) * std::sin(phi) + (pcaY - referenceY) * std::cos(phi);

  // 5. angle turned along the travel, hit by hit (no 2 pi jumps), and the arc length of every hit
  //    from the pca
  std::vector<double> turned(points.size());
  double previousAngle = std::atan2(points[0].y - centreY, points[0].x - centreX);
  turned[0] = 0.0;
  for (size_t i = 1; i < points.size(); ++i) {
    const double angle = std::atan2(points[i].y - centreY, points[i].x - centreX);
    const double step = k4rt::wrapAngle(sense * (angle - previousAngle));
    turned[i] = turned[i - 1] + step;
    previousAngle = angle;
  }
  const double angleNearest = std::atan2(points[nearestHit].y - centreY, points[nearestHit].x - centreX);
  const double pcaStep = k4rt::wrapAngle(sense * (anglePca - angleNearest));
  const double turnedPca = turned[nearestHit] + pcaStep;

  // 6. weighted line fit z = zPca + tanLambda * s (weights 1 / sigma_z^2; a disk hit measures z
  //    exactly, its sigma_z^2 is 0 or rounding noise: floor (1 um)^2), about the weighted means
  //    (no cancellation when one weight is much larger)
  const double minimumVarianceZ = 1e-6; // mm^2
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
  const double meanS = sumWS / sumW;
  const double meanZ = sumWZ / sumW;
  double sumWDSDS = 0.0, sumWDSDZ = 0.0;
  for (size_t i = 0; i < points.size(); ++i) {
    const double deviationS = arcLength[i] - meanS;
    sumWDSDS += weight[i] * deviationS * deviationS;
    sumWDSDZ += weight[i] * deviationS * (points[i].z - meanZ);
  }
  if (!(sumWDSDS / sumW > 1e-6))
    return false; // the hits do not spread along the arc (weighted spread below 1 um)
  const double tanLambda = sumWDSDZ / sumWDSDS;
  const double zPca = meanZ - tanLambda * meanS;

  // 7. covariance: d0, phi, omega from the circle covariance (gradients in x0, y0, R, cm), with
  //    lower limits; z0, tanLambda from the line fit
  const double centreToReferenceSquaredCm2 = centreToReference * centreToReference / 100.0;
  const Eigen::Vector3d gradientD0(-centreToReferenceX / centreToReference, -centreToReferenceY / centreToReference,
                                   -1.0); // cm/cm
  const Eigen::Vector3d gradientPhi(centreToReferenceY / 10.0 / centreToReferenceSquaredCm2,
                                    -centreToReferenceX / 10.0 / centreToReferenceSquaredCm2, 0.0); // rad/cm
  const Eigen::Vector3d gradientOmega(0.0, 0.0, 1.0 / (circle.radius * circle.radius * 10.0));      // (1/mm)/cm
  const Eigen::Matrix3d& covariance = circle.covariance;
  const double varianceD0 = std::max(gradientD0.dot(covariance * gradientD0) * 100.0, 0.01); // mm^2
  const double variancePhi = std::max(gradientPhi.dot(covariance * gradientPhi), 0.01);
  const double varianceOmega = std::max(gradientOmega.dot(covariance * gradientOmega), std::pow(0.1 * omega, 2));

  using Parameter = edm4hep::TrackParams;
  for (int row = 0; row < 6; ++row)
    for (int column = 0; column <= row; ++column)
      state.setCovMatrix(0.0f, static_cast<Parameter>(row), static_cast<Parameter>(column));
  state.setCovMatrix(static_cast<float>(varianceD0), Parameter::d0, Parameter::d0);
  state.setCovMatrix(static_cast<float>(variancePhi), Parameter::phi, Parameter::phi);
  state.setCovMatrix(static_cast<float>(varianceOmega), Parameter::omega, Parameter::omega);
  // line fit, with S = sum w (s - meanS)^2: var(tanLambda) = 1/S, var(z0) = 1/sumW + meanS^2/S,
  // cov(z0, tanLambda) = -meanS/S
  state.setCovMatrix(static_cast<float>(1.0 / sumW + meanS * meanS / sumWDSDS), Parameter::z0, Parameter::z0);
  state.setCovMatrix(static_cast<float>(1.0 / sumWDSDS), Parameter::tanLambda, Parameter::tanLambda);
  state.setCovMatrix(static_cast<float>(-meanS / sumWDSDS), Parameter::z0, Parameter::tanLambda);

  state.D0 = static_cast<float>(d0);
  state.phi = static_cast<float>(phi);
  state.omega = static_cast<float>(omega);
  state.Z0 = static_cast<float>(zPca - referenceZ);
  state.tanLambda = static_cast<float>(tanLambda);
  state.referencePoint = edm4hep::Vector3f(referenceX, referenceY, referenceZ);
  return true;
}

// track state with these parameters and a default diagonal covariance (D0, Z0 in mm^2,
// phi, tanLambda in rad^2, omega in (1/mm)^2)
edm4hep::TrackState stateWithDefaultCovariance(double d0, double phi, double omega, double z0, double tanLambda,
                                               int location) {
  edm4hep::TrackState state;
  state.D0 = d0;
  state.phi = phi;
  state.omega = omega;
  state.Z0 = z0;
  state.tanLambda = tanLambda;
  state.location = location;
  for (int row = 0; row < 6; ++row)
    for (int column = 0; column <= row; ++column)
      state.setCovMatrix(0.0, static_cast<edm4hep::TrackParams>(row), static_cast<edm4hep::TrackParams>(column));
  state.setCovMatrix(1.0, edm4hep::TrackParams::d0, edm4hep::TrackParams::d0);
  state.setCovMatrix(0.01, edm4hep::TrackParams::phi, edm4hep::TrackParams::phi);
  state.setCovMatrix(1e-10, edm4hep::TrackParams::omega, edm4hep::TrackParams::omega);
  state.setCovMatrix(1.0, edm4hep::TrackParams::z0, edm4hep::TrackParams::z0);
  state.setCovMatrix(0.01, edm4hep::TrackParams::tanLambda, edm4hep::TrackParams::tanLambda);
  return state;
}

// The two crossings of a circle (centre at the distance d from the origin, radius R) with the
// solenoid cylinder (radius Rs): on the line from the origin to the centre at
// h = (Rs^2 + d^2 - R^2) / (2 d), and +-l = sqrt(Rs^2 - h^2) sideways.
void circleCylinderCrossings(double centreX, double centreY, double radius, double centreDistance,
                             double solenoidRadius, double crossingX[2], double crossingY[2]) {
  const double along =
      (solenoidRadius * solenoidRadius + centreDistance * centreDistance - radius * radius) / (2.0 * centreDistance);
  const double sideways = std::sqrt(std::max(0.0, solenoidRadius * solenoidRadius - along * along));
  const double unitX = centreX / centreDistance, unitY = centreY / centreDistance;
  const double perpendicularX = -unitY, perpendicularY = unitX;
  crossingX[0] = along * unitX + sideways * perpendicularX;
  crossingY[0] = along * unitY + sideways * perpendicularY;
  crossingX[1] = along * unitX - sideways * perpendicularX;
  crossingY[1] = along * unitY - sideways * perpendicularY;
}

// angle on the circle from the hit (at angleOfHit about the centre) back to the point, in [0, 2pi),
// against the direction of travel (the track came from the IP side)
double backwardArcAngle(double centreX, double centreY, double angleOfHit, bool clockwise, double pointX,
                        double pointY) {
  const double angleOfPoint = std::atan2(pointY - centreY, pointX - centreX);
  return clockwise ? k4rt::wrapAngleZeroToTwoPi(angleOfPoint - angleOfHit)
                   : k4rt::wrapAngleZeroToTwoPi(angleOfHit - angleOfPoint);
}

// next choice of k of the numberOfLayers layers (increasing indices, lexicographic order); false
// after the last
bool nextLayerChoice(std::vector<int>& chosen, int numberOfLayers) {
  const int size = static_cast<int>(chosen.size());
  int i = size - 1;
  while (i >= 0 && chosen[i] == numberOfLayers - size + i)
    --i;
  if (i < 0)
    return false;
  ++chosen[i];
  for (int j = i + 1; j < size; ++j)
    chosen[j] = chosen[j - 1] + 1;
  return true;
}

// next choice of one hit per chosen layer, like an odometer (the last layer turns first); false
// after the last
bool nextHitChoice(std::vector<int>& position, const std::vector<const std::vector<size_t>*>& hitsOfLayer) {
  for (int i = static_cast<int>(position.size()) - 1; i >= 0; --i) {
    ++position[i];
    if (position[i] < static_cast<int>(hitsOfLayer[i]->size()))
      return true;
    position[i] = 0;
  }
  return false;
}

} // namespace
