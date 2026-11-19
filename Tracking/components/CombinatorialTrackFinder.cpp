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
#include <cmath>
#include <cstdint>
#include <exception>
#include <iomanip>

DECLARE_COMPONENT(CombinatorialTrackFinder)

using CTF::kMmPerCm;
using CTF::TrackCandidate;

namespace {

// key of a hit: (collectionID, index); merged input collections all start their indices at 0
template <typename Hit>
uint64_t hitKey(const Hit& hit) {
  const auto id = hit.id();
  return (static_cast<uint64_t>(id.collectionID) << 32) | static_cast<uint32_t>(id.index);
}

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
  if (setupCellIDDecoding().isFailure() || setupSurfaces().isFailure())
    return StatusCode::FAILURE;
  setupField();
  if (readSearchOptions().isFailure())
    return StatusCode::FAILURE;
  if (m_strategyChoice == Strategy::Exhaustive)
    setupExhaustiveSearch();
  if (m_strategyChoice == Strategy::SeedAndExtend && setupSeedFinder().isFailure())
    return StatusCode::FAILURE;
  info() << "Found " << m_surfaces.size() << " detector surfaces in " << m_surfacesByLayer.size() << " layers"
         << endmsg;
  if (m_ipState == IPState::SolenoidCrossing)
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

// surface map (one detector or merged), surfaces by layer, region of every layer key, and the
// layer geometry
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
    if (inserted && region != CTF::regionFromKey(layerKey))
      ++keysWithOtherRegion;
  }
  if (keysWithOtherRegion > 0)
    info() << keysWithOtherRegion << " layer key(s) whose region differs from the key-range rule; "
           << "using the readout region field" << endmsg;

  m_layerGeometry.build(m_surfacesByLayer);
  info() << "Layer geometry (inside-out by closest approach to the IP):" << endmsg;
  for (const auto& layer : m_layerGeometry.layers()) {
    info() << "  layer key " << std::setw(6) << layer.layerKey << "  " << std::setw(7) << CTF::kindName(layer.kind)
           << "  r=[" << std::fixed << std::setprecision(2) << layer.rMin_cm << "," << layer.rMax_cm << "]"
           << "  z=[" << layer.zMin_cm << "," << layer.zMax_cm << "] cm"
           << "  dMin=" << layer.closestApproachCm << "  nSurf=" << layer.numberOfSurfaces << std::defaultfloat
           << std::setprecision(6) << endmsg;
  }
  for (const auto& problem : m_layerGeometry.validate())
    warning() << "LayerGeometry: " << problem << endmsg;
  return StatusCode::SUCCESS;
}

// the field map: Bz at the IP (pT of the tracks) and at any point (the searches)
void CombinatorialTrackFinder::setupField() {
  m_field = m_detector->field();
  const dd4hep::Direction fieldAtCentre = m_field.magneticField(dd4hep::Position(0, 0, 0));
  info() << "Magnetic field at the detector centre: (" << fieldAtCentre.x() / dd4hep::tesla << ", "
         << fieldAtCentre.y() / dd4hep::tesla << ", " << fieldAtCentre.z() / dd4hep::tesla << ") T" << endmsg;

  // adapter: Bz of the field map at a point (cm), for the tools without DD4hep
  m_bzAt = [this](double xCm, double yCm, double zCm) {
    return m_field.magneticField(dd4hep::Position(xCm, yCm, zCm)).z() / dd4hep::tesla;
  };
}

// the options shared by both strategies; text options -> enums; DoubletLayers checked
StatusCode CombinatorialTrackFinder::readSearchOptions() {
  m_commonConfig = CTF::CommonSearchConfig{.minTrackHits = m_minTrackHits.value(),
                                           .minPtGeV = m_minPt.value(),
                                           .maxPairHitDistanceCm = m_maxPairHitDistance.value(),
                                           .bFieldTesla = m_field.magneticField(dd4hep::Position(0, 0, 0)).z() /
                                                          dd4hep::tesla, // Bz at the IP
                                           .sigmaHitDefaultCm = m_sigmaHitDefault.value(),
                                           .maxHitEdepKeV = m_maxHitEdepKeV.value()};

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
  return StatusCode::SUCCESS;
}

// Exhaustive: the settings, filled once from the options, and the search
void CombinatorialTrackFinder::setupExhaustiveSearch() {
  // the doublets of the pair coincidence, two layer numbers each
  std::vector<std::pair<int, int>> doubletLayerPairs;
  std::map<int, int> partnerLayer;
  const std::vector<int>& doubletLayers = m_doubletLayers.value();
  for (size_t i = 0; i + 1 < doubletLayers.size(); i += 2) {
    doubletLayerPairs.emplace_back(doubletLayers[i], doubletLayers[i + 1]);
    partnerLayer[doubletLayers[i]] = doubletLayers[i + 1];
    partnerLayer[doubletLayers[i + 1]] = doubletLayers[i];
  }
  CTF::ExhaustiveSearch::Config config{
      .common = m_commonConfig,
      .roadCuts = {.maxConsecutiveDeltaPhi = m_maxConsecutiveDeltaPhi.value(),
                   .maxCombinationPhiSpread = m_maxCombinationPhiSpread.value(),
                   .maxConsecutiveHitDistanceCm = m_maxConsecutiveHitDistance.value(),
                   .maxPairHitDistanceCm = m_commonConfig.maxPairHitDistanceCm,
                   .maxCombinationTanLambdaDeviation = m_maxCombinationTanLambdaDeviation.value()},
      .maxCombinatorialHits = m_maxCombinatorialHits.value(),
      .maxCircleChi2NDF = m_maxCircleChi2NDF.value(),
      .maxPtGeV = m_maxPt.value(),
      .hitIsolationCutCm = m_hitIsolationCut.value(),
      .outlierSigma = m_outlierSigma.value(),
      .maxOutlierIterations = m_maxOutlierIterations.value(),
      .maxHitsPerLayer = m_maxHitsPerLayer.value(),
      .edepScoreWeight = m_edepScoreWeight.value(),
      .edepNormKeV = m_edepNormKeV.value(),
      .doHitClustering = m_doHitClustering.value(),
      .doCrowdedLayerHitSelection = m_doCrowdedLayerHitSelection.value(),
      .requirePairCoincidence = m_requirePairCoincidence.value(),
      .doubletLayerPairs = std::move(doubletLayerPairs),
      .partnerLayer = std::move(partnerLayer),
      .pairCoincidenceMaxDistanceCm = m_pairCoincidenceMaxDistance.value(),
      .neighbourTrackMaxDistanceCm = m_neighbourTrackMaxDist.value(),
      .neighbourTrackMinHits = m_neighbourTrackMinHits.value(),
      .neighbourTrackMaxChi2NDF = m_neighbourTrackMaxChi2NDF.value()};
  m_exhaustiveSearch = std::make_unique<CTF::ExhaustiveSearch>(std::move(config), m_bzAt);
}

// SeedAndExtend: the material of every layer and the seed finder
StatusCode CombinatorialTrackFinder::setupSeedFinder() {
  CTF::SeedLayerPolicy policyChoice = CTF::SeedLayerPolicy::Innermost;
  const std::string& policy = m_seedLayerPolicy.value();
  if (policy == "innermost") {
    policyChoice = CTF::SeedLayerPolicy::Innermost;
  } else if (policy == "outermost") {
    policyChoice = CTF::SeedLayerPolicy::Outermost;
  } else if (policy == "consecutive") {
    policyChoice = CTF::SeedLayerPolicy::Consecutive;
  } else {
    error() << "Unknown SeedLayerPolicy '" << policy << "' (expected innermost, outermost or consecutive)" << endmsg;
    return StatusCode::FAILURE;
  }
  const CTF::SeedFinder::Config seedConfig{.common = m_commonConfig,
                                           .policy = policyChoice,
                                           .maxCandidatesPerLayer = m_maxCandidatesPerLayer.value(),
                                           .shortTrackMaxHits = m_shortTrackMaxHits.value(),
                                           .shadowHitWindow_cm = m_shadowHitWindow.value(),
                                           .kinkMaxChi2Triplet = m_kinkMaxChi2Triplet.value(),
                                           .seedKinkMaxChi2NDF = m_seedKinkMaxChi2NDF.value(),
                                           .trackKinkMaxChi2NDF = m_trackKinkMaxChi2NDF.value()};

  // layer material, for the multiple scattering in the kink chi2
  m_layerGeometry.measureMaterial(*m_detector, m_surfacesByLayer);
  auto& materialLine = info();
  materialLine << "SeedAndExtend: B = " << seedConfig.common.bFieldTesla << " T; layer x/X0 [%]:";
  for (const auto& layer : m_layerGeometry.layers())
    materialLine << " " << layer.layerKey << ":" << std::fixed << std::setprecision(2) << 100.0 * layer.xOverX0;
  materialLine << std::defaultfloat << std::setprecision(6) << endmsg;

  m_seedFinder = std::make_unique<CTF::SeedFinder>(m_layerGeometry, seedConfig);
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

// SolenoidCrossing: the solenoid boundary (radius first: the z scan runs at half of it) and the
// inner field (Bz at the origin)
void CombinatorialTrackFinder::findSolenoid() {
  m_solenoid.radiusMm = findSolenoidBoundary(false);
  info() << "Solenoid Bz boundary (radial): R = " << std::fixed << std::setprecision(4) << m_solenoid.radiusMm
         << " mm  (" << m_solenoid.radiusMm / kMmPerCm << " cm)" << endmsg;
  m_solenoid.halfLengthMm = findSolenoidBoundary(true);
  info() << "Solenoid Bz boundary (endcap): |Z| = " << std::fixed << std::setprecision(4) << m_solenoid.halfLengthMm
         << " mm  (" << m_solenoid.halfLengthMm / kMmPerCm << " cm)" << endmsg;
  m_solenoid.innerFieldTesla = m_field.magneticField(dd4hep::Position(0.0, 0.0, 0.0)).z() / dd4hep::tesla;
}

// ============================================================================================
// the search in one event
// ============================================================================================

// One event: every input hit is described once (CTF::EventHit), the strategy finds the tracks
//   Exhaustive: CTF::ExhaustiveSearch, one track per iteration;
//   SeedAndExtend: the SeedFinder builds all candidates once; they are taken best first, skipping
//   those sharing a hit with an accepted track;
// then each track is written with its states, hits and reco-sim links.
class CombinatorialTrackFinder::EventSearch {
public:
  EventSearch(const CombinatorialTrackFinder& finder, const edm4hep::TrackerHitPlaneCollection& hits,
              edm4hep::TrackCollection& outputTracks, edm4hep::TrackerHitPlaneCollection& outputHits,
              const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks,
              edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinks);

  void run();

private:
  void describeHits();
  std::vector<TrackCandidate> runExhaustive() const;
  std::vector<TrackCandidate> runSeedAndExtend();

  // ---- output ----
  void buildTrack(const TrackCandidate& candidate, int trackNumber);
  std::optional<edm4hep::TrackState> addHelixStates(const TrackCandidate& candidate, int trackNumber,
                                                    edm4hep::MutableTrack& track) const;
  void addSolenoidCrossingState(const TrackCandidate& candidate,
                                const std::optional<edm4hep::TrackState>& firstHitState,
                                edm4hep::MutableTrack& track) const;
  bool allHitsOutsideSolenoid(const TrackCandidate& candidate) const;
  void copyHitToOutput(const edm4hep::TrackerHitPlane& hit, edm4hep::MutableTrack& track);

  // ---- the event ----
  const CombinatorialTrackFinder& m_finder;
  const edm4hep::TrackerHitPlaneCollection& m_inputHits;
  edm4hep::TrackCollection& m_outputTracks;
  edm4hep::TrackerHitPlaneCollection& m_outputHits;
  edm4hep::TrackerHitSimTrackerHitLinkCollection& m_outputLinks;
  std::unordered_map<uint64_t, std::vector<edm4hep::SimTrackerHit>> m_recoToSimHits; ///< hitKey -> its SimTrackerHits
  std::vector<CTF::EventHit> m_hits; ///< every input hit as the searches see it (index = hit index)
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
      const double pT = (std::abs(state.omega) > 1e-12) ? CTF::ptFromOmega(bField, state.omega) : 0.0;
      debug() << "  track " << i + 1 << ": " << track.trackerHits_size() << " hits, chi2/ndf "
              << track.getChi2() / std::max(1, track.getNdf()) << ", pT " << pT << " GeV, tanLambda " << state.tanLambda
              << endmsg;
      break;
    }
  }
}

CombinatorialTrackFinder::EventSearch::EventSearch(const CombinatorialTrackFinder& finder,
                                                   const edm4hep::TrackerHitPlaneCollection& hits,
                                                   edm4hep::TrackCollection& outputTracks,
                                                   edm4hep::TrackerHitPlaneCollection& outputHits,
                                                   const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks,
                                                   edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinks)
    : m_finder(finder), m_inputHits(hits), m_outputTracks(outputTracks), m_outputHits(outputHits),
      m_outputLinks(outputLinks), m_hits(hits.size()) {
  for (const auto& link : recoSimLinks)
    m_recoToSimHits[hitKey(link.getFrom())].push_back(link.getTo());
  describeHits();
}

// what the searches need of every hit, computed once
void CombinatorialTrackFinder::EventSearch::describeHits() {
  const bool needsPlanes = (m_finder.m_strategyChoice == Strategy::SeedAndExtend);
  for (size_t i = 0; i < m_inputHits.size(); ++i) {
    const auto hit = m_inputHits[i];
    const dd4hep::rec::Surface* surface = m_finder.surfaceOf(hit.getCellID());
    if (!surface)
      continue;
    CTF::EventHit& description = m_hits[i];
    description.usable = true;
    description.cellID = hit.getCellID();
    const auto& position = hit.getPosition(); // mm
    description.positionMm = position;
    description.transverseRadiusMm = std::sqrt(position[0] * position[0] + position[1] * position[1]);
    description.absoluteZMm = std::abs(position[2]);
    description.azimuth = std::atan2(position[1], position[0]);
    description.layerKey = m_finder.layerKeyOf(hit.getCellID());
    description.region = m_finder.regionOfLayerKey(description.layerKey);
    description.duMm = hit.getDu();
    description.dvMm = hit.getDv();
    description.varianceZMm2 = hit.getCovMatrix(edm4hep::Cartesian::z, edm4hep::Cartesian::z);
    description.edepGeV = hit.getEDep();
    if (needsPlanes) { // the sensor plane: origin (cm) and unit normal
      const auto surfaceOrigin = surface->origin();
      const auto surfaceNormal = surface->normal(surfaceOrigin);
      description.hasPlane = true;
      description.planeOriginCm = Eigen::Vector3d(surfaceOrigin.x(), surfaceOrigin.y(), surfaceOrigin.z());
      description.planeNormal = Eigen::Vector3d(surfaceNormal.x(), surfaceNormal.y(), surfaceNormal.z()).normalized();
    }
  }
}

void CombinatorialTrackFinder::EventSearch::run() {
  std::vector<TrackCandidate> tracks;
  if (m_finder.m_strategyChoice == Strategy::SeedAndExtend)
    tracks = runSeedAndExtend();
  else
    tracks = runExhaustive();
  for (size_t i = 0; i < tracks.size(); ++i)
    buildTrack(tracks[i], static_cast<int>(i) + 1);
}

// the exhaustive search, its statistics into the Gaudi counters
std::vector<TrackCandidate> CombinatorialTrackFinder::EventSearch::runExhaustive() const {
  CTF::ExhaustiveStatistics statistics;
  std::vector<TrackCandidate> tracks = m_finder.m_exhaustiveSearch->run(m_hits, statistics);
  m_finder.m_combinationsTried += statistics.combinationsFitted;
  m_finder.m_rejectedRoadCuts += statistics.rejectedRoadCuts;
  m_finder.m_rejectedDistanceCuts += statistics.rejectedDistanceCuts;
  m_finder.m_rejectedTanLambdaCut += statistics.rejectedTanLambdaCut;
  m_finder.m_rejectedChi2Cut += statistics.rejectedChi2Cut;
  m_finder.m_rejectedPtCuts += statistics.rejectedPtCuts;
  return tracks;
}

// all candidates at once, from the hits free at the start; taken best first, each only if all its
// hits are still free, until the stop rule of the exhaustive search (too few hits or layers left)
std::vector<TrackCandidate> CombinatorialTrackFinder::EventSearch::runSeedAndExtend() {
  std::vector<TrackCandidate> tracks;
  std::vector<bool> used(m_hits.size(), false);
  const double maxHitEdepKeV = m_finder.m_commonConfig.maxHitEdepKeV;
  CTF::HitPool pool;
  if (!CTF::collectHitPool(m_hits, used, maxHitEdepKeV, pool))
    return tracks;
  const std::vector<TrackCandidate> candidates = m_finder.m_seedFinder->run(m_hits, pool.hits);

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
    tracks.push_back(candidates[next]);
    for (size_t hitIndex : candidates[next].hits)
      used[hitIndex] = true;
    ++next;
    if (!CTF::collectHitPool(m_hits, used, maxHitEdepKeV, pool))
      break;
  }
  return tracks;
}

// ============================================================================================
// output: the track and its states
// ============================================================================================

// the candidate as an EDM4hep track: helix states, hits and reco-sim links copied
void CombinatorialTrackFinder::EventSearch::buildTrack(const TrackCandidate& candidate, int trackNumber) {
  std::vector<edm4hep::TrackerHitPlane> trackHits;
  trackHits.reserve(candidate.hits.size());
  for (size_t hitIndex : candidate.hits)
    trackHits.push_back(m_inputHits[hitIndex]);

  auto track = m_outputTracks.create();
  const std::optional<edm4hep::TrackState> firstHitState = addHelixStates(candidate, trackNumber, track);
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
    track.setType(1); // 1 = looper
  for (const auto& hit : trackHits)
    copyHitToOutput(hit, track);
  if (m_finder.m_ipState == IPState::SolenoidCrossing)
    addSolenoidCrossingState(candidate, firstHitState, track);

  m_finder.m_hitsPerTrack += trackHits.size();
}

// the helix as AtFirstHit, AtLastHit and (IPStateMode = Helix) AtIP; returns the AtFirstHit state
std::optional<edm4hep::TrackState>
CombinatorialTrackFinder::EventSearch::addHelixStates(const TrackCandidate& candidate, int trackNumber,
                                                      edm4hep::MutableTrack& track) const {
  struct StoredState {
    CTF::StatePoint where;
    int location;
    bool wanted;
  };
  const StoredState storedStates[] = {
      {CTF::StatePoint::FirstHit, edm4hep::TrackState::AtFirstHit, true},
      {CTF::StatePoint::LastHit, edm4hep::TrackState::AtLastHit, true},
      {CTF::StatePoint::Origin, edm4hep::TrackState::AtIP, m_finder.m_ipState == IPState::Helix},
  };
  std::optional<edm4hep::TrackState> firstHitState;
  for (const auto& stored : storedStates) {
    if (!stored.wanted)
      continue;
    edm4hep::TrackState state;
    if (!CTF::makeStateFromCircle(candidate.circle, m_hits, candidate.hits, stored.where, state)) {
      m_finder.warning() << "Could not build track state " << stored.location << " of track " << trackNumber
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
    const TrackCandidate& candidate, const std::optional<edm4hep::TrackState>& firstHitState,
    edm4hep::MutableTrack& track) const {
  if (!allHitsOutsideSolenoid(candidate) || !firstHitState || firstHitState->omega == 0.0f)
    return;
  CTF::OuterHelix outer;
  outer.centreXCm = candidate.circle.centreX;
  outer.centreYCm = candidate.circle.centreY;
  outer.radiusCm = candidate.circle.radius;
  outer.clockwise = (firstHitState->omega > 0.0f); // EDM4hep: omega > 0 = clockwise from +z
  outer.ptGeV =
      CTF::ptFromRadius(CTF::fieldAtCentroid(m_hits, candidate.hits, m_finder.m_bzAt), candidate.circle.radius);
  outer.tanLambda = firstHitState->tanLambda;
  const CTF::EventHit& innermostHit = m_hits[candidate.hits[0]]; // anchors the crossing
  const auto& position = innermostHit.positionMm;
  outer.innermostHitCm = Eigen::Vector3d(position[0] / kMmPerCm, position[1] / kMmPerCm, position[2] / kMmPerCm);
  outer.innermostHitInEndcap = (m_finder.regionOf(innermostHit.cellID) != 0);
  const std::optional<edm4hep::TrackState> atIP = CTF::solenoidCrossingState(outer, m_finder.m_solenoid);
  if (atIP)
    track.addToTrackStates(*atIP);
  m_finder.m_solenoidCrossingBuilt += atIP.has_value();
}

bool CombinatorialTrackFinder::EventSearch::allHitsOutsideSolenoid(const TrackCandidate& candidate) const {
  for (size_t hitIndex : candidate.hits) {
    const auto& position = m_hits[hitIndex].positionMm;
    const Eigen::Vector3d positionCm(position[0] / kMmPerCm, position[1] / kMmPerCm, position[2] / kMmPerCm);
    if (m_finder.m_solenoid.contains(positionCm))
      return false;
  }
  return true;
}

// a copy of the input hit in the output, attached to the track, with all its reco-sim links (a hit
// can have several SimHits, e.g. a muon and a delta ray in one cell)
void CombinatorialTrackFinder::EventSearch::copyHitToOutput(const edm4hep::TrackerHitPlane& hit,
                                                            edm4hep::MutableTrack& track) {
  auto outputHit = hit.clone(); // every data member (TrackerHitPlane has no relations)
  m_outputHits.push_back(outputHit);
  const auto simHits = m_recoToSimHits.find(hitKey(hit));
  if (simHits != m_recoToSimHits.end()) {
    for (const auto& simHit : simHits->second) {
      auto link = m_outputLinks.create();
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
    return CTF::kEndcapKeyOffset + key;
  return -CTF::kEndcapKeyOffset - key;
}

// region of a layer key: from the readout's region field (collected in initialize()), else from
// the key itself
int CombinatorialTrackFinder::regionOfLayerKey(int layerKey) const {
  const auto entry = m_regionOfLayerKey.find(layerKey);
  if (entry != m_regionOfLayerKey.end())
    return entry->second;
  return CTF::regionFromKey(layerKey);
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
    position = dd4hep::Position(0.5 * m_solenoid.radiusMm / kMmPerCm, 0.0, distanceCm);
  return m_field.magneticField(position).z() / dd4hep::tesla;
}

bool CombinatorialTrackFinder::outsideCentralField(double distanceMm, bool alongZ, double bzCentral) const {
  return std::abs(fieldAlongScan(distanceMm, alongZ) - bzCentral) > 0.5 * std::abs(bzCentral);
}
