/**
 * SeedFinder: the SeedAndExtend search of CombinatorialTrackFinder (overview in SeedFinder.h).
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

#include "SeedFinder.h"
#include "HelixExtrapolator.h"
#include "TrackingMath.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <unordered_map>

namespace k4rt {

namespace {

  // ---- fixed constants of the search ----
  constexpr double kRaySweepFirstDegree = 0.5; // seed windows: rays from 0.5 to 180 degrees ...
  constexpr double kRaySweepLastDegree = 180.0;
  constexpr double kRaySweepStepDegree = 0.25;       // ... in 0.25 degree steps
  constexpr double kMinimumHitErrorCm = 1e-4;        // smallest hit error in the kink chi2 (no du / dv)
  constexpr double kPrecheckRelativeMargin = 1.01;   // 3-hit pre-check: z cut + 1% ...
  constexpr double kPrecheckAbsoluteMarginCm = 1e-6; // ... + 1e-6 cm (rounding)
  constexpr int kMaxExtensionSteps = 64;             // safety limit, far more layers than any tracker has
  constexpr int kMaxHoleSearchSteps = 12;            // hole check: layers followed outward beyond the last hit
  constexpr double kFillChi2NdfMargin = 1.0;         // fill: new kink chi2/ndf below the current one + this
  constexpr double kSwapMinChi2Gain = 2.0;           // swap: a hit must lower the kink chi2 by more than this
  // shadow hits: layers crossed from 1 cm before the first to 1 cm after the last hit, within 1 cm
  // of the layer edges; hits more than 5 cm from the crossing are not looked at
  constexpr double kShadowArcMarginCm = 1.0;
  constexpr double kShadowLayerToleranceCm = 1.0;
  constexpr double kShadowSearchRadiusCm = 5.0;
  constexpr double kMergeWindowFactor = 3.0; // merging: on the other helix within 3x the search window
  // join of two pieces: its triplets may exceed the triplet cut 4x (skipped layers' material is missing)
  constexpr double kJoinTripletCutFactor = 4.0;

  /// The full fit of an ordered hit list.
  struct TrackFit {
    CircleFit circle;
    Helix helix;
    std::vector<double> arcLengths; ///< of every hit (cm)
    KinkResult kink;
  };

  struct Seed {
    std::vector<size_t> hits;
    double kinkChi2NDF; ///< the lower, the earlier it is extended
  };

  /// A hit as the search sees it (computed once per event): kink-chi2 view and sensor plane.
  struct CachedHit {
    KinkPoint kinkPoint;
    bool hasPlane = false;
    Eigen::Vector3d planeOrigin = Eigen::Vector3d::Zero(); ///< cm
    Eigen::Vector3d planeNormal = Eigen::Vector3d::Zero(); ///< unit vector
  };

  /// Where the hits of one accepted piece lie relative to the helix of another.
  enum class Placement { NotNear, Overlapping, Before, After };

  // helpers, defined at the end of the file
  Eigen::Vector3d positionCm(const edm4hep::TrackerHitPlane& hit);
  void completeSubsets(const std::vector<int>& layers, size_t first, size_t end, int size, std::vector<int>& current,
                       std::set<std::vector<int>>& out);
  void windowSubsets(const std::vector<int>& layers, size_t first, size_t length, int size,
                     std::set<std::vector<int>>& out);
  bool circumcircle(const Eigen::Vector3d& pointA, const Eigen::Vector3d& pointB, const Eigen::Vector3d& pointC,
                    double& centreX, double& centreY, double& radius);
  bool tripletPrecheck(const Eigen::Vector3d& pointA, const Eigen::Vector3d& pointB, const Eigen::Vector3d& pointC,
                       double maxZResidual);
  bool hitOnTrajectory(const Eigen::Vector3d& pointA, const Eigen::Vector3d& pointB, const Eigen::Vector3d& pointC,
                       const Eigen::Vector3d& next, double window, double slope);

} // namespace

SeedFinder::SeedFinder(const LayerGeometry& geometry, SeedConfig config, PlaneFn plane)
    : m_geometry(geometry), m_config(std::move(config)), m_plane(std::move(plane)) {}

std::vector<std::string> SeedFinder::prepare() {
  std::vector<std::string> problems;
  m_setsBySize.clear();
  m_passes.clear();
  for (int seedSize : m_config.passes) {
    if (seedSize < 3)
      problems.push_back("a seed size < 3 was raised to 3");
    m_passes.push_back(std::max(3, seedSize));
  }
  for (int seedSize : m_passes) {
    if (!m_setsBySize.count(seedSize))
      m_setsBySize[seedSize] = buildSets(seedSize);
  }
  // inside-out: windows ordered by the closest approach to the IP of their layers, in turn
  for (auto& [seedSize, sets] : m_setsBySize) {
    std::stable_sort(sets.begin(), sets.end(), [this](const std::vector<int>& left, const std::vector<int>& right) {
      for (size_t i = 0; i < std::min(left.size(), right.size()); ++i) {
        if (closestApproachCm(left[i]) != closestApproachCm(right[i]))
          return closestApproachCm(left[i]) < closestApproachCm(right[i]);
      }
      return left.size() < right.size();
    });
  }
  return problems;
}

std::vector<std::vector<int>> SeedFinder::buildSets(int seedSize) const {
  // the order in which a straight ray from the origin crosses the layers is the inside-out order
  // a track in that direction sees them
  std::set<std::vector<int>> sets;
  const size_t windowLength = static_cast<size_t>(seedSize + std::max(0, m_config.layerSkip));
  for (double degree = kRaySweepFirstDegree; degree < kRaySweepLastDegree; degree += kRaySweepStepDegree) {
    const std::vector<int> crossed = m_geometry.layersCrossedByRay(degree * M_PI / 180.0);
    if (static_cast<int>(crossed.size()) < seedSize)
      continue;
    if (m_config.policy == SeedLayerPolicy::Innermost) {
      windowSubsets(crossed, 0, windowLength, seedSize, sets);
    } else if (m_config.policy == SeedLayerPolicy::Outermost) {
      // innermost on the reversed list, turned back
      const std::vector<int> reversed(crossed.rbegin(), crossed.rend());
      std::set<std::vector<int>> reversedSets;
      windowSubsets(reversed, 0, windowLength, seedSize, reversedSets);
      for (std::vector<int> set : reversedSets) {
        std::reverse(set.begin(), set.end());
        sets.insert(set);
      }
    } else {
      for (size_t first = 0; first + seedSize <= crossed.size(); ++first) // a window at every layer
        windowSubsets(crossed, first, windowLength, seedSize, sets);
    }
  }
  return std::vector<std::vector<int>>(sets.begin(), sets.end());
}

double SeedFinder::closestApproachCm(int layerKey) const {
  const LayerDescriptor* layer = m_geometry.find(layerKey);
  return layer ? layer->closestApproachCm : 1e9;
}

// The search in one event (best seeds first):
//   (a) seed windows inside-out; seeds best kink chi2 first; an accepted track claims its hits at once
//   (b) hits are attached and tracks judged by the kink chi2 with multiple scattering
//   (c) confirmation: fill crossed empty layers, swap in better hits, drop the worst kink; until stable
//   (d) sequential pieces of one trajectory are merged; overlapping arcs are flagged loopers
//   (e) passes of 4-hit, then 3-hit seeds; the last pass claims tracks of >= MinTrackHits
class SeedFinder::EventSearch {
public:
  EventSearch(const SeedFinder& finder, const edm4hep::TrackerHitPlaneCollection& hits,
              const std::vector<std::pair<size_t, int>>& pool);

  void searchPasses();                   ///< (a) + (e)
  void mergePieces();                    ///< (d)
  std::vector<TrackCandidate> results(); ///< the candidates not merged into another, longest first

private:
  // ---- hits ----
  const CachedHit& cachedHit(size_t hitIndex);
  const KinkPoint& kinkPoint(size_t hitIndex) { return cachedHit(hitIndex).kinkPoint; }
  bool anyClaimed(const std::vector<size_t>& hitIndices) const;
  void collectUnclaimedHits();
  bool fitTrack(const std::vector<size_t>& track, TrackFit& fit);
  // ---- (a) seeds ----
  std::vector<Seed> seedsOfWindow(const std::vector<int>& layerSet, int seedSize);
  void growSeeds(int level, const std::vector<int>& layerSet, int seedSize, std::vector<size_t>& chain,
                 std::vector<Seed>& seeds);
  bool chainStillOnTrack(int seedSize, const std::vector<size_t>& chain);
  void keepSeedIfGood(const std::vector<size_t>& chain, int seedSize, std::vector<Seed>& seeds);
  // ---- (b) extension ----
  const LayerDescriptor* nextCrossedLayer(const Helix& helix, double fromArcLength, int direction, double tolerance,
                                          double maxArcLength, const std::vector<int>& visited, Crossing& crossing,
                                          double& distance) const;
  std::vector<std::pair<double, size_t>> nearHits(const Helix& helix, const Crossing& crossing,
                                                  const LayerDescriptor& layer, double window,
                                                  const std::vector<size_t>& track);
  void extend(std::vector<size_t>& track, int direction, std::vector<int>& visited, TrackFit& fit);
  bool attachBestHit(std::vector<size_t>& track, int direction, TrackFit& fit, const Crossing& crossing,
                     const LayerDescriptor& layer, double window);
  // ---- (c) confirmation ----
  void confirm(std::vector<size_t>& track, TrackFit& fit);
  bool fillMissingLayers(std::vector<size_t>& track, TrackFit& fit);
  bool swapBetterHits(std::vector<size_t>& track, TrackFit& fit);
  bool dropWorstKink(std::vector<size_t>& track, TrackFit& fit);
  // ---- acceptance ----
  bool passesTrackCuts(const std::vector<size_t>& track, const TrackFit& fit, int claimMinHits);
  bool hitNear(const LayerDescriptor& layer, const Crossing& crossing);
  int countHoles(const std::vector<size_t>& track, const TrackFit& fit);
  void accept(const std::vector<size_t>& track, const TrackFit& fit);
  void claimShadowHits(const TrackFit& fit);
  // ---- (d) merging ----
  bool mergeOnePair();
  Placement placeOnHelix(size_t longer, size_t shorter);
  bool tryMerge(size_t longer, size_t shorter, bool shorterComesFirst);

  const SeedFinder& m_finder;
  const SeedConfig& m_config;
  const LayerGeometry& m_geometry;
  const edm4hep::TrackerHitPlaneCollection& m_hits;
  const std::vector<std::pair<size_t, int>>& m_pool;
  const KinkChi2Calculator m_kinkCalculator;

  // per hit of the collection (index = hit index)
  std::vector<int> m_layerKeyOfHit;    ///< 0 for hits not in the pool
  std::vector<CachedHit> m_cachedHits; ///< filled on first use
  std::vector<bool> m_isCached;
  std::vector<bool> m_claimed; ///< used by an accepted track

  std::map<int, std::vector<size_t>> m_unclaimedByLayer;         ///< unclaimed pool hits per layer
  std::unordered_map<int, std::vector<size_t>> m_allHitsByLayer; ///< all pool hits per layer
  std::vector<TrackCandidate> m_accepted;                        ///< accepted candidates,
  std::vector<TrackFit> m_acceptedFits;                          ///< their fits,
  std::vector<bool> m_merged;                                    ///< and whether merged into another
};

std::vector<TrackCandidate> SeedFinder::run(const edm4hep::TrackerHitPlaneCollection& hits,
                                            const std::vector<std::pair<size_t, int>>& pool) const {
  EventSearch search(*this, hits, pool);
  search.searchPasses();
  search.mergePieces();
  return search.results();
}

SeedFinder::EventSearch::EventSearch(const SeedFinder& finder, const edm4hep::TrackerHitPlaneCollection& hits,
                                     const std::vector<std::pair<size_t, int>>& pool)
    : m_finder(finder), m_config(finder.m_config), m_geometry(finder.m_geometry), m_hits(hits), m_pool(pool),
      m_kinkCalculator(finder.m_config.common.bFieldTesla), m_layerKeyOfHit(hits.size(), 0), m_cachedHits(hits.size()),
      m_isCached(hits.size(), false), m_claimed(hits.size(), false) {
  // every pool hit by layer, claimed or not (a hole is a crossing with no hit at all nearby)
  for (const auto& [hitIndex, layerKey] : pool) {
    m_layerKeyOfHit[hitIndex] = layerKey;
    m_allHitsByLayer[layerKey].push_back(hitIndex);
  }
}

void SeedFinder::EventSearch::searchPasses() {
  const std::vector<int>& passes = m_finder.m_passes;
  for (size_t passIndex = 0; passIndex < passes.size(); ++passIndex) {
    const int seedSize = passes[passIndex];
    // the last pass claims tracks of >= MinTrackHits, the others only longer ones
    int claimMinHits = std::max(m_config.claimMinHits, m_config.common.minTrackHits);
    if (passIndex + 1 == passes.size())
      claimMinHits = std::max(3, m_config.common.minTrackHits);
    collectUnclaimedHits();
    for (const std::vector<int>& layerSet : m_finder.m_setsBySize.at(seedSize)) {
      bool allLayersHaveHits = true; // only windows with unclaimed hits on every layer
      for (int layerKey : layerSet) {
        if (!m_unclaimedByLayer.count(layerKey)) {
          allLayersHaveHits = false;
          break;
        }
      }
      if (!allLayersHaveHits)
        continue;

      for (const Seed& seed : seedsOfWindow(layerSet, seedSize)) {
        if (anyClaimed(seed.hits))
          continue; // claimed by a track accepted after this seed was made
        std::vector<size_t> track = seed.hits;
        TrackFit fit;
        if (!fitTrack(track, fit))
          continue;
        std::vector<int> visited(layerSet.begin(), layerSet.end());
        // innermost seeds grow outward, outermost inward, consecutive both ways
        if (m_config.policy != SeedLayerPolicy::Outermost)
          extend(track, +1, visited, fit);
        if (m_config.policy != SeedLayerPolicy::Innermost)
          extend(track, -1, visited, fit);
        confirm(track, fit);
        if (passesTrackCuts(track, fit, claimMinHits))
          accept(track, fit);
      }
    }
  }
}

// ============================================================================================
// (a) seeds
// ============================================================================================

std::vector<Seed> SeedFinder::EventSearch::seedsOfWindow(const std::vector<int>& layerSet, int seedSize) {
  std::vector<Seed> seeds;
  std::vector<size_t> chain;
  growSeeds(0, layerSet, seedSize, chain, seeds);
  std::stable_sort(seeds.begin(), seeds.end(),
                   [](const Seed& left, const Seed& right) { return left.kinkChi2NDF < right.kinkChi2NDF; });
  return seeds;
}

// every chain of one unclaimed hit per window layer (recursive, one layer per call), with cheap
// checks before the fit; complete chains go to keepSeedIfGood
void SeedFinder::EventSearch::growSeeds(int level, const std::vector<int>& layerSet, int seedSize,
                                        std::vector<size_t>& chain, std::vector<Seed>& seeds) {
  if (level == seedSize) {
    keepSeedIfGood(chain, seedSize, seeds);
    return;
  }
  for (size_t hitIndex : m_unclaimedByLayer.at(layerSet[level])) {
    if (m_claimed[hitIndex])
      continue;
    // consecutive seed hits at most MaxPairHitDistance apart
    if (level > 0 &&
        (kinkPoint(hitIndex).position - kinkPoint(chain.back()).position).norm() > m_config.common.maxPairHitDistanceCm)
      continue;
    chain.push_back(hitIndex);
    if (level == 0 || chainStillOnTrack(seedSize, chain))
      growSeeds(level + 1, layerSet, seedSize, chain, seeds);
    chain.pop_back();
  }
}

// seeds of 4+ hits: the first three pass the 3-hit pre-check (twice the z cut), each later hit is
// on the helix through the first, a middle and the previous hit
bool SeedFinder::EventSearch::chainStillOnTrack(int seedSize, const std::vector<size_t>& chain) {
  if (seedSize <= 3)
    return true;
  if (chain.size() == 3)
    return tripletPrecheck(kinkPoint(chain[0]).position, kinkPoint(chain[1]).position, kinkPoint(chain[2]).position,
                           2.0 * m_config.maxZResid_cm);
  if (chain.size() >= 4) {
    const size_t last = chain.size() - 1;
    return hitOnTrajectory(kinkPoint(chain[0]).position, kinkPoint(chain[(last - 1) / 2]).position,
                           kinkPoint(chain[last - 1]).position, kinkPoint(chain[last]).position, m_config.window_cm,
                           m_config.windowSlope);
  }
  return true;
}

// a complete chain is a seed if its hits are unclaimed, it passes the pre-check (3 hits), its fit
// works, its hits follow the travel order and its kink chi2/ndf is small
void SeedFinder::EventSearch::keepSeedIfGood(const std::vector<size_t>& chain, int seedSize, std::vector<Seed>& seeds) {
  if (anyClaimed(chain))
    return;
  if (seedSize == 3 && !tripletPrecheck(kinkPoint(chain[0]).position, kinkPoint(chain[1]).position,
                                        kinkPoint(chain[2]).position, m_config.maxZResid_cm))
    return;
  TrackFit fit;
  if (!fitTrack(chain, fit))
    return;
  for (size_t hit = 1; hit < fit.arcLengths.size(); ++hit) {
    if (!(fit.arcLengths[hit] > fit.arcLengths[hit - 1]))
      return;
  }
  if (fit.kink.chi2ndf() > m_config.seedKinkMaxChi2NDF)
    return;
  seeds.push_back({chain, fit.kink.chi2ndf()});
}

// ============================================================================================
// (b) extension
// ============================================================================================

// layer by layer along the helix, the best hit of each crossed layer is attached; stops after more
// than maxMissedLayers crossed layers in a row without a hit
void SeedFinder::EventSearch::extend(std::vector<size_t>& track, int direction, std::vector<int>& visited,
                                     TrackFit& fit) {
  int missed = 0;
  for (int step = 0; step < kMaxExtensionSteps; ++step) {
    const double fromArcLength = direction > 0 ? fit.arcLengths.back() : fit.arcLengths.front();
    // the next layer within half a turn (pi R): beyond it the track curls back
    Crossing crossing;
    double distance = 0.0;
    const LayerDescriptor* layer = nextCrossedLayer(fit.helix, fromArcLength, direction, m_config.window_cm,
                                                    M_PI * fit.circle.radius, visited, crossing, distance);
    if (!layer)
      break;
    visited.push_back(layer->layerKey);
    // the window grows with the extrapolation distance (the prediction gets less precise)
    const double window = m_config.window_cm + m_config.windowSlope * distance;
    if (attachBestHit(track, direction, fit, crossing, *layer, window))
      missed = 0;
    else if (++missed > m_config.maxMissedLayers)
      break;
  }
}

// the next unvisited layer the helix crosses from fromArcLength in the direction (within
// maxArcLength), with the crossing and its distance; nullptr if none
const LayerDescriptor* SeedFinder::EventSearch::nextCrossedLayer(const Helix& helix, double fromArcLength,
                                                                 int direction, double tolerance, double maxArcLength,
                                                                 const std::vector<int>& visited, Crossing& crossing,
                                                                 double& distance) const {
  const LayerDescriptor* next = nullptr;
  distance = std::numeric_limits<double>::infinity();
  for (const auto& layer : m_geometry.layers()) {
    if (std::find(visited.begin(), visited.end(), layer.layerKey) != visited.end())
      continue;
    const Crossing candidate = helix.crossLayer(layer, fromArcLength, direction, tolerance, maxArcLength);
    if (!candidate.ok)
      continue;
    const double candidateDistance = direction * (candidate.arcLength - fromArcLength);
    if (candidateDistance < distance) {
      distance = candidateDistance;
      crossing = candidate;
      next = &layer;
    }
  }
  return next;
}

// The hit with the smallest chi2 increase of the refit that includes it (below the triplet cut).
// Judging only the kink against the current curvature would fail across long gaps between layers.
bool SeedFinder::EventSearch::attachBestHit(std::vector<size_t>& track, int direction, TrackFit& fit,
                                            const Crossing& crossing, const LayerDescriptor& layer, double window) {
  std::vector<size_t> bestTrack;
  TrackFit bestFit;
  double bestChi2Increase = std::numeric_limits<double>::infinity();
  const double sigmaRadius = std::sqrt(std::max(0.0, fit.circle.covariance(2, 2)));
  for (const auto& [distance, hitIndex] : nearHits(fit.helix, crossing, layer, window, track)) {
    // first the kink against the current curvature, with its error (strict where it is well measured)
    const size_t size = track.size();
    double tripletChi2 = 0.0;
    if (direction > 0)
      tripletChi2 = m_kinkCalculator.triplet(kinkPoint(track[size - 2]), kinkPoint(track[size - 1]),
                                             kinkPoint(hitIndex), fit.circle.radius, fit.helix.rotation, sigmaRadius);
    else
      tripletChi2 = m_kinkCalculator.triplet(kinkPoint(hitIndex), kinkPoint(track[0]), kinkPoint(track[1]),
                                             fit.circle.radius, fit.helix.rotation, sigmaRadius);
    if (tripletChi2 > m_config.kinkMaxChi2Triplet)
      continue;
    // then the chi2 increase of the refit
    std::vector<size_t> trial = track;
    if (direction > 0)
      trial.push_back(hitIndex);
    else
      trial.insert(trial.begin(), hitIndex);
    TrackFit trialFit;
    if (!fitTrack(trial, trialFit))
      continue;
    const double chi2Increase = trialFit.kink.chi2 - fit.kink.chi2;
    if (chi2Increase < bestChi2Increase) {
      bestChi2Increase = chi2Increase;
      bestTrack = std::move(trial);
      bestFit = std::move(trialFit);
    }
  }
  if (bestTrack.empty() || bestChi2Increase > m_config.kinkMaxChi2Triplet)
    return false;
  track = std::move(bestTrack);
  fit = std::move(bestFit);
  return true;
}

// hits of the layer near the crossing, nearest first (at most MaxCandidatesPerLayer), the distance
// measured on each hit's own sensor plane
std::vector<std::pair<double, size_t>> SeedFinder::EventSearch::nearHits(const Helix& helix, const Crossing& crossing,
                                                                         const LayerDescriptor& layer, double window,
                                                                         const std::vector<size_t>& track) {
  std::vector<std::pair<double, size_t>> near;
  auto layerHits = m_unclaimedByLayer.find(layer.layerKey);
  if (layerHits == m_unclaimedByLayer.end())
    return near;
  for (size_t hitIndex : layerHits->second) {
    if (m_claimed[hitIndex] || std::find(track.begin(), track.end(), hitIndex) != track.end())
      continue;
    const Eigen::Vector3d position = kinkPoint(hitIndex).position;
    if ((position - crossing.position).norm() > window + layer.depthCm())
      continue; // quick test, with room for the layer depth
    double distance = (position - crossing.position).norm();
    const CachedHit& cached = cachedHit(hitIndex);
    if (cached.hasPlane) {
      const Crossing onSensor = helix.refineToPlane(cached.planeOrigin, cached.planeNormal, crossing.arcLength);
      if (onSensor.ok)
        distance = (position - onSensor.position).norm();
    }
    if (distance <= window)
      near.emplace_back(distance, hitIndex);
  }
  std::sort(near.begin(), near.end());
  if (static_cast<int>(near.size()) > m_config.maxCandidatesPerLayer)
    near.resize(m_config.maxCandidatesPerLayer);
  return near;
}

// ============================================================================================
// (c) confirmation
// ============================================================================================

void SeedFinder::EventSearch::confirm(std::vector<size_t>& track, TrackFit& fit) {
  for (int iteration = 0; iteration < m_config.confirmIterations; ++iteration) {
    // all three steps run every iteration, in this order
    const bool filled = fillMissingLayers(track, fit);
    const bool swapped = swapBetterHits(track, fit);
    const bool dropped = dropWorstKink(track, fit);
    if (!filled && !swapped && !dropped)
      break;
  }
}

// a layer the helix crosses between the first and last hit, without a hit of the track, gets the
// hit with the lowest kink chi2/ndf, if no triplet of the refit exceeds the triplet cut
bool SeedFinder::EventSearch::fillMissingLayers(std::vector<size_t>& track, TrackFit& fit) {
  bool changed = false;
  std::set<int> present;
  for (size_t hitIndex : track)
    present.insert(m_layerKeyOfHit[hitIndex]);
  for (const auto& layer : m_geometry.layers()) {
    if (present.count(layer.layerKey))
      continue;
    const Crossing crossing = fit.helix.crossLayer(layer, fit.arcLengths.front(), +1, m_config.window_cm,
                                                   fit.arcLengths.back() - fit.arcLengths.front());
    if (!crossing.ok)
      continue;
    std::vector<size_t> bestTrack;
    TrackFit bestFit;
    double bestChi2NDF = fit.kink.chi2ndf() + kFillChi2NdfMargin;
    for (const auto& [distance, hitIndex] : nearHits(fit.helix, crossing, layer, m_config.window_cm, track)) {
      // the hit goes where its arc length puts it (copy loop, not vector::insert: gcc 16 false warning)
      const double arcLength = fit.helix.arcLengthNear(kinkPoint(hitIndex).position, crossing.arcLength);
      size_t insertAt = 0;
      while (insertAt < fit.arcLengths.size() && fit.arcLengths[insertAt] < arcLength)
        ++insertAt;
      std::vector<size_t> trial;
      trial.reserve(track.size() + 1);
      for (size_t i = 0; i < track.size(); ++i) {
        if (i == insertAt)
          trial.push_back(hitIndex);
        trial.push_back(track[i]);
      }
      if (insertAt == track.size())
        trial.push_back(hitIndex);
      TrackFit trialFit;
      if (!fitTrack(trial, trialFit) || !(trialFit.kink.chi2ndf() < bestChi2NDF))
        continue;
      double largestTripletChi2 = 0;
      for (double tripletChi2 : trialFit.kink.perTriplet)
        largestTripletChi2 = std::max(largestTripletChi2, tripletChi2);
      if (largestTripletChi2 <= m_config.kinkMaxChi2Triplet) {
        bestChi2NDF = trialFit.kink.chi2ndf();
        bestTrack = std::move(trial);
        bestFit = std::move(trialFit);
      }
    }
    if (!bestTrack.empty()) {
      track = std::move(bestTrack);
      fit = std::move(bestFit);
      present.insert(layer.layerKey);
      changed = true;
    }
  }
  return changed;
}

// a different hit of the same layer (within the window) that lowers the kink chi2 clearly
bool SeedFinder::EventSearch::swapBetterHits(std::vector<size_t>& track, TrackFit& fit) {
  bool changed = false;
  for (size_t i = 0; i < track.size(); ++i) {
    const int layerKey = m_layerKeyOfHit[track[i]];
    auto layerHits = m_unclaimedByLayer.find(layerKey);
    if (!m_geometry.find(layerKey) || layerHits == m_unclaimedByLayer.end())
      continue;
    for (size_t hitIndex : layerHits->second) {
      if (hitIndex == track[i] || m_claimed[hitIndex] || std::find(track.begin(), track.end(), hitIndex) != track.end())
        continue;
      if ((kinkPoint(hitIndex).position - kinkPoint(track[i]).position).norm() > m_config.window_cm)
        continue;
      std::vector<size_t> trial = track;
      trial[i] = hitIndex;
      TrackFit trialFit;
      if (fitTrack(trial, trialFit) && trialFit.kink.chi2 < fit.kink.chi2 - kSwapMinChi2Gain) {
        track = std::move(trial);
        fit = std::move(trialFit);
        changed = true;
      }
    }
  }
  return changed;
}

// worst triplet beyond the cut: remove the one of its three hits giving the lowest kink chi2/ndf
// (if lower than now)
bool SeedFinder::EventSearch::dropWorstKink(std::vector<size_t>& track, TrackFit& fit) {
  if (track.size() <= 3)
    return false;
  int worst = -1;
  double worstChi2 = 0;
  for (size_t i = 0; i < fit.kink.perTriplet.size(); ++i) {
    if (fit.kink.perTriplet[i] > worstChi2) {
      worstChi2 = fit.kink.perTriplet[i];
      worst = static_cast<int>(i);
    }
  }
  if (worst < 0 || worstChi2 <= m_config.kinkMaxChi2Triplet)
    return false;
  std::vector<size_t> bestTrack;
  TrackFit bestFit;
  double bestChi2NDF = fit.kink.chi2ndf();
  for (int removed = worst - 1; removed <= worst + 1; ++removed) {
    if (removed < 0 || removed >= static_cast<int>(track.size()))
      continue;
    std::vector<size_t> trial = track;
    trial.erase(trial.begin() + removed);
    TrackFit trialFit;
    if (fitTrack(trial, trialFit) && trialFit.kink.chi2ndf() < bestChi2NDF) {
      bestChi2NDF = trialFit.kink.chi2ndf();
      bestTrack = std::move(trial);
      bestFit = std::move(trialFit);
    }
  }
  if (bestTrack.empty())
    return false;
  track = std::move(bestTrack);
  fit = std::move(bestFit);
  return true;
}

// ============================================================================================
// acceptance
// ============================================================================================

// enough hits and a small kink chi2/ndf; a short track also needs MinComboPT and few holes
bool SeedFinder::EventSearch::passesTrackCuts(const std::vector<size_t>& track, const TrackFit& fit, int claimMinHits) {
  if (static_cast<int>(track.size()) < claimMinHits || fit.kink.chi2ndf() > m_config.trackKinkMaxChi2NDF)
    return false;
  if (static_cast<int>(track.size()) > m_config.shortTrackMaxHits)
    return true;
  if (m_config.common.minPtGeV > 0.0 &&
      ptFromRadius(m_config.common.bFieldTesla, fit.circle.radius) < m_config.common.minPtGeV)
    return false;
  return countHoles(track, fit) <= m_config.shortTrackMaxHoles;
}

// holes inside the track's span, then outward until it leaves the tracker (or half a turn); the
// crossings must be inside the layer (no tolerance), so passing between layers is no hole
int SeedFinder::EventSearch::countHoles(const std::vector<size_t>& track, const TrackFit& fit) {
  int holes = 0;
  std::vector<int> visited;
  for (size_t hitIndex : track)
    visited.push_back(m_layerKeyOfHit[hitIndex]);
  const double firstArcLength = fit.arcLengths.front();
  const double lastArcLength = fit.arcLengths.back();
  for (const auto& layer : m_geometry.layers()) { // inside the span
    if (std::find(visited.begin(), visited.end(), layer.layerKey) != visited.end())
      continue;
    const Crossing crossing = fit.helix.crossLayer(layer, firstArcLength, +1, -1.0, lastArcLength - firstArcLength);
    if (crossing.ok && !hitNear(layer, crossing))
      ++holes;
  }
  double fromArcLength = lastArcLength; // outward, layer by layer
  for (int step = 0; step < kMaxHoleSearchSteps && holes <= m_config.shortTrackMaxHoles; ++step) {
    Crossing crossing;
    double distance = 0.0;
    const double maxArcLength = M_PI * fit.circle.radius - (fromArcLength - firstArcLength);
    const LayerDescriptor* layer =
        nextCrossedLayer(fit.helix, fromArcLength, +1, -1.0, maxArcLength, visited, crossing, distance);
    if (!layer)
      break;
    visited.push_back(layer->layerKey);
    if (!hitNear(*layer, crossing))
      ++holes;
    fromArcLength = crossing.arcLength;
  }
  return holes;
}

// any pool hit (claimed or not) near the crossing of the layer?
bool SeedFinder::EventSearch::hitNear(const LayerDescriptor& layer, const Crossing& crossing) {
  auto layerHits = m_allHitsByLayer.find(layer.layerKey);
  if (layerHits == m_allHitsByLayer.end())
    return false;
  for (size_t hitIndex : layerHits->second) {
    if ((kinkPoint(hitIndex).position - crossing.position).norm() <= m_config.window_cm + layer.depthCm())
      return true;
  }
  return false;
}

// an accepted track is stored; its hits (and overlapping-sensor double hits) are claimed
void SeedFinder::EventSearch::accept(const std::vector<size_t>& track, const TrackFit& fit) {
  TrackCandidate candidate;
  candidate.hits = track;
  candidate.circle = fit.circle;
  candidate.kinkChi2 = fit.kink.chi2;
  candidate.kinkNdf = fit.kink.ndf;
  candidate.looper = (fit.arcLengths.back() - fit.arcLengths.front()) > M_PI * fit.circle.radius;
  for (size_t hitIndex : track)
    m_claimed[hitIndex] = true;
  if (m_config.shadowHitWindow_cm > 0.0)
    claimShadowHits(fit);
  m_accepted.push_back(std::move(candidate));
  m_acceptedFits.push_back(fit);
  m_merged.push_back(false);
}

// unclaimed hits within ShadowHitWindow of the trajectory on their own sensor (overlapping-sensor
// double hits) are claimed, so they cannot seed a duplicate
void SeedFinder::EventSearch::claimShadowHits(const TrackFit& fit) {
  const double firstArcLength = fit.arcLengths.front();
  const double lastArcLength = fit.arcLengths.back();
  for (const auto& layer : m_geometry.layers()) {
    auto layerHits = m_allHitsByLayer.find(layer.layerKey);
    if (layerHits == m_allHitsByLayer.end())
      continue;
    const Crossing crossing =
        fit.helix.crossLayer(layer, firstArcLength - kShadowArcMarginCm, +1, kShadowLayerToleranceCm,
                             lastArcLength - firstArcLength + 2.0 * kShadowArcMarginCm);
    if (!crossing.ok)
      continue;
    for (size_t hitIndex : layerHits->second) {
      if (m_claimed[hitIndex])
        continue;
      const Eigen::Vector3d position = kinkPoint(hitIndex).position;
      if ((position - crossing.position).norm() > kShadowSearchRadiusCm)
        continue;
      const CachedHit& cached = cachedHit(hitIndex);
      if (!cached.hasPlane)
        continue;
      const Crossing onSensor = fit.helix.refineToPlane(cached.planeOrigin, cached.planeNormal, crossing.arcLength);
      if (onSensor.ok && (position - onSensor.position).norm() <= m_config.shadowHitWindow_cm)
        m_claimed[hitIndex] = true;
    }
  }
}

// ============================================================================================
// (d) merging
// ============================================================================================

void SeedFinder::EventSearch::mergePieces() {
  while (mergeOnePair()) { // until no pair can be merged
  }
}

// merges the first pair (longer piece, shorter piece) that can be merged; true if one was.
// Pieces overlapping on the same helix are flagged loopers.
bool SeedFinder::EventSearch::mergeOnePair() {
  for (size_t longer = 0; longer < m_accepted.size(); ++longer) {
    if (m_merged[longer])
      continue;
    for (size_t shorter = 0; shorter < m_accepted.size(); ++shorter) {
      if (shorter == longer || m_merged[shorter] || m_accepted[shorter].hits.size() > m_accepted[longer].hits.size())
        continue;
      const Placement placement = placeOnHelix(longer, shorter);
      if (placement == Placement::NotNear)
        continue;
      if (placement == Placement::Overlapping) {
        m_accepted[longer].looper = true;
        m_accepted[shorter].looper = true;
        continue;
      }
      if (tryMerge(longer, shorter, placement == Placement::Before))
        return true;
    }
  }
  return false;
}

// are all hits of the shorter piece close to the helix of the longer one, all before or all after it?
Placement SeedFinder::EventSearch::placeOnHelix(size_t longer, size_t shorter) {
  const TrackFit& fit = m_acceptedFits[longer];
  const double firstArcLength = fit.arcLengths.front();
  const double lastArcLength = fit.arcLengths.back();
  bool allBefore = true;
  bool allAfter = true;
  for (size_t hitIndex : m_accepted[shorter].hits) {
    const Eigen::Vector3d position = kinkPoint(hitIndex).position;
    // placed on the helix from the nearer end of the longer piece
    const bool nearerToLast = (position - fit.helix.positionAt(lastArcLength)).norm() <
                              (position - fit.helix.positionAt(firstArcLength)).norm();
    const double arcLength = fit.helix.arcLengthNear(position, nearerToLast ? lastArcLength : firstArcLength);
    const double extrapolation = std::max(0.0, std::max(arcLength - lastArcLength, firstArcLength - arcLength));
    const double window = kMergeWindowFactor * (m_config.window_cm + m_config.windowSlope * extrapolation);
    if ((position - fit.helix.positionAt(arcLength)).norm() > window || extrapolation > M_PI * fit.circle.radius)
      return Placement::NotNear;
    allBefore = allBefore && arcLength < firstArcLength;
    allAfter = allAfter && arcLength > lastArcLength;
  }
  if (!allBefore && !allAfter)
    return Placement::Overlapping;
  return allBefore ? Placement::Before : Placement::After;
}

// the pieces joined in travel order and refitted; kept below the track cut, with no triplet beyond
// the triplet cut (4x at the join)
bool SeedFinder::EventSearch::tryMerge(size_t longer, size_t shorter, bool shorterComesFirst) {
  const std::vector<size_t>& firstPiece = shorterComesFirst ? m_accepted[shorter].hits : m_accepted[longer].hits;
  const std::vector<size_t>& secondPiece = shorterComesFirst ? m_accepted[longer].hits : m_accepted[shorter].hits;
  std::vector<size_t> joined = firstPiece;
  joined.insert(joined.end(), secondPiece.begin(), secondPiece.end());
  TrackFit joinedFit;
  if (!fitTrack(joined, joinedFit) || joinedFit.kink.chi2ndf() > m_config.trackKinkMaxChi2NDF)
    return false;
  const size_t join = firstPiece.size(); // first hit of the second piece
  for (size_t triplet = 1; triplet + 1 < joined.size(); ++triplet) {
    const bool atJoin = (triplet + 1 == join) || (triplet == join);
    const double cut = (atJoin ? kJoinTripletCutFactor : 1.0) * m_config.kinkMaxChi2Triplet;
    if (joinedFit.kink.perTriplet[triplet] > cut)
      return false;
  }
  TrackCandidate& merged = m_accepted[longer];
  merged.hits = joined;
  merged.circle = joinedFit.circle;
  merged.kinkChi2 = joinedFit.kink.chi2;
  merged.kinkNdf = joinedFit.kink.ndf;
  merged.looper =
      merged.looper || (joinedFit.arcLengths.back() - joinedFit.arcLengths.front()) > M_PI * joinedFit.circle.radius;
  m_acceptedFits[longer] = std::move(joinedFit);
  m_merged[shorter] = true;
  return true;
}

std::vector<TrackCandidate> SeedFinder::EventSearch::results() {
  std::vector<TrackCandidate> candidates;
  for (size_t i = 0; i < m_accepted.size(); ++i) {
    if (!m_merged[i])
      candidates.push_back(std::move(m_accepted[i]));
  }
  std::stable_sort(candidates.begin(), candidates.end(), [](const TrackCandidate& left, const TrackCandidate& right) {
    return left.hits.size() > right.hits.size();
  });
  return candidates;
}

// ============================================================================================
// hits
// ============================================================================================

const CachedHit& SeedFinder::EventSearch::cachedHit(size_t hitIndex) {
  CachedHit& cached = m_cachedHits[hitIndex];
  if (m_isCached[hitIndex])
    return cached;
  m_isCached[hitIndex] = true;
  const auto hit = m_hits[hitIndex];
  KinkPoint& point = cached.kinkPoint;
  point.position = positionCm(hit);
  point.sigmaTransverseCm = std::max(kMinimumHitErrorCm, kCmPerMm * hit.getDu());
  point.sigmaZCm = std::max(kMinimumHitErrorCm, kCmPerMm * hit.getDv());
  const LayerDescriptor* layer = m_geometry.find(m_layerKeyOfHit[hitIndex]);
  point.xOverX0 = layer ? layer->xOverX0 : 0.0;
  cached.hasPlane = m_finder.m_plane(hit.getCellID(), cached.planeOrigin, cached.planeNormal);
  if (cached.hasPlane)
    point.normal = cached.planeNormal;
  else if (layer && !layer->isBarrel())
    point.normal = Eigen::Vector3d::UnitZ(); // a disk faces along z
  else
    point.normal = Eigen::Vector3d(point.position.x(), point.position.y(), 0).normalized(); // barrel: outward
  return cached;
}

bool SeedFinder::EventSearch::anyClaimed(const std::vector<size_t>& hitIndices) const {
  for (size_t hitIndex : hitIndices) {
    if (m_claimed[hitIndex])
      return true;
  }
  return false;
}

void SeedFinder::EventSearch::collectUnclaimedHits() {
  m_unclaimedByLayer.clear();
  for (const auto& [hitIndex, layerKey] : m_pool) {
    if (!m_claimed[hitIndex])
      m_unclaimedByLayer[layerKey].push_back(hitIndex);
  }
}

// circle (fitCircle), helix and kink chi2 of an ordered hit list
bool SeedFinder::EventSearch::fitTrack(const std::vector<size_t>& track, TrackFit& fit) {
  if (track.size() < 3)
    return false;
  std::vector<Eigen::Vector3d> positions;
  std::vector<KinkPoint> kinkPoints;
  for (size_t hitIndex : track) {
    positions.push_back(kinkPoint(hitIndex).position);
    kinkPoints.push_back(kinkPoint(hitIndex));
  }
  const std::optional<CircleFit> circle = fitCircle(m_hits, track, m_config.common.sigmaHitDefaultCm);
  if (!circle)
    return false;
  fit.circle = *circle;
  if (!Helix::fromHits(fit.circle.centreX, fit.circle.centreY, fit.circle.radius, positions, fit.helix,
                       &fit.arcLengths))
    return false;
  fit.kink = m_kinkCalculator.track(kinkPoints, fit.circle.radius, fit.helix.rotation);
  return true;
}

// ============================================================================================
// helper formulas (declared at the top)
// ============================================================================================

namespace {

  Eigen::Vector3d positionCm(const edm4hep::TrackerHitPlane& hit) {
    const auto& position = hit.getPosition();
    return {position[0] * kCmPerMm, position[1] * kCmPerMm, position[2] * kCmPerMm};
  }

  // every way of completing `current` to `size` layers, taking layers in order from layers[first..end)
  // (recursive: one more layer per call)
  void completeSubsets(const std::vector<int>& layers, size_t first, size_t end, int size, std::vector<int>& current,
                       std::set<std::vector<int>>& out) {
    if (static_cast<int>(current.size()) == size) {
      out.insert(current);
      return;
    }
    for (size_t next = first; next < end; ++next) {
      current.push_back(layers[next]);
      completeSubsets(layers, next + 1, end, size, current, out);
      current.pop_back();
    }
  }

  // all subsets of `size` layers of layers[first..first+length) that contain layers[first], in order
  void windowSubsets(const std::vector<int>& layers, size_t first, size_t length, int size,
                     std::set<std::vector<int>>& out) {
    if (first + size > layers.size())
      return;
    const size_t end = std::min(layers.size(), first + length);
    std::vector<int> current{layers[first]};
    completeSubsets(layers, first + 1, end, size, current, out);
  }

  // Circle through three points in the transverse plane: with D = 2 (ax (by - cy) + bx (cy - ay) +
  // cx (ay - by)) and a2 = ax^2 + ay^2 etc., the centre is ((a2 (by - cy) + b2 (cy - ay) +
  // c2 (ay - by)) / D, (a2 (cx - bx) + b2 (ax - cx) + c2 (bx - ax)) / D). false if collinear (D = 0).
  bool circumcircle(const Eigen::Vector3d& pointA, const Eigen::Vector3d& pointB, const Eigen::Vector3d& pointC,
                    double& centreX, double& centreY, double& radius) {
    const double ax = pointA.x(), ay = pointA.y(), bx = pointB.x(), by = pointB.y(), cx = pointC.x(), cy = pointC.y();
    const double twiceArea = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if (std::abs(twiceArea) < 1e-12)
      return false;
    const double squaredA = ax * ax + ay * ay;
    const double squaredB = bx * bx + by * by;
    const double squaredC = cx * cx + cy * cy;
    centreX = (squaredA * (by - cy) + squaredB * (cy - ay) + squaredC * (ay - by)) / twiceArea;
    centreY = (squaredA * (cx - bx) + squaredB * (ax - cx) + squaredC * (bx - ax)) / twiceArea;
    radius = std::hypot(ax - centreX, ay - centreY);
    return true;
  }

  // 3-hit seed checks without a fit: for 3 points the fitted circle is the circumcircle, so the
  // travel order and the z residuals can be tested here (1% margin; the full checks still run).
  // Collinear points pass to the full checks unchanged.
  bool tripletPrecheck(const Eigen::Vector3d& pointA, const Eigen::Vector3d& pointB, const Eigen::Vector3d& pointC,
                       double maxZResidual) {
    double centreX = 0.0, centreY = 0.0, radius = 0.0;
    if (!circumcircle(pointA, pointB, pointC, centreX, centreY, radius))
      return true;
    Helix helix;
    std::vector<double> arcLengths;
    if (!Helix::fromHits(centreX, centreY, radius, {pointA, pointB, pointC}, helix, &arcLengths))
      return true;
    if (!(arcLengths[1] > arcLengths[0]) || !(arcLengths[2] > arcLengths[1]))
      return false; // not in travel order
    const Eigen::Vector3d* points[3] = {&pointA, &pointB, &pointC};
    double largestZResidual = 0.0;
    for (int hit = 0; hit < 3; ++hit) {
      const double zOnHelix = helix.referenceZ + helix.tanLambda * arcLengths[hit];
      largestZResidual = std::max(largestZResidual, std::abs(points[hit]->z() - zOnHelix));
    }
    return largestZResidual <= maxZResidual * kPrecheckRelativeMargin + kPrecheckAbsoluteMarginCm;
  }

  // is `next` on the helix through a, b, c (travel order) within window + slope * (arc from c)?
  // collinear anchors: accepted (left to the seed fit)
  bool hitOnTrajectory(const Eigen::Vector3d& pointA, const Eigen::Vector3d& pointB, const Eigen::Vector3d& pointC,
                       const Eigen::Vector3d& next, double window, double slope) {
    double centreX = 0.0, centreY = 0.0, radius = 0.0;
    if (!circumcircle(pointA, pointB, pointC, centreX, centreY, radius))
      return true;
    Helix helix;
    std::vector<double> arcLengths;
    if (!Helix::fromHits(centreX, centreY, radius, {pointA, pointB, pointC}, helix, &arcLengths))
      return true;
    const double arcLengthOfNext = helix.arcLengthNear(next, arcLengths.back());
    if (!(arcLengthOfNext > arcLengths.back()))
      return false; // behind the chain
    const double tolerance = window + slope * (arcLengthOfNext - arcLengths.back());
    return (next - helix.positionAt(arcLengthOfNext)).norm() <= tolerance;
  }

} // namespace

} // namespace k4rt
