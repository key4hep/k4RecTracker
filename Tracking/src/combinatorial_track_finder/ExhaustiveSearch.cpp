/**
 * ExhaustiveSearch: the Exhaustive strategy of CombinatorialTrackFinder (overview in ExhaustiveSearch.h).
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

#include "ExhaustiveSearch.h"
#include "LayerGeometry.h"
#include "TrackingMath.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>

namespace CTF {

namespace {

  // ---- fixed constants ----
  // road cuts: within 15 degrees of the beam (either side) the phi cuts are skipped (low-pT tracks
  // turn a lot between endcap layers)
  constexpr double kNearBeamThetaMinDeg = 15.0;
  constexpr double kPairTolerancePhi = 0.15;      // two hits of a doublet are coincident within
  constexpr double kPairToleranceEta = 0.10;      // these |dphi| (rad) and |deta|
  constexpr int kMinCompleteDoublets = 2;         // pair coincidence: complete doublets needed
  constexpr double kGuidedGeometricWeight = 0.75; // crowded layers: score = 0.75 * geometric term
  constexpr double kGuidedEdepWeight = 0.25;      //                       + 0.25 * edep term

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

  // helpers, defined at the end of the file
  bool nextLayerChoice(std::vector<int>& chosen, int numberOfLayers);
  long completionsAfter(const std::vector<const std::vector<size_t>*>& hitsOfChosenLayer, size_t depth);

} // namespace

ExhaustiveSearch::ExhaustiveSearch(Config config, FieldFunction bzAt)
    : m_config(std::move(config)), m_bzAt(std::move(bzAt)) {}

// One track per iteration, until no candidate is left: collect the free hits, find the best
// combination, refine it (crowded layers, neighbour tracks, outliers), accept it.
class ExhaustiveSearch::EventSearch {
public:
  EventSearch(const ExhaustiveSearch& search, const std::vector<EventHit>& hits, ExhaustiveStatistics& statistics);

  std::vector<TrackCandidate> run();

private:
  // ---- the free hits of one iteration ----
  void prepareIteration();
  void computeEdepScores();
  void findIsolationNeighbours();

  // ---- the best combination ----
  bool searchCombinations(TrackCandidate& best);
  HitCluster pairCoincidenceCluster(const HitCluster& allHits) const;
  std::vector<HitCluster> buildClusters(const HitCluster& allHits) const;
  void joinLayers(UnionFind& groups, const HitCluster& allHits, const std::vector<size_t>& poolPosition, int step,
                  double phiLimit, double distanceLimit) const;
  void searchCluster(const HitCluster& cluster, const RoadCuts& roadCuts, SearchResult& result) const;
  void chooseHits(size_t depth, const std::vector<const std::vector<size_t>*>& hitsOfChosenLayer,
                  const RoadCuts& roadCuts, std::vector<size_t>& combination, SearchResult& result) const;
  bool passesPreFitCuts(const std::vector<size_t>& combination, const RoadCuts& roadCuts) const;
  void offerCombination(const std::vector<size_t>& combination, SearchResult& result) const;
  double meanEdepScore(const std::vector<size_t>& combination) const;
  bool isBetterCandidate(int numberOfHits, const CircleFit& fit, double edepScoreKeV, const TrackCandidate& best) const;
  static bool isBetterFallback(int numberOfHits, const CircleFit& fit, const TrackCandidate& fallback);

  // ---- refinement of the chosen combination ----
  void selectCrowdedLayerHits(TrackCandidate& best) const;
  std::vector<size_t> cleanLayerHits() const;
  bool addBestHitOfCrowdedLayers(const std::function<double(size_t)>& geometricTerm, std::vector<size_t>& guided) const;
  void applyGuidedChoice(const std::vector<size_t>& guided, TrackCandidate& best) const;
  bool rejectedAsNeighbour(const TrackCandidate& best);
  bool nearAcceptedTrack(size_t hitIndex, double distanceSquaredMm2) const;
  void removeOutliers(TrackCandidate& best) const;
  void accept(const TrackCandidate& candidate);

  // ---- hit tests ----
  bool isolated(const std::vector<size_t>& combination) const;
  bool coincident(size_t first, size_t second) const;
  int completeDoublets(const std::vector<size_t>& combination) const;

  // ---- the event ----
  const Config& m_config;
  const FieldFunction& m_bzAt;
  const std::vector<EventHit>& m_hits;
  ExhaustiveStatistics& m_statistics;
  std::vector<bool> m_used;                                       ///< by an accepted track (or given up as a neighbour)
  std::vector<std::array<double, 3>> m_acceptedTrackHitPositions; ///< mm, for the neighbour-track cuts
  std::vector<TrackCandidate> m_accepted;                         ///< the tracks found, in order
  bool m_pairCoincidence = false;                                 ///< RequirePairCoincidence

  // ---- the current iteration ----
  HitPool m_pool;                                         ///< the free hits (hitsByLayer after MaxHitsPerCompositeID)
  std::vector<double> m_edepScore;                        ///< keV per hit (0 when not used)
  std::vector<std::vector<size_t>> m_isolationNeighbours; ///< pool hits closer than HitIsolationCut, other layer
  bool m_isolationNeighboursFound = false;                ///< found in the first iteration
  double m_absoluteFieldTesla = 0.0;                      ///< |Bz| at the pool centroid, for the pT
  double m_radiusAtMaxPtCm = 0.0;                         ///< radius of a track with pT = MaxComboPT
  int m_maxHitsPerCombination = 0;
  int m_minHitsPerCombination = 0;
};

std::vector<TrackCandidate> ExhaustiveSearch::run(const std::vector<EventHit>& hits,
                                                  ExhaustiveStatistics& statistics) const {
  EventSearch search(*this, hits, statistics);
  return search.run();
}

ExhaustiveSearch::EventSearch::EventSearch(const ExhaustiveSearch& search, const std::vector<EventHit>& hits,
                                           ExhaustiveStatistics& statistics)
    : m_config(search.m_config), m_bzAt(search.m_bzAt), m_hits(hits), m_statistics(statistics),
      m_used(hits.size(), false), m_pairCoincidence(search.m_config.requirePairCoincidence) {}

std::vector<TrackCandidate> ExhaustiveSearch::EventSearch::run() {
  while (true) {
    if (!collectHitPool(m_hits, m_used, m_config.common.maxHitEdepKeV, m_pool))
      break;
    prepareIteration();
    TrackCandidate best;
    if (!searchCombinations(best))
      break;
    if (m_config.doCrowdedLayerHitSelection)
      selectCrowdedLayerHits(best);
    if (rejectedAsNeighbour(best))
      continue;
    removeOutliers(best);
    accept(best);
  }
  return std::move(m_accepted);
}

// the accepted track: its hits are used, and remembered for the neighbour-track cuts
void ExhaustiveSearch::EventSearch::accept(const TrackCandidate& candidate) {
  for (size_t hitIndex : candidate.hits)
    m_used[hitIndex] = true;
  if (m_config.neighbourTrackMaxDistanceCm > 0.0) {
    for (size_t hitIndex : candidate.hits) {
      const auto& position = m_hits[hitIndex].positionMm;
      m_acceptedTrackHitPositions.push_back({position.x, position.y, position.z});
    }
  }
  m_accepted.push_back(candidate);
}

// ============================================================================================
// the free hits of one iteration
// ============================================================================================

// Exhaustive: field and pT cap, energy-deposit scores, size range of a combination, isolation
void ExhaustiveSearch::EventSearch::prepareIteration() {
  // field at the pool centroid: MaxComboPT as a radius cap for the k = 3 ranking
  m_absoluteFieldTesla = std::max(std::abs(fieldAtCentroid(m_hits, m_pool.hits, m_bzAt)), 1e-6);
  m_radiusAtMaxPtCm = radiusFromPt(m_absoluteFieldTesla, m_config.maxPtGeV);

  computeEdepScores();

  int availableHits = 0; // after MaxHitsPerCompositeID
  for (int layerKey : m_pool.layerOrder)
    availableHits += static_cast<int>(m_pool.hitsByLayer[layerKey].size());
  m_maxHitsPerCombination = std::min(availableHits, m_config.maxCombinatorialHits);
  m_minHitsPerCombination = std::max(3, m_config.common.minTrackHits);

  // the isolation neighbours once per event: later pools only lose hits (isolated() skips them)
  if (m_config.hitIsolationCutCm > 0.0 && !m_isolationNeighboursFound) {
    findIsolationNeighbours();
    m_isolationNeighboursFound = true;
  }
}

// Energy-deposit score: the right hit of a crowded layer usually deposits like the hits of the
// single-hit layers (the baseline), delta rays and secondaries much more or less.
//   score = distance of edep from [baseline, 2 x baseline] (keV, lower = better)
// MaxHitsPerCompositeID keeps only the best-scoring hits of crowded layers.
void ExhaustiveSearch::EventSearch::computeEdepScores() {
  std::vector<double> edepKeV(m_hits.size(), 0.0);
  for (size_t hitIndex : m_pool.hits)
    edepKeV[hitIndex] = m_hits[hitIndex].edepGeV * 1e6; // GeV -> keV

  // baseline: mean edep of the single-hit layers, else the median of all hits
  double baselineKeV = 0.0;
  int singleHitLayers = 0;
  for (const auto& [layerKey, layerHits] : m_pool.hitsByLayer) {
    if (layerHits.size() == 1) {
      baselineKeV += edepKeV[layerHits[0]];
      ++singleHitLayers;
    }
  }
  if (singleHitLayers > 0) {
    baselineKeV /= static_cast<double>(singleHitLayers);
  } else {
    std::vector<double> sortedEdeps;
    sortedEdeps.reserve(m_pool.hits.size());
    for (size_t hitIndex : m_pool.hits)
      sortedEdeps.push_back(edepKeV[hitIndex]);
    std::sort(sortedEdeps.begin(), sortedEdeps.end());
    const size_t middle = sortedEdeps.size() / 2;
    if (sortedEdeps.size() % 2 == 0)
      baselineKeV = 0.5 * (sortedEdeps[middle - 1] + sortedEdeps[middle]);
    else
      baselineKeV = sortedEdeps[middle];
  }

  m_edepScore.assign(m_hits.size(), 0.0);
  if (!(m_config.edepScoreWeight > 0.0))
    return;
  for (size_t hitIndex : m_pool.hits) {
    const double below = std::max(0.0, baselineKeV - edepKeV[hitIndex]);
    const double above = std::max(0.0, edepKeV[hitIndex] - 2.0 * baselineKeV);
    m_edepScore[hitIndex] = below + above;
  }
  const int maxHitsPerLayer = m_config.maxHitsPerLayer;
  if (maxHitsPerLayer <= 0)
    return;
  for (auto& [layerKey, layerHits] : m_pool.hitsByLayer) {
    if (static_cast<int>(layerHits.size()) <= maxHitsPerLayer)
      continue;
    std::sort(layerHits.begin(), layerHits.end(),
              [this](size_t left, size_t right) { return m_edepScore[left] < m_edepScore[right]; });
    layerHits.resize(maxHitsPerLayer);
  }
}

// HitIsolationCut: for every pool hit, the pool hits of other layers closer than the cut (cm)
// (same-layer hits, e.g. secondaries next to the track, do not count)
void ExhaustiveSearch::EventSearch::findIsolationNeighbours() {
  const double cut = m_config.hitIsolationCutCm;
  m_isolationNeighbours.assign(m_hits.size(), {});
  for (size_t hitIndex : m_pool.hits) {
    const auto position = m_hits[hitIndex].positionMm;
    const Eigen::Vector3d reference(position[0] / kMmPerCm, position[1] / kMmPerCm, position[2] / kMmPerCm);
    for (size_t otherIndex : m_pool.hits) {
      if (otherIndex == hitIndex || m_hits[otherIndex].layerKey == m_hits[hitIndex].layerKey)
        continue;
      const auto otherPosition = m_hits[otherIndex].positionMm;
      const Eigen::Vector3d other(otherPosition[0] / kMmPerCm, otherPosition[1] / kMmPerCm,
                                  otherPosition[2] / kMmPerCm);
      if ((other - reference).norm() < cut)
        m_isolationNeighbours[hitIndex].push_back(otherIndex);
    }
  }
}

// ============================================================================================
// Exhaustive: the best combination
// ============================================================================================

// all C(N,k) combinations for k = MaxCombinatorialHits down to MinTrackHits, after road, isolation
// and quality cuts, ranked by their circle fit
bool ExhaustiveSearch::EventSearch::searchCombinations(TrackCandidate& best) {
  const RoadCuts roadCuts(m_config.roadCuts, m_hits);

  HitCluster allHits; // all pool hits, layer by layer from the inside out
  allHits.layerKeys = m_pool.layerOrder;
  for (int layerKey : m_pool.layerOrder)
    allHits.hitsPerLayer.push_back(m_pool.hitsByLayer.at(layerKey));
  const std::vector<HitCluster> clusters = buildClusters(allHits);

  // pair coincidence: hits with a coincident partner first, then all hits if a doublet stays
  // incomplete; otherwise one pass over all hits
  std::vector<HitCluster> firstPassClusters;
  if (m_pairCoincidence)
    firstPassClusters.push_back(pairCoincidenceCluster(allHits));
  const int numberOfDoublets = static_cast<int>(m_config.doubletLayerPairs.size());
  SearchResult result;
  const int numberOfPasses = m_pairCoincidence ? 2 : 1;
  for (int pass = 0; pass < numberOfPasses; ++pass) {
    if (pass == 1 && result.found && completeDoublets(result.best.hits) >= numberOfDoublets)
      break;
    const bool coincidentHitsOnly = (m_pairCoincidence && pass == 0);
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
HitCluster ExhaustiveSearch::EventSearch::pairCoincidenceCluster(const HitCluster& allHits) const {
  HitCluster cluster;
  cluster.layerKeys = allHits.layerKeys;
  cluster.hitsPerLayer.assign(allHits.layerKeys.size(), std::vector<size_t>{});
  for (size_t layerIndex = 0; layerIndex < allHits.layerKeys.size(); ++layerIndex) {
    for (size_t hitIndex : allHits.hitsPerLayer[layerIndex]) {
      // partner layer key: the same key with the layer number swapped
      const int layerKey = m_hits[hitIndex].layerKey;
      const int layerNumber = layerNumberOfKey(layerKey);
      const auto partner = m_config.partnerLayer.find(layerNumber);
      if (partner == m_config.partnerLayer.end())
        continue; // not in a doublet
      const int absolutePartnerKey = std::abs(layerKey) - layerNumber + partner->second;
      const int partnerKey = (layerKey < 0) ? -absolutePartnerKey : absolutePartnerKey;
      const auto partnerHits = m_pool.hitsByLayer.find(partnerKey);
      if (partnerHits == m_pool.hitsByLayer.end())
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
std::vector<HitCluster> ExhaustiveSearch::EventSearch::buildClusters(const HitCluster& allHits) const {
  const int numberOfLayers = static_cast<int>(allHits.layerKeys.size());
  if (!(m_config.doHitClustering && numberOfLayers >= 2))
    return {allHits};

  // union-find on the pool position of a hit; clusters come in the (hash) order of their
  // representative position, the order the results always had
  std::vector<size_t> poolPosition(m_hits.size(), 0);
  for (size_t position = 0; position < m_pool.hits.size(); ++position)
    poolPosition[m_pool.hits[position]] = position;
  UnionFind groups(m_pool.hits.size());
  const double phiCut = m_config.roadCuts.maxConsecutiveDeltaPhi;
  const double distanceCut = m_config.roadCuts.maxConsecutiveHitDistanceCm; // cm
  joinLayers(groups, allHits, poolPosition, 1, phiCut, distanceCut);
  joinLayers(groups, allHits, poolPosition, 2, (phiCut > 0.0) ? phiCut * 2.0 : -1.0,
             (distanceCut > 0.0) ? distanceCut * 2.0 : -1.0); // across one missing layer, cuts doubled

  std::unordered_map<size_t, std::vector<size_t>> hitsOfGroup; // representative -> hits
  for (int layerIndex = 0; layerIndex < numberOfLayers; ++layerIndex) {
    for (size_t hitIndex : allHits.hitsPerLayer[layerIndex])
      hitsOfGroup[groups.find(poolPosition[hitIndex])].push_back(hitIndex);
  }

  std::vector<HitCluster> clusters; // with at least m_minHitsPerCombination layers, inner -> outer
  for (auto& [representative, groupHits] : hitsOfGroup) {
    std::unordered_map<int, std::vector<size_t>> hitsOfLayer;
    for (size_t hitIndex : groupHits)
      hitsOfLayer[m_hits[hitIndex].layerKey].push_back(hitIndex);
    if (static_cast<int>(hitsOfLayer.size()) < m_minHitsPerCombination)
      continue;
    HitCluster cluster;
    for (int layerKey : m_pool.layerOrder) {
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
void ExhaustiveSearch::EventSearch::joinLayers(UnionFind& groups, const HitCluster& allHits,
                                               const std::vector<size_t>& poolPosition, int step, double phiLimit,
                                               double distanceLimit) const {
  const int numberOfLayers = static_cast<int>(allHits.hitsPerLayer.size());
  for (int layerIndex = 0; layerIndex + step < numberOfLayers; ++layerIndex) {
    for (size_t innerHit : allHits.hitsPerLayer[layerIndex]) {
      const auto innerPosition = m_hits[innerHit].positionMm;
      for (size_t outerHit : allHits.hitsPerLayer[layerIndex + step]) {
        const auto outerPosition = m_hits[outerHit].positionMm;
        if (phiLimit > 0.0) {
          const double deltaPhi = wrapAngle(m_hits[outerHit].azimuth - m_hits[innerHit].azimuth);
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
void ExhaustiveSearch::EventSearch::searchCluster(const HitCluster& cluster, const RoadCuts& roadCuts,
                                                  SearchResult& result) const {
  const int numberOfLayers = static_cast<int>(cluster.layerKeys.size());
  const int largestSize = std::min(m_maxHitsPerCombination, numberOfLayers);
  std::vector<size_t> combination;
  for (int size = largestSize; size >= m_minHitsPerCombination; --size) {
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
      combination.assign(size, 0);
      chooseHits(0, hitsOfChosenLayer, roadCuts, combination, result);
    } while (nextLayerChoice(chosenLayers, numberOfLayers));
  }
}

// one hit per chosen layer, from layer `depth` on, depth first (the order of an odometer: the last
// layer turns first); a partial combination beyond the distance cuts is not completed (all its
// completions count as rejected), a complete one goes to the other cuts and the fit
void ExhaustiveSearch::EventSearch::chooseHits(size_t depth,
                                               const std::vector<const std::vector<size_t>*>& hitsOfChosenLayer,
                                               const RoadCuts& roadCuts, std::vector<size_t>& combination,
                                               SearchResult& result) const {
  if (depth == hitsOfChosenLayer.size()) {
    if (passesPreFitCuts(combination, roadCuts))
      offerCombination(combination, result);
    return;
  }
  for (size_t hitIndex : *hitsOfChosenLayer[depth]) {
    combination[depth] = hitIndex;
    if (!roadCuts.newHitPassesDistanceCuts(combination, depth)) {
      m_statistics.rejectedDistanceCuts += completionsAfter(hitsOfChosenLayer, depth);
      continue;
    }
    chooseHits(depth + 1, hitsOfChosenLayer, roadCuts, combination, result);
  }
}

// the cuts before the circle fit, cheapest first (the counters record the rejections)
bool ExhaustiveSearch::EventSearch::passesPreFitCuts(const std::vector<size_t>& combination,
                                                     const RoadCuts& roadCuts) const {
  if (!roadCuts.passesRoadCuts(combination)) {
    ++m_statistics.rejectedRoadCuts;
    return false;
  }
  if (!isolated(combination))
    return false; // HitIsolationCut
  if (!roadCuts.passesTanLambdaCut(combination)) {
    ++m_statistics.rejectedTanLambdaCut;
    return false;
  }
  if (m_pairCoincidence && completeDoublets(combination) < kMinCompleteDoublets)
    return false;
  return true;
}

// no hit of the combination has an isolation neighbour outside the combination
bool ExhaustiveSearch::EventSearch::isolated(const std::vector<size_t>& combination) const {
  if (m_config.hitIsolationCutCm <= 0.0)
    return true;
  for (size_t hitIndex : combination) {
    for (size_t neighbour : m_isolationNeighbours[hitIndex]) {
      if (m_used[neighbour])
        continue; // no longer in the pool
      if (std::find(combination.begin(), combination.end(), neighbour) == combination.end())
        return false;
    }
  }
  return true;
}

// two hits of a doublet are coincident within |dphi| < 0.15, |deta| < 0.10 and PairCoincidenceMaxDist
bool ExhaustiveSearch::EventSearch::coincident(size_t first, size_t second) const {
  const double maxDistanceMm = m_config.pairCoincidenceMaxDistanceCm * kMmPerCm;
  const auto firstPosition = m_hits[first].positionMm;
  const auto secondPosition = m_hits[second].positionMm;
  const double dx = secondPosition[0] - firstPosition[0];
  const double dy = secondPosition[1] - firstPosition[1];
  const double dz = secondPosition[2] - firstPosition[2];
  if ((dx * dx + dy * dy + dz * dz) > (maxDistanceMm * maxDistanceMm))
    return false;
  const double deltaPhi = wrapAngle(m_hits[second].azimuth - m_hits[first].azimuth);
  if (std::abs(deltaPhi) > kPairTolerancePhi)
    return false;
  const double firstRadius = m_hits[first].transverseRadiusMm;
  const double secondRadius = m_hits[second].transverseRadiusMm;
  if (firstRadius < 1e-6 || secondRadius < 1e-6)
    return false;
  const double firstEta = std::asinh(firstPosition[2] / firstRadius);
  const double secondEta = std::asinh(secondPosition[2] / secondRadius);
  return std::abs(secondEta - firstEta) <= kPairToleranceEta;
}

// complete doublets: hits in both layers of a doublet, in the same region, coincident
int ExhaustiveSearch::EventSearch::completeDoublets(const std::vector<size_t>& combination) const {
  int completeCount = 0;
  for (const auto& [firstLayer, secondLayer] : m_config.doubletLayerPairs) {
    bool complete = false;
    for (size_t firstHit : combination) {
      const int firstKey = m_hits[firstHit].layerKey;
      if (layerNumberOfKey(firstKey) != firstLayer)
        continue;
      for (size_t secondHit : combination) {
        const int secondKey = m_hits[secondHit].layerKey;
        if (secondHit == firstHit || layerNumberOfKey(secondKey) != secondLayer)
          continue;
        if (regionFromKey(secondKey) == regionFromKey(firstKey) && coincident(firstHit, secondHit)) {
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
void ExhaustiveSearch::EventSearch::offerCombination(const std::vector<size_t>& combination,
                                                     SearchResult& result) const {
  ++m_statistics.combinationsFitted;
  const std::optional<CircleFit> fit = fitCircle(m_hits, combination, m_config.common.sigmaHitDefaultCm);
  if (!fit)
    return;
  if (fit->chi2ndf >= m_config.maxCircleChi2NDF) {
    ++m_statistics.rejectedChi2Cut;
    return;
  }
  const int numberOfHits = static_cast<int>(combination.size());
  const double pT = ptFromRadius(m_absoluteFieldTesla, fit->radius); // GeV
  if (pT > m_config.maxPtGeV) {
    ++m_statistics.rejectedPtCuts;
    if (isBetterFallback(numberOfHits, *fit, result.fallback)) {
      result.fallback.hits = combination;
      result.fallback.circle = *fit;
      result.fallback.edepScoreKeV = 0.0;
    }
    return;
  }
  if (m_config.common.minPtGeV > 0.0 && pT < m_config.common.minPtGeV) {
    ++m_statistics.rejectedPtCuts; // below MinComboPT: no fallback
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
double ExhaustiveSearch::EventSearch::meanEdepScore(const std::vector<size_t>& combination) const {
  if (!(m_config.edepScoreWeight > 0.0))
    return 0.0;
  double sum = 0.0;
  for (size_t hitIndex : combination)
    sum += m_edepScore[hitIndex];
  return sum / static_cast<double>(combination.size());
}

// ranking (lower score is better): more hits always win; k >= 4: chi2/ndf + w * edep / EdepNormKeV;
// k = 3 (chi2 is 0): w * edep / EdepNormKeV - min(R, R_cap) / R_cap (straighter and MIP-like
// wins), with w = 0 the larger radius
bool ExhaustiveSearch::EventSearch::isBetterCandidate(int numberOfHits, const CircleFit& fit, double edepScoreKeV,
                                                      const TrackCandidate& best) const {
  if (numberOfHits != best.numberOfHits())
    return numberOfHits > best.numberOfHits();
  const double weight = m_config.edepScoreWeight;
  const double edepNorm = std::max(m_config.edepNormKeV, 1e-9);
  if (numberOfHits == 3) {
    if (!(weight > 0.0))
      return fit.radius > best.circle.radius;
    const double radiusTerm = std::min(fit.radius, m_radiusAtMaxPtCm) / m_radiusAtMaxPtCm; // saturates at R_cap
    const double bestRadiusTerm = std::min(best.circle.radius, m_radiusAtMaxPtCm) / m_radiusAtMaxPtCm;
    const double score = weight * edepScoreKeV / edepNorm - radiusTerm;
    const double bestScore = weight * best.edepScoreKeV / edepNorm - bestRadiusTerm;
    return score < bestScore;
  }
  const double score = fit.chi2ndf + (weight > 0.0 ? weight * edepScoreKeV / edepNorm : 0.0);
  const double bestScore = best.circle.chi2ndf + (weight > 0.0 ? weight * best.edepScoreKeV / edepNorm : 0.0);
  return score < bestScore;
}

// the best fallback: more hits, then lower chi2/ndf (k >= 4) or lower pT (k = 3)
bool ExhaustiveSearch::EventSearch::isBetterFallback(int numberOfHits, const CircleFit& fit,
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
void ExhaustiveSearch::EventSearch::selectCrowdedLayerHits(TrackCandidate& best) const {
  int cleanLayers = 0;
  int crowdedLayers = 0;
  for (const auto& [layerKey, layerHits] : m_pool.hitsByLayer) {
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
    const std::optional<CircleFit> cleanCircle = fitCircle(m_hits, cleanHits, m_config.common.sigmaHitDefaultCm);
    if (!cleanCircle || !(cleanCircle->radius > 0))
      return;
    const double distanceNorm = std::max(m_config.common.sigmaHitDefaultCm, 1e-9); // cm
    auto distanceToCircle = [&](size_t hitIndex) {
      const auto position = m_hits[hitIndex].positionMm; // mm
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
    const auto innerPosition = m_hits[cleanHits[0]].positionMm; // mm
    const auto outerPosition = m_hits[cleanHits[1]].positionMm;
    const double innerX = innerPosition.x * kCmPerMm, innerY = innerPosition.y * kCmPerMm,
                 innerZ = innerPosition.z * kCmPerMm;
    const double outerX = outerPosition.x * kCmPerMm, outerY = outerPosition.y * kCmPerMm,
                 outerZ = outerPosition.z * kCmPerMm;
    const double lineX = outerX - innerX, lineY = outerY - innerY, lineZ = outerZ - innerZ;
    const double lineLength = std::sqrt(lineX * lineX + lineY * lineY + lineZ * lineZ);
    if (!(lineLength > 1e-9))
      return;                                                                         // the two clean hits coincide
    const double kinkNorm = std::max(m_config.roadCuts.maxConsecutiveDeltaPhi, 1e-4); // rad
    auto kinkAngle = [&](size_t hitIndex) {              // angle at the inner clean hit between the line and the hit
      const auto position = m_hits[hitIndex].positionMm; // mm
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
std::vector<size_t> ExhaustiveSearch::EventSearch::cleanLayerHits() const {
  std::vector<size_t> cleanHits;
  for (const auto& [layerKey, layerHits] : m_pool.hitsByLayer) {
    if (layerHits.size() == 1)
      cleanHits.push_back(layerHits[0]);
  }
  return cleanHits;
}

// adds the lowest-score hit of every crowded layer (by layer key) to `guided`; false if a crowded
// layer has no candidate
bool ExhaustiveSearch::EventSearch::addBestHitOfCrowdedLayers(const std::function<double(size_t)>& geometricTerm,
                                                              std::vector<size_t>& guided) const {
  const double edepNorm = std::max(m_config.edepNormKeV, 1e-9); // keV
  for (const auto& [layerKey, layerHits] : m_pool.hitsByLayer) {
    if (layerHits.size() <= 1)
      continue;
    size_t bestHit = SIZE_MAX;
    double bestScore = std::numeric_limits<double>::max();
    for (size_t hitIndex : layerHits) {
      const double score =
          kGuidedGeometricWeight * geometricTerm(hitIndex) + kGuidedEdepWeight * (m_edepScore[hitIndex] / edepNorm);
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
void ExhaustiveSearch::EventSearch::applyGuidedChoice(const std::vector<size_t>& guided, TrackCandidate& best) const {
  const std::optional<CircleFit> fit = fitCircle(m_hits, guided, m_config.common.sigmaHitDefaultCm);
  if (!fit || fit->chi2ndf > m_config.maxCircleChi2NDF)
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
bool ExhaustiveSearch::EventSearch::rejectedAsNeighbour(const TrackCandidate& best) {
  if (!(m_config.neighbourTrackMaxDistanceCm > 0.0) || m_acceptedTrackHitPositions.empty())
    return false;
  const double distanceMm = m_config.neighbourTrackMaxDistanceCm * kMmPerCm;
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
  const bool enoughHits = (best.numberOfHits() >= m_config.neighbourTrackMinHits);
  const bool goodFit = (best.circle.chi2ndf <= m_config.neighbourTrackMaxChi2NDF);
  if (enoughHits && goodFit)
    return false;
  for (size_t hitIndex : best.hits) {
    if (nearAcceptedTrack(hitIndex, distanceSquared))
      m_used[hitIndex] = true;
  }
  return true;
}

bool ExhaustiveSearch::EventSearch::nearAcceptedTrack(size_t hitIndex, double distanceSquaredMm2) const {
  const auto position = m_hits[hitIndex].positionMm; // mm
  for (const auto& accepted : m_acceptedTrackHitPositions) {
    const double dx = position.x - accepted[0];
    const double dy = position.y - accepted[1];
    const double dz = position.z - accepted[2];
    if (dx * dx + dy * dy + dz * dz < distanceSquaredMm2)
      return true;
  }
  return false;
}

// drop the worst hit while its residual is above OutlierSigma * SigmaHitDefault, and refit
void ExhaustiveSearch::EventSearch::removeOutliers(TrackCandidate& best) const {
  const double outlierCutCm = m_config.outlierSigma * m_config.common.sigmaHitDefaultCm;
  for (int iteration = 0; iteration < m_config.maxOutlierIterations; ++iteration) {
    if (best.numberOfHits() <= m_minHitsPerCombination)
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
    const std::optional<CircleFit> refit = fitCircle(m_hits, best.hits, m_config.common.sigmaHitDefaultCm);
    if (refit)
      best.circle = *refit;
  }
}

// ============================================================================================
// the road cuts
// ============================================================================================

bool RoadCuts::passesRoadCuts(const std::vector<size_t>& hitIndices) const {
  const int numberOfHits = static_cast<int>(hitIndices.size());

  // near the beam (polar angle of the mean point, atan2(mean r, mean |z|)) the phi cuts are skipped
  double sumR = 0.0;
  double sumAbsZ = 0.0;
  for (size_t hitIndex : hitIndices) {
    sumR += m_hits[hitIndex].transverseRadiusMm;
    sumAbsZ += m_hits[hitIndex].absoluteZMm;
  }
  const double thetaCombination = std::atan2(sumR / numberOfHits, sumAbsZ / numberOfHits);
  // (atan2 of |z|: the angle to the beam on either side, 0 to 90 degrees)
  const double nearBeamLow = kNearBeamThetaMinDeg * M_PI / 180.0;
  if (thetaCombination < nearBeamLow)
    return true;

  // cut 1: |dphi| of consecutive layers (wrapped to (-pi, pi]); barrel <-> endcap pairs skipped
  if (m_config.maxConsecutiveDeltaPhi > 0.0) {
    for (int i = 0; i + 1 < numberOfHits; ++i) {
      const EventHit& inner = m_hits[hitIndices[i]];
      const EventHit& outer = m_hits[hitIndices[i + 1]];
      if ((inner.region == 0) != (outer.region == 0))
        continue;
      const double deltaPhi = std::abs(wrapAngle(outer.azimuth - inner.azimuth));
      if (deltaPhi > m_config.maxConsecutiveDeltaPhi)
        return false;
    }
  }

  // cut 2: phi spread, max - min of every phi relative to the innermost hit (window may cross +-pi)
  if (m_config.maxCombinationPhiSpread > 0.0) {
    const double innermostAzimuth = m_hits[hitIndices[0]].azimuth;
    double lowest = 0.0;
    double highest = 0.0;
    for (int i = 1; i < numberOfHits; ++i) {
      const double difference = wrapAngle(m_hits[hitIndices[i]].azimuth - innermostAzimuth);
      lowest = std::min(lowest, difference);
      highest = std::max(highest, difference);
    }
    if (highest - lowest > m_config.maxCombinationPhiSpread)
      return false;
  }

  // cut 3: in each endcap r can only grow with |z| for a track from the IP
  std::vector<std::pair<double, double>> positiveEndcapHits; // (|z|, r) in mm
  std::vector<std::pair<double, double>> negativeEndcapHits;
  for (size_t hitIndex : hitIndices) {
    const EventHit& hit = m_hits[hitIndex];
    if (hit.region == 1)
      positiveEndcapHits.emplace_back(hit.absoluteZMm, hit.transverseRadiusMm);
    if (hit.region == -1)
      negativeEndcapHits.emplace_back(hit.absoluteZMm, hit.transverseRadiusMm);
  }
  return radiusGrowsWithZ(positiveEndcapHits) && radiusGrowsWithZ(negativeEndcapHits);
}

bool RoadCuts::radiusGrowsWithZ(std::vector<std::pair<double, double>>& endcapHits) {
  if (endcapHits.size() < 2)
    return true;
  std::sort(endcapHits.begin(), endcapHits.end(),
            [](const auto& left, const auto& right) { return left.first < right.first; });
  for (size_t i = 0; i + 1 < endcapHits.size(); ++i) {
    if (endcapHits[i + 1].second - endcapHits[i].second < 0.0)
      return false;
  }
  return true;
}

double RoadCuts::distanceCm(size_t first, size_t second) const {
  const EventHit& firstHit = m_hits[first];
  const EventHit& secondHit = m_hits[second];
  const double dx = secondHit.positionMm.x - firstHit.positionMm.x, dy = secondHit.positionMm.y - firstHit.positionMm.y,
               dz = secondHit.positionMm.z - firstHit.positionMm.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz) * kCmPerMm;
}

bool RoadCuts::newHitPassesDistanceCuts(const std::vector<size_t>& hitIndices, size_t newHit) const {
  if (m_config.maxConsecutiveHitDistanceCm > 0.0 && newHit > 0) { // the hit before it
    if (distanceCm(hitIndices[newHit - 1], hitIndices[newHit]) > m_config.maxConsecutiveHitDistanceCm)
      return false;
  }
  if (m_config.maxPairHitDistanceCm > 0.0) { // every earlier hit
    for (size_t earlier = 0; earlier < newHit; ++earlier) {
      if (distanceCm(hitIndices[earlier], hitIndices[newHit]) > m_config.maxPairHitDistanceCm)
        return false;
    }
  }
  return true;
}

bool RoadCuts::passesTanLambdaCut(const std::vector<size_t>& hitIndices) const {
  const double maxDeviation = m_config.maxCombinationTanLambdaDeviation;
  const size_t numberOfHits = hitIndices.size();
  if (maxDeviation <= 0.0 || numberOfHits < 3)
    return true;

  // mixed barrel + endcap combinations have very different dz/dr per pair by geometry: they pass
  bool hasBarrel = false;
  bool hasEndcap = false;
  for (size_t hitIndex : hitIndices) {
    if (m_hits[hitIndex].region == 0)
      hasBarrel = true;
    else
      hasEndcap = true;
  }
  if (hasBarrel && hasEndcap)
    return true;

  // dip of each consecutive pair: dz / transverse distance (close to the arc length for near hits)
  std::vector<double> pairTanLambdas;
  pairTanLambdas.reserve(numberOfHits - 1);
  for (size_t i = 0; i + 1 < numberOfHits; ++i) {
    const EventHit& inner = m_hits[hitIndices[i]];
    const EventHit& outer = m_hits[hitIndices[i + 1]];
    const double transverseDistance =
        std::sqrt((outer.positionMm.x - inner.positionMm.x) * (outer.positionMm.x - inner.positionMm.x) +
                  (outer.positionMm.y - inner.positionMm.y) * (outer.positionMm.y - inner.positionMm.y));
    if (transverseDistance < 1.0) // mm
      continue;                   // too close in xy for a slope
    pairTanLambdas.push_back((outer.positionMm.z - inner.positionMm.z) / transverseDistance);
  }
  if (pairTanLambdas.size() < 2)
    return true;

  // every pair close to the median (not pulled by one wrong hit)
  std::vector<double> sorted = pairTanLambdas;
  std::sort(sorted.begin(), sorted.end());
  const double median = sorted[sorted.size() / 2];
  for (double pairTanLambda : pairTanLambdas) {
    if (std::abs(pairTanLambda - median) > maxDeviation)
      return false;
  }
  return true;
}

// ============================================================================================
// helpers (declared at the top)
// ============================================================================================

namespace {

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

  // number of complete combinations that share a partial combination ending at layer `depth`:
  // the product of the hit counts of the layers after it
  long completionsAfter(const std::vector<const std::vector<size_t>*>& hitsOfChosenLayer, size_t depth) {
    long completions = 1;
    for (size_t later = depth + 1; later < hitsOfChosenLayer.size(); ++later)
      completions *= static_cast<long>(hitsOfChosenLayer[later]->size());
    return completions;
  }

} // namespace

} // namespace CTF
