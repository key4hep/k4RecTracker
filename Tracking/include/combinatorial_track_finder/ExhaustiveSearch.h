// ExhaustiveSearch - the Exhaustive strategy of CombinatorialTrackFinder (sparse events, e.g. a muon
// system). Per event, one track per iteration: all k-hit combinations of the free hits (k = most
// hits down to the fewest) pass the road cuts and a circle fit, the best is refined (crowded layers,
// neighbour tracks, outliers) and accepted, its hits are removed, and the search repeats.
// Positions in mm, distance cuts in cm. No Gaudi: the field comes in as a function.
#ifndef TRACKING_EXHAUSTIVESEARCH_H
#define TRACKING_EXHAUSTIVESEARCH_H

#include "TrackCandidate.h"
#include "TrackingMath.h"

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace CTF {

/// What the search did in one event (added to the Gaudi counters of CombinatorialTrackFinder).
struct ExhaustiveStatistics {
  long combinationsFitted = 0;
  long rejectedRoadCuts = 0;
  long rejectedDistanceCuts = 0;
  long rejectedTanLambdaCut = 0;
  long rejectedChi2Cut = 0;
  long rejectedPtCuts = 0;
};

/// Quick geometric tests that reject a combination which cannot be one track, before the circle
/// fit. A combination = hit indices, inner to outer layer.
class RoadCuts {
public:
  /// Settings of the road cuts (from the CombinatorialTrackFinder options, no defaults here); a
  /// value <= 0 disables that cut.
  struct Config {
    double maxConsecutiveDeltaPhi;           ///< rad, consecutive layers (MaxConsecDeltaPhi)
    double maxCombinationPhiSpread;          ///< rad, max - min phi (MaxComboPhiSpread)
    double maxConsecutiveHitDistanceCm;      ///< consecutive hits (MaxConsecutiveHitDistance)
    double maxPairHitDistanceCm;             ///< any two hits (MaxPairHitDistance)
    double maxCombinationTanLambdaDeviation; ///< hit-pair dz/dr from the median (MaxComboTanLambdaDev)
  };

  /// hits: the hits of the event (used by reference: must outlive this object)
  RoadCuts(const Config& config, const std::vector<EventHit>& hits) : m_config(config), m_hits(hits) {}

  /// |dphi| of consecutive layers, total phi spread (both skipped near the beam), and in each
  /// endcap a radius that does not decrease with |z|.
  bool passesRoadCuts(const std::vector<size_t>& hitIndices) const;
  /// 3D distance of the hit hitIndices[newHit] from the hit before it (consecutive hits) and from
  /// every earlier hit (any two hits): run for each new hit, all pairs are tested.
  bool newHitPassesDistanceCuts(const std::vector<size_t>& hitIndices, size_t newHit) const;
  /// The dip of each consecutive pair within maxCombinationTanLambdaDeviation of the median (pure
  /// barrel or pure endcap combinations only).
  bool passesTanLambdaCut(const std::vector<size_t>& hitIndices) const;

private:
  /// endcap hits as (|z|, r): does r never decrease with |z|? (sorts the list)
  static bool radiusGrowsWithZ(std::vector<std::pair<double, double>>& endcapHits);
  double distanceCm(size_t first, size_t second) const;

  Config m_config;
  const std::vector<EventHit>& m_hits;
};

class ExhaustiveSearch {
public:
  /// Settings of the exhaustive search (from the CombinatorialTrackFinder options, no defaults here).
  struct Config {
    CommonSearchConfig common;
    RoadCuts::Config roadCuts;
    int maxCombinatorialHits;                           ///< MaxCombinatorialHits
    double maxCircleChi2NDF;                            ///< NHitMaxChi2NDF
    double maxPtGeV;                                    ///< MaxComboPT
    double hitIsolationCutCm;                           ///< HitIsolationCut (0 = off)
    double outlierSigma;                                ///< OutlierSigma
    int maxOutlierIterations;                           ///< MaxOutlierIterations
    int maxHitsPerLayer;                                ///< MaxHitsPerCompositeID (0 = all)
    double edepScoreWeight;                             ///< EdepScoreWeight (0 = off)
    double edepNormKeV;                                 ///< EdepNormKeV
    bool doHitClustering;                               ///< DoHitClustering
    bool doCrowdedLayerHitSelection;                    ///< DoCrowdedLayerHitSelection
    bool requirePairCoincidence;                        ///< RequirePairCoincidence
    std::vector<std::pair<int, int>> doubletLayerPairs; ///< from DoubletLayers
    std::map<int, int> partnerLayer;                    ///< layer number -> the other layer of its doublet
    double pairCoincidenceMaxDistanceCm;                ///< PairCoincidenceMaxDist
    double neighbourTrackMaxDistanceCm;                 ///< NeighbourTrackMaxDist (0 = off)
    int neighbourTrackMinHits;                          ///< NeighbourTrackMinHits
    double neighbourTrackMaxChi2NDF;                    ///< NeighbourTrackMaxChi2NDF
  };

  ExhaustiveSearch(Config config, FieldFunction bzAt);

  /// The accepted tracks of one event, in the order they were found.
  std::vector<TrackCandidate> run(const std::vector<EventHit>& hits, ExhaustiveStatistics& statistics) const;

private:
  class EventSearch; ///< the search in one event (ExhaustiveSearch.cpp)

  Config m_config;
  FieldFunction m_bzAt;
};

} // namespace CTF

#endif // TRACKING_EXHAUSTIVESEARCH_H
