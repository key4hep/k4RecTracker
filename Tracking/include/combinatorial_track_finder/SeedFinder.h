/**
 * SeedFinder: the SeedAndExtend search of CombinatorialTrackFinder.
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

// SeedFinder - the SeedAndExtend strategy (best seeds first), for events too dense for the
// exhaustive search.
// 1. Seed windows from the geometry: straight rays from the origin swept in theta; the layers a ray
//    crosses give the windows (SeedLayerPolicy: first, last or every one; one layer may be skipped).
//    4-hit seeds first, then 3-hit seeds.
// 2. Seeds with the best kink chi2 (multiple scattering included) are extended first, layer by
//    layer along their helix; a hit is attached only if its triplet kink chi2 stays small.
// 3. An accepted track claims its hits at once. Every track gets a confirmation pass; short tracks
//    need a hole check and a minimum pT; sequential pieces of one trajectory are merged.
#ifndef TRACKING_SEEDFINDER_H
#define TRACKING_SEEDFINDER_H

#include "LayerGeometry.h"
#include "TrackCandidate.h"
#include "TrackingMath.h"

#include "edm4hep/TrackerHitPlaneCollection.h"

#include <Eigen/Dense>

#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace k4rt {

/// Where seeds are made, and so the extension direction (option SeedLayerPolicy).
enum class SeedLayerPolicy {
  Innermost,  ///< first layers a ray from the IP crosses; extended outward
  Outermost,  ///< last layers; extended inward
  Consecutive ///< any run of layers (displaced tracks); extended both ways
};

struct SeedConfig {
  /// shared with the exhaustive search (here minPt is required of short tracks only: a 3-hit
  /// track has one constraint, so random low-pT triplets pass any kink cut)
  CommonSearchConfig common;

  // --- from the CombinatorialTrackFinder options ---
  SeedLayerPolicy policy = SeedLayerPolicy::Innermost;
  int maxCandidatesPerLayer = 3;    ///< hits tried per crossed layer, nearest first
  int shortTrackMaxHits = 4;        ///< tracks up to this many hits need the hole check and minimum pT
  double shadowHitWindow_cm = 0.1;  ///< unclaimed hits this close to an accepted track are claimed (0 = off)
  double kinkMaxChi2Triplet = 40.0; ///< a hit is attached below this triplet kink chi2 (2 dof)
  double seedKinkMaxChi2NDF = 9.0;  ///< seeds are kept below this kink chi2/ndf
  double trackKinkMaxChi2NDF = 5.0; ///< a finished track is kept (and claims its hits) below this

  // --- fixed values ---
  std::vector<int> passes = {4, 3}; ///< seed sizes: 4-hit seeds, then 3-hit seeds on the hits left
  int claimMinHits = 4;             ///< in a pass other than the last, a track claims with >= 4 hits
  int confirmIterations = 2;        ///< confirmation passes (fill, swap, drop outliers)
  int layerSkip = 1;                ///< layers that may be skipped inside a seed window
  double maxZResid_cm = 0.5;        ///< max |z - line| of seed hits about their z-vs-arc-length line
  double window_cm = 1.0;           ///< hit search radius at a predicted crossing ...
  double windowSlope = 0.02;        ///< ... plus this many cm per cm of extrapolation
  int maxMissedLayers = 2;          ///< stop extending after this many crossed layers without a hit
  /// short tracks: at most this many crossed layers without a hit (inside the span and outward;
  /// inward is not counted, a displaced track has no hits inside its origin)
  int shortTrackMaxHoles = 1;
};

class SeedFinder {
public:
  /// Sensor plane of a hit: origin (cm) and unit normal; false if unknown.
  using PlaneFn = std::function<bool(uint64_t cellID, Eigen::Vector3d& origin, Eigen::Vector3d& normal)>;

  SeedFinder(const LayerGeometry& geometry, SeedConfig config, PlaneFn plane);

  /// Seed layer sets from the geometry (call once, in initialize()); returns the problems found.
  std::vector<std::string> prepare();
  const std::vector<int>& passes() const { return m_passes; } ///< seed size of each pass
  size_t numSets(int seedSize) const {                        ///< seed layer sets of this size
    auto entry = m_setsBySize.find(seedSize);
    return entry == m_setsBySize.end() ? 0 : entry->second.size();
  }

  /// Candidates among `pool` = (hit index, layer key) of the usable hits, longest first.
  std::vector<TrackCandidate> run(const edm4hep::TrackerHitPlaneCollection& hits,
                                  const std::vector<std::pair<size_t, int>>& pool) const;

private:
  class EventSearch; ///< the search in one event (SeedFinder.cpp)

  std::vector<std::vector<int>> buildSets(int seedSize) const; ///< from the ray sweep
  double closestApproachCm(int layerKey) const;                ///< for the inside-out window order

  const LayerGeometry& m_geometry;
  SeedConfig m_config;
  PlaneFn m_plane;
  std::map<int, std::vector<std::vector<int>>> m_setsBySize; ///< seed layer sets per seed size
  std::vector<int> m_passes;                                 ///< seed size of each pass
};

} // namespace k4rt

#endif // TRACKING_SEEDFINDER_H
