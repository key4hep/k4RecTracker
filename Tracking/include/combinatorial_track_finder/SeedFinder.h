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

#include <Eigen/Dense>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace CTF {

/// Where seeds are made, and so the extension direction (option SeedLayerPolicy).
enum class SeedLayerPolicy {
  Innermost,  ///< first layers a ray from the IP crosses; extended outward
  Outermost,  ///< last layers; extended inward
  Consecutive ///< any run of layers (displaced tracks); extended both ways
};

class SeedFinder {
public:
  /// Settings of the search (from the CombinatorialTrackFinder options, no defaults here); the
  /// fixed values of the search are constants in SeedFinder.cpp.
  struct Config {
    /// shared with the exhaustive search (here minPt is required of short tracks only: a 3-hit
    /// track has one constraint, so random low-pT triplets pass any kink cut)
    CommonSearchConfig common;
    SeedLayerPolicy policy;    ///< SeedLayerPolicy
    int maxCandidatesPerLayer; ///< MaxCandidatesPerLayer: hits tried per crossed layer, nearest first
    int shortTrackMaxHits;     ///< ShortTrackMaxHits: tracks up to this many hits need the hole check and minimum pT
    double
        shadowHitWindow_cm; ///< ShadowHitWindow: unclaimed hits this close to an accepted track are claimed (0 = off)
    double kinkMaxChi2Triplet;  ///< KinkMaxChi2Triplet: a hit is attached below this triplet kink chi2 (2 dof)
    double seedKinkMaxChi2NDF;  ///< SeedKinkMaxChi2NDF: seeds are kept below this kink chi2/ndf
    double trackKinkMaxChi2NDF; ///< TrackKinkMaxChi2NDF: a finished track is kept (and claims its hits) below this
  };

  SeedFinder(const LayerGeometry& geometry, Config config);

  /// Seed layer sets from the geometry (call once, in initialize()); returns the problems found.
  std::vector<std::string> prepare();
  const std::vector<int>& passes() const { return m_passes; } ///< seed size of each pass
  size_t numSets(int seedSize) const {                        ///< seed layer sets of this size
    auto entry = m_setsBySize.find(seedSize);
    return entry == m_setsBySize.end() ? 0 : entry->second.size();
  }

  /// Candidates among the pool hits (indices into hits, sensor planes filled), longest first.
  std::vector<TrackCandidate> run(const std::vector<EventHit>& hits, const std::vector<size_t>& pool) const;

private:
  class EventSearch; ///< the search in one event (SeedFinder.cpp)

  std::vector<std::vector<int>> buildSets(int seedSize) const; ///< from the ray sweep
  double closestApproachCm(int layerKey) const;                ///< for the inside-out window order

  const LayerGeometry& m_geometry;
  Config m_config;
  std::map<int, std::vector<std::vector<int>>> m_setsBySize; ///< seed layer sets per seed size
  std::vector<int> m_passes;                                 ///< seed size of each pass
};

} // namespace CTF

#endif // TRACKING_SEEDFINDER_H
