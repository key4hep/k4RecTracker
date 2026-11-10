/**
 * SeedFinder: the SeedAndExtend search of CombinatorialTrackFinder.
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

// SeedFinder - the SeedAndExtend strategy of CombinatorialTrackFinder (best seeds first), for
// events too dense for the exhaustive C(N,k) search over all layers.
//
// 1. Seed layer sets are derived once from LayerGeometry: straight rays from the origin are
//    swept in theta; the layers a ray crosses give the seed windows. SeedLayerPolicy picks the
//    first window ("innermost"), the last ("outermost") or every one ("consecutive"); one layer
//    may be skipped inside a window. Seeds of 4 hits are tried first, then of 3.
// 2. Within a window the seeds best in kink chi2 (multiple scattering included, TrackingMath.h)
//    are extended first, layer by layer along their helix (HelixExtrapolator): outward for
//    innermost seeds, inward for outermost seeds, both ways for consecutive ones. A hit is
//    attached only if the kink chi2 of its triplet stays small.
// 3. An accepted track claims its hits at once; later seeds that use them are skipped.
//    Each track gets a confirmation pass, short tracks must pass a hole check and a minimum pT,
//    and sequential pieces of one trajectory are merged.
#ifndef TRACKING_SEEDFINDER_H
#define TRACKING_SEEDFINDER_H

#include "HelixExtrapolator.h"
#include "LayerGeometry.h"
#include "TrackingMath.h"

#include "edm4hep/TrackerHitPlane.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include <Eigen/Dense>

#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace k4rt {

struct SeedConfig {
  /// shared with the exhaustive search: minimum hits; minimum pT (here required of tracks with
  /// at most shortTrackMaxHits hits: a 3-hit track has one constraint, so random low-pT triplets
  /// pass any kink cut); largest distance between consecutive seed hits; field
  CommonSearchConfig common;

  // --- set from the CombinatorialTrackFinder options ---
  std::string policy = "innermost"; ///< innermost | outermost | consecutive (also sets the extension direction)
  int maxCandPerLayer = 3;          ///< MaxCandidatesPerLayer: hits tried per crossed layer, nearest first
  /// ShortTrackMaxHits: a track with at most this many hits must pass the hole check and the minimum pT
  int shortTrackMaxHits = 4;
  /// ShadowHitWindow: after a track is accepted, unclaimed hits on its trajectory within this
  /// distance (cm, on their own sensor plane) are claimed too (overlapping-sensor double hits,
  /// which would otherwise seed a parallel duplicate). 0 = off.
  double shadowHitWindow_cm = 0.1;
  // kink chi2 cuts (multiple scattering included, TrackingMath.h)
  double kinkMaxChi2Triplet = 40.0; ///< a hit is attached only if its triplet's kink chi2 (2 dof) is below this
  double seedKinkMaxChi2NDF = 9.0;  ///< seed acceptance (kink chi2/ndf)
  double trackKinkMaxChi2NDF = 5.0; ///< a finished track is kept (and claims its hits) below this

  // --- fixed values (no option) ---
  std::vector<int> passes = {4, 3}; ///< seed sizes: 4-hit seeds first, then 3-hit seeds on the hits left
  int claimMinHits = 4;             ///< in the 4-hit pass a candidate claims its hits if it has >= 4
  int confirmIterations = 2;        ///< confirmation passes after extension (fill, swap, drop outliers)
  int layerSkip = 1;                ///< layers that may be skipped inside one seed window
  double maxZResid_cm = 0.5;        ///< max |z - line| of seed hits about their z-vs-arc-length line
  double window_cm = 1.0;           ///< hit search radius at a predicted crossing ...
  double windowSlope = 0.02;        ///< ... growing by this many cm per cm of extrapolation
  int maxMissedLayers = 2;          ///< stop extending after this many crossed layers without a hit
  /// Hole check of short tracks: at most this many layers crossed by the helix (between the first
  /// and last hit, and outward until it leaves the tracker) without a hit near the crossing.
  /// Inward holes are not counted: a displaced track has none inside its origin.
  int shortTrackMaxHoles = 1;
};

/// A found track candidate. Hits are ordered along the direction of travel.
struct SECandidate {
  std::vector<size_t> gIdx;
  double chi2ndf = 0, x0 = 0, y0 = 0, R = 0;
  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  std::vector<double> residuals;
  int nSeedHits = 0;
  double kinkChi2 = -1.0; ///< multiple-scattering-aware kink chi2 (TrackingMath.h)
  int kinkNdf = 0;
  bool looper = false; ///< turns by more than half a circle, or overlaps another arc of the same circle
  bool merged = false; ///< absorbed into another candidate (internal)
};

struct SeedStats {
  long seedLayerSets = 0;
  long rawTriplets = 0; ///< hit combinations enumerated, before any check
  long seeds = 0;       ///< combinations passing the cheap geometric pre-check
  long seedsAfterZ = 0, extended = 0, candidates = 0;
  long skippedClaimed = 0, rejectedTracks = 0, merged = 0, loopers = 0;
  long confirmAdded = 0, confirmSwapped = 0, confirmRemoved = 0, rejectedHoles = 0, rejectedShortPt = 0;
  long shadowClaimed = 0;
};

class SeedFinder {
public:
  /// Circle fit: (ordered hits) -> x0, y0, R (cm), chi2/ndf, 3x3 cov, residuals. This is the
  /// component's fitCircleNHits, so there is one fit, not two.
  using FitFn = std::function<bool(const std::vector<edm4hep::TrackerHitPlane>&, double&, double&, double&, double&,
                                   Eigen::Matrix3d&, std::vector<double>&)>;
  /// Sensor plane of a hit: origin (cm) and unit normal. false if unknown.
  using PlaneFn = std::function<bool(uint64_t cellID, Eigen::Vector3d&, Eigen::Vector3d&)>;

  SeedFinder(const LayerGeometry& geo, SeedConfig cfg, FitFn fit, PlaneFn plane);

  /// Derive the seed layer sets from the geometry. Call once (initialize()).
  /// Returns problems found (e.g. an unknown policy).
  std::vector<std::string> prepare();
  const std::vector<std::vector<int>>& seedLayerSets() const { return m_sets; } ///< first pass
  const std::vector<int>& passes() const { return m_passes; }
  size_t numSets(int n) const {
    auto it = m_setsByN.find(n);
    return it == m_setsByN.end() ? 0 : it->second.size();
  }
  const SeedConfig& config() const { return m_cfg; }

  /// Find candidates among `pool` = (global hit index, compositeID) of usable hits.
  std::vector<SECandidate> run(const edm4hep::TrackerHitPlaneCollection& hits,
                               const std::vector<std::pair<size_t, int>>& pool, SeedStats& stats) const;

private:
  struct LayerHits { // the unclaimed hits of one layer (indices into the hit collection)
    std::vector<size_t> gIdx;
  };
  /// The search in one event: holds the event's hits, the claimed hits and the accepted
  /// candidates, and does the steps of run() as its methods (defined in SeedFinder.cpp).
  class EventSearch;

  /// Seed layer sets of size n from the ray sweep (policy consecutive/innermost/outermost).
  std::vector<std::vector<int>> buildSets(int n, std::vector<std::string>& problems) const;

  const LayerGeometry& m_geo;
  SeedConfig m_cfg;
  FitFn m_fit;
  PlaneFn m_plane;
  std::vector<std::vector<int>> m_sets;                   ///< sets of the first pass
  std::map<int, std::vector<std::vector<int>>> m_setsByN; ///< seed layer sets per seed size
  std::vector<int> m_passes;                              ///< seed size of each pass
};

} // namespace k4rt

#endif // TRACKING_SEEDFINDER_H
