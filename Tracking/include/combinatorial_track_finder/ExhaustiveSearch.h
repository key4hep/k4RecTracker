// ExhaustiveSearch - the parts of CombinatorialTrackFinder's exhaustive search that need no event
// state: the per-hit numbers (HitGeometry, computed once per event) and the road cuts (RoadCuts)
// that reject a combination of hits before its fit. Positions in mm, distance cuts in cm.
#ifndef TRACKING_EXHAUSTIVESEARCH_H
#define TRACKING_EXHAUSTIVESEARCH_H

#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

namespace k4rt {

/// A hit as the road cuts see it (positions in mm).
struct HitGeometry {
  double x = 0.0, y = 0.0, z = 0.0;
  double transverseRadius = 0.0; ///< sqrt(x^2 + y^2)
  double absoluteZ = 0.0;        ///< |z|
  double azimuth = 0.0;          ///< atan2(y, x), rad
  int layerKey = 0;              ///< see LayerGeometry.h
  int region = 0;                ///< regionOfLayer(layerKey): 0 barrel, +1 / -1 endcaps
};

/// Settings of the road cuts; a value <= 0 disables that cut.
struct RoadCutConfig {
  double maxConsecDeltaPhi = -1.0;    ///< rad, consecutive layers (MaxConsecDeltaPhi)
  double maxComboPhiSpread = -1.0;    ///< rad, max - min phi (MaxComboPhiSpread)
  double maxConsecHitDist = -1.0;     ///< cm, consecutive hits (MaxConsecutiveHitDistance)
  double maxPairHitDist = -1.0;       ///< cm, any two hits (MaxPairHitDistance)
  double maxComboTanLambdaDev = -1.0; ///< hit-pair dz/dr from the median of all pairs (MaxComboTanLambdaDev)
  // fixed values: near the beam (mean polar angle outside 15-160 degrees) the phi cuts are skipped
  // (low-pT tracks turn a lot between endcap layers); pairs closer than 1 mm in xy skip the dz/dr cut
  double nearBeamThetaMinDeg = 15.0;
  double nearBeamThetaMaxDeg = 160.0;
  double minPairTransverseDistanceMm = 1.0;
};

/// Region of a layer key (0 barrel, +1 / -1 endcaps): from the readout's region field collected at
/// initialize() (regionOfLayerKey), else from the key itself (regionFromKey).
int regionOfLayer(int layerKey, const std::unordered_map<int, int>* regionOfLayerKey);

/// Quick geometric tests that reject a combination which cannot be one track, before the circle
/// fit. A combination = hit indices, inner to outer layer.
class RoadCuts {
public:
  /// hitGeometry: one entry per hit of the event (used by reference: must outlive this object)
  RoadCuts(const RoadCutConfig& config, const std::vector<HitGeometry>& hitGeometry)
      : m_config(config), m_hitGeometry(hitGeometry) {}

  /// |dphi| of consecutive layers, total phi spread (both skipped near the beam), and in each
  /// endcap a radius that does not decrease with |z|.
  bool passesRoadCuts(const std::vector<size_t>& hitIndices) const;
  /// 3D distance of consecutive hits and of any two hits.
  bool passesDistanceCuts(const std::vector<size_t>& hitIndices) const;
  /// The dip of each consecutive pair within maxComboTanLambdaDev of the median (pure barrel or
  /// pure endcap combinations only).
  bool passesTanLambdaCut(const std::vector<size_t>& hitIndices) const;

private:
  /// endcap hits as (|z|, r): does r never decrease with |z|? (sorts the list)
  static bool radiusGrowsWithZ(std::vector<std::pair<double, double>>& endcapHits);
  double distanceCm(size_t first, size_t second) const;

  RoadCutConfig m_config;
  const std::vector<HitGeometry>& m_hitGeometry;
};

} // namespace k4rt

#endif // TRACKING_EXHAUSTIVESEARCH_H
