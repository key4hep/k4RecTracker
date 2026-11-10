// ExhaustiveSearch - the parts of CombinatorialTrackFinder's exhaustive search that do not
// depend on the per-event search state: the hit and combination records, and the road cuts
// (class RoadCuts) that reject a combination of hits before it is fitted.
//
// Units: hit positions in mm (EDM4hep), circle parameters and distance cuts in cm, B in tesla.
#ifndef TRACKING_EXHAUSTIVESEARCH_H
#define TRACKING_EXHAUSTIVESEARCH_H

#include "LayerGeometry.h"
#include "TrackingMath.h"

#include "edm4hep/TrackerHitPlane.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include <Eigen/Dense>

#include <cstddef>
#include <limits>
#include <unordered_map>
#include <vector>

namespace k4rt {

/// One unused hit offered to the search, with its decoded layer identity.
struct HitInfo {
  size_t index; ///< index into the full input TrackerHitPlaneCollection
  edm4hep::TrackerHitPlane hit;
  int layerID;
  int typeID;
  int compositeID;

  HitInfo(size_t idx, const edm4hep::TrackerHitPlane& h, int layer, int type, int composite)
      : index(idx), hit(h), layerID(layer), typeID(type), compositeID(composite) {}
};

/// A fitted hit combination and the track parameters derived from it.
struct CombResult {
  std::vector<size_t> hiiVec;  // indices into allHitInfo
  std::vector<size_t> gIdxVec; // indices into the full hit collection
  double chi2ndf = std::numeric_limits<double>::max();
  double avgProxCm = std::numeric_limits<double>::max(); // avg edep-deviation score (keV)
  double x0 = 0, y0 = 0, radius = 0;
  Eigen::Matrix3d fitCov = Eigen::Matrix3d::Zero();
  std::vector<double> residuals;
  int nHits = 0;
  double pT = 0.0; // GeV, from the radius and the field at the centroid of the hits
  // SeedAndExtend only: multiple-scattering-aware kink chi2, looper flag
  double kinkChi2 = -1.0;
  int kinkNdf = 0;
  bool looper = false;
};

/// Settings of the road cuts. A cut value <= 0 disables that cut.
struct RoadCutConfig {
  double maxConsecDeltaPhi = -1.0;    ///< rad: |dphi| between hits of consecutive layers (MaxConsecDeltaPhi)
  double maxComboPhiSpread = -1.0;    ///< rad: max - min phi of all hits (MaxComboPhiSpread)
  double maxConsecHitDist = -1.0;     ///< cm: 3D distance between consecutive hits (MaxConsecutiveHitDistance)
  double maxPairHitDist = -1.0;       ///< cm: 3D distance between any two hits (MaxPairHitDistance)
  double maxComboTanLambdaDev = -1.0; ///< largest deviation of a hit pair's dz/dr from the median of all pairs

  // fixed values (no option)
  /// near the beam (mean polar angle of the hits below / above these, degrees) the phi cuts are
  /// skipped: low-pT tracks there turn by a large angle in phi between endcap layers
  double nearBeamThetaMinDeg = 15.0;
  double nearBeamThetaMaxDeg = 160.0;
  /// hit pairs closer than this in the transverse plane (mm) are not used in the dz/dr cut
  double minPairTransverseDistanceMm = 1.0;
};

/// Region of a compositeID: 0 barrel, +1 +endcap, -1 -endcap.
///
/// The answer comes from the readout's region field, collected per compositeID at initialize()
/// (pass that map here); without a map, or for an unknown cid, from the key (regionOfKey).
int regionOfCid(int cid, const std::unordered_map<int, int>* cidRegion);

/// The road cuts of the exhaustive search: quick geometric tests that reject a combination of
/// hits which cannot belong to one track, before the (slower) circle fit. A combination is given
/// as indices into the event's hit collection, ordered from the inner to the outer layer.
class RoadCuts {
public:
  /// compositeOfHit: hit index -> compositeID (layer key); cidRegion: compositeID -> region
  /// (see regionOfCid). The collections are used by reference and must outlive this object.
  RoadCuts(const RoadCutConfig& config, const edm4hep::TrackerHitPlaneCollection& hits,
           const std::unordered_map<size_t, int>& compositeOfHit, const std::unordered_map<int, int>* cidRegion)
      : m_config(config), m_hits(hits), m_compositeOfHit(compositeOfHit), m_cidRegion(cidRegion) {}

  /// Azimuthal road: |dphi| between consecutive layers, the total phi spread, and in each endcap
  /// a radius that does not decrease with |z|. The phi cuts are skipped near the beam.
  bool passesRoadCuts(const std::vector<size_t>& hitIndices) const;

  /// 3D distance between consecutive hits, and between any two hits.
  bool passesDistanceCuts(const std::vector<size_t>& hitIndices) const;

  /// The dip of every pair of consecutive hits (dz over their transverse distance) must stay
  /// within maxComboTanLambdaDev of the median of all pairs. Mixed barrel + endcap combinations pass.
  bool passesTanLambdaCut(const std::vector<size_t>& hitIndices) const;

private:
  /// endcap hits as (|z|, r) in mm: does r never decrease when |z| grows? (sorts the list by |z|)
  static bool radiusGrowsWithZ(std::vector<std::pair<double, double>>& endcapHits);
  /// 3D distance (cm) between two hits (positions in mm)
  double distanceCm(size_t first, size_t second) const;

  RoadCutConfig m_config;
  const edm4hep::TrackerHitPlaneCollection& m_hits;
  const std::unordered_map<size_t, int>& m_compositeOfHit;
  const std::unordered_map<int, int>* m_cidRegion;
};

} // namespace k4rt

#endif // TRACKING_EXHAUSTIVESEARCH_H
