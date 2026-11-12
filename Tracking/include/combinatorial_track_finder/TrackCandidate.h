// TrackCandidate - a track candidate of CombinatorialTrackFinder (its hits and their transverse
// circle) and the circle fit itself. Both strategies produce TrackCandidates.
// Units: hit positions in mm (EDM4hep); circle, covariance and residuals in cm.
#ifndef TRACKING_TRACKCANDIDATE_H
#define TRACKING_TRACKCANDIDATE_H

#include "edm4hep/TrackerHitPlaneCollection.h"

#include <Eigen/Dense>

#include <cstddef>
#include <limits>
#include <optional>
#include <vector>

namespace k4rt {

/// A circle fitted to the transverse positions of hits.
struct CircleFit {
  double centreX = 0.0;                                 ///< cm
  double centreY = 0.0;                                 ///< cm
  double radius = 0.0;                                  ///< cm
  double chi2ndf = std::numeric_limits<double>::max();  ///< ndf = hits - 3 (at least 1)
  Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero(); ///< of (centreX, centreY, radius), cm^2
  std::vector<double> residuals;                        ///< distance of every hit from the circle, cm
};

/// Circle fit to the hits hitIndices (at least 3): weighted algebraic first guess, then
/// Gauss-Newton on the real hit distances. Hit weight from du (defaultSigmaCm if du is not set).
/// No value when degenerate (fewer than 3 hits, no real circle, radius <= 0).
std::optional<CircleFit> fitCircle(const edm4hep::TrackerHitPlaneCollection& hits,
                                   const std::vector<size_t>& hitIndices, double defaultSigmaCm);

struct TrackCandidate {
  /// hit indices: inner to outer layer (Exhaustive) or along the travel (SeedAndExtend)
  std::vector<size_t> hits;
  CircleFit circle;
  double edepScoreKeV = std::numeric_limits<double>::max(); ///< Exhaustive: mean energy-deposit score (ranking)
  double kinkChi2 = -1.0; ///< SeedAndExtend: kink chi2 with multiple scattering (TrackingMath.h)
  int kinkNdf = 0;
  bool looper = false; ///< SeedAndExtend: turns > half a circle, or overlaps another arc of the circle

  int numberOfHits() const { return static_cast<int>(hits.size()); }
};

} // namespace k4rt

#endif // TRACKING_TRACKCANDIDATE_H
