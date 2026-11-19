// TrackCandidate - what both search strategies of CombinatorialTrackFinder share: the hits of one
// event as the searches see them (EventHit, filled once per event), the pool of hits still free,
// the circle fit, the track candidate, and the EDM4hep track states built from a candidate.
// Units: hit positions in mm (EDM4hep); circle, covariance and residuals in cm; states in mm.
#ifndef TRACKING_TRACKCANDIDATE_H
#define TRACKING_TRACKCANDIDATE_H

#include "edm4hep/TrackState.h"
#include "edm4hep/Vector3d.h"

#include <Eigen/Dense>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <vector>

namespace CTF {

/// One input hit as the searches see it, filled once per event (index = hit index in the collection).
struct EventHit {
  bool usable = false; ///< on a known detector surface (the other hits are never used)
  uint64_t cellID = 0;
  edm4hep::Vector3d positionMm{};  ///< as in the input hit
  double transverseRadiusMm = 0.0; ///< sqrt(x^2 + y^2)
  double absoluteZMm = 0.0;        ///< |z|
  double azimuth = 0.0;            ///< atan2(y, x), rad
  int layerKey = 0;                ///< see LayerGeometry.h
  int region = 0;                  ///< 0 barrel, +1 / -1 endcaps
  double duMm = 0.0;               ///< measurement error along u (0 = not set)
  double dvMm = 0.0;               ///< measurement error along v (0 = not set)
  double varianceZMm2 = 0.0;       ///< zz element of the position covariance
  double edepGeV = 0.0;            ///< energy deposit
  // SeedAndExtend only: the sensor plane of the hit
  bool hasPlane = false;
  Eigen::Vector3d planeOriginCm = Eigen::Vector3d::Zero();
  Eigen::Vector3d planeNormal = Eigen::Vector3d::Zero(); ///< unit vector
};

/// The usable hits not yet taken by a track, by layer (one iteration of a search).
struct HitPool {
  std::vector<size_t> hits;                       ///< in collection order
  std::map<int, std::vector<size_t>> hitsByLayer; ///< layer key -> pool hits
  std::vector<int> layerOrder;                    ///< layer keys, by layer number (innermost first), then by key
};

/// Fills the pool with the unused usable hits below maxHitEdepKeV (<= 0: no limit); false (stop
/// the search) with fewer than 3 unused hits, pool hits or layers.
bool collectHitPool(const std::vector<EventHit>& hits, const std::vector<bool>& used, double maxHitEdepKeV,
                    HitPool& pool);

/// Bz (T) of the field map at a point (cm).
using FieldFunction = std::function<double(double xCm, double yCm, double zCm)>;

/// Bz (T) at the centroid of the hits hitIndices (at least one).
double fieldAtCentroid(const std::vector<EventHit>& hits, const std::vector<size_t>& hitIndices,
                       const FieldFunction& bzAt);

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
std::optional<CircleFit> fitCircle(const std::vector<EventHit>& hits, const std::vector<size_t>& hitIndices,
                                   double defaultSigmaCm);

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

// ---- EDM4hep track states ----

/// Where a stored helix state is placed.
enum class StatePoint { FirstHit, LastHit, Origin };

/// The fitted helix (circle + line of z versus transverse arc length) of the hits trackHits as an
/// EDM4hep track state at the point `where`; false (state unchanged) for degenerate input.
bool makeStateFromCircle(const CircleFit& circle, const std::vector<EventHit>& hits,
                         const std::vector<size_t>& trackHits, StatePoint where, edm4hep::TrackState& state);

/// A track state with these parameters (mm, rad, 1/mm) and a default diagonal covariance.
edm4hep::TrackState stateWithDefaultCovariance(double d0, double phi, double omega, double z0, double tanLambda,
                                               int location);

/// The solenoid as a cylinder closed by two disks (from the field map, CombinatorialTrackFinder).
struct Solenoid {
  double radiusMm = 0.0;
  double halfLengthMm = 0.0;
  double innerFieldTesla = 0.0; ///< Bz at the origin
  bool contains(const Eigen::Vector3d& positionCm) const;
};

/// The outer helix of a track measured outside the solenoid (SolenoidCrossing AtIP state), cm.
struct OuterHelix {
  double centreXCm = 0.0, centreYCm = 0.0, radiusCm = 0.0;  ///< fitted circle
  bool clockwise = false;                                   ///< seen from +z
  double ptGeV = 0.0;                                       ///< radius and field at the hits' centroid
  double tanLambda = 0.0;                                   ///< kept across the solenoid boundary
  Eigen::Vector3d innermostHitCm = Eigen::Vector3d::Zero(); ///< selects the crossing, anchors z
  bool innermostHitInEndcap = false;                        ///< then: the solenoid's end face
};

/// SolenoidCrossing: the outer helix followed back from its innermost hit to the solenoid boundary,
/// matched there to an inner helix (inner field, same pT and tanLambda) whose perigee is the AtIP
/// state; no value without a crossing.
std::optional<edm4hep::TrackState> solenoidCrossingState(const OuterHelix& outer, const Solenoid& solenoid);

} // namespace CTF

#endif // TRACKING_TRACKCANDIDATE_H
