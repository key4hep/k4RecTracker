/**
 * HelixExtrapolator: helix model and layer crossings of CombinatorialTrackFinder.
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

// The class Helix answers the question of the seed-and-extend search: where does a fitted track
// next cross a barrel cylinder, an endcap disk, or one sensor plane, going forward or backward?
//
// The helix is defined GEOMETRICALLY, without charge or field sign:
//   - a circle in the transverse plane: centre (cx, cy) and radius R,
//   - the rotation sense `rot` of travel about the centre (+1: the angle about the centre grows
//     along the track), taken from the order of the hits,
//   - a straight line in z: z = zRef + tanL * s, where s is the transverse arc length measured
//     from a reference point (the first hit) along the direction of travel.
// A point at arc length s is at the angle phiRef + rot * s / R about the centre.
//
// Uniform-field model (a helix). Units: cm. No Gaudi/DD4hep dependency, so it is unit-testable on
// its own (Tracking/test/testCombinatorialTrackFinder/test_HelixExtrapolator.cpp).
#ifndef TRACKING_HELIXEXTRAPOLATOR_H
#define TRACKING_HELIXEXTRAPOLATOR_H

#include "LayerGeometry.h"

#include <Eigen/Dense>

#include <limits>
#include <vector>

namespace k4rt {

/// Where a helix crosses a surface: ok = false if it does not.
struct Crossing {
  bool ok = false;
  double s = 0; ///< arc length at the crossing (cm)
  Eigen::Vector3d pos = Eigen::Vector3d::Zero();
};

class Helix {
public:
  double cx = 0, cy = 0, R = 0; ///< transverse circle, cm
  int rot = +1;                 ///< +1: the angle about the centre increases along the travel
  double phiRef = 0;            ///< angle of the reference point about the centre
  double zRef = 0;              ///< z at the reference point, cm
  double tanL = 0;              ///< dz / ds (s = transverse arc length)

  /// Build a helix from a fitted circle (centre, radius) and hits ordered ALONG THE TRACK (cm).
  /// The reference point is the first hit; rot comes from the hit order; zRef and tanL from a
  /// straight-line fit of z versus s, with s unwrapped hit to hit (so tracks turning by more than
  /// pi are handled). The arc length of every hit is returned in sOfHits if given.
  /// Returns false with fewer than 2 hits or R <= 0.
  static bool fromHits(double centreX, double centreY, double radius, const std::vector<Eigen::Vector3d>& orderedHitsCm,
                       Helix& out, std::vector<double>* sOfHits = nullptr);

  /// Angle about the centre at arc length s.
  double angleAt(double s) const { return phiRef + rot * s / R; }
  /// Position at arc length s from the reference point.
  Eigen::Vector3d at(double s) const;
  /// Unit 3D direction of travel at arc length s.
  Eigen::Vector3d dirAt(double s) const;

  /// Arc length of the point of the helix closest in angle to p, taken as the one nearest to
  /// sGuess (used to place a hit on the helix).
  double arcLengthNear(const Eigen::Vector3d& p, double sGuess) const;

  /// First crossing of the transverse radius rCylinder, moving from sFrom in the direction dir
  /// (+1 forward, -1 backward), at least sEps and at most sMax away in |s|.
  Crossing crossCylinder(double rCylinder, double sFrom, int dir, double sEps = 1e-6,
                         double sMax = std::numeric_limits<double>::infinity()) const;

  /// Crossing of the plane z = zPlane (unique: z changes monotonically with s). Same s window.
  Crossing crossPlaneZ(double zPlane, double sFrom, int dir, double sEps = 1e-6,
                       double sMax = std::numeric_limits<double>::infinity()) const;

  /// Crossing of a layer within its extent (+- tolerance): the cylinder at its mean radius for a
  /// barrel layer, the plane at its mean z for a disk. Approximate by construction (a barrel
  /// layer is a polygon of flat staves): use refineToPlane() against the actual sensor.
  Crossing crossLayer(const LayerDescriptor& layer, double sFrom, int dir, double toleranceCm,
                      double sMax = std::numeric_limits<double>::infinity()) const;

  /// Exact crossing of the plane through `origin` with unit normal `normal`, by Newton iteration
  /// starting from sStart. ok = false if it does not converge.
  Crossing refineToPlane(const Eigen::Vector3d& origin, const Eigen::Vector3d& normal, double sStart,
                         int maxIterations = 20, double toleranceCm = 1e-9) const;
};

} // namespace k4rt

#endif // TRACKING_HELIXEXTRAPOLATOR_H
