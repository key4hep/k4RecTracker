/**
 * HelixExtrapolator: helix model and layer crossings of CombinatorialTrackFinder.
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

// HelixExtrapolator - forward/backward propagation of a circle+linear-z helix to layers:
// where does a fitted seed next cross a barrel cylinder, an endcap disk, or one sensor plane?
//
// The helix is defined GEOMETRICALLY, without charge or field sign: a circle (cx, cy, R),
// the rotation sense `rot` of travel about its centre (taken from the hit order), and
// z = zRef + tanL * s where s is the transverse arc length from a reference point. That
// keeps it immune to the charge/field-sign conventions that are already duplicated three
// ways in the component.
//
// Uniform-field model (a helix). Units: cm. No Gaudi/DD4hep dependency, so it is
// unit-testable on its own (Tracking/test/testCombinatorialTrackFinder/test_HelixExtrapolator.cpp).
#ifndef TRACKING_HELIXEXTRAPOLATOR_H
#define TRACKING_HELIXEXTRAPOLATOR_H

#include "LayerGeometry.h"

#include <Eigen/Dense>

#include <limits>
#include <vector>

namespace k4rt {

struct Helix {
  double cx = 0, cy = 0, R = 0; ///< transverse circle, cm
  int rot = +1;                 ///< +1: angle about the centre increases along travel
  double phiRef = 0;            ///< angle of the reference point about the centre
  double zRef = 0;              ///< z at the reference point, cm
  double tanL = 0;              ///< dz / ds (s = transverse arc length)

  double angleAt(double s) const { return phiRef + rot * s / R; }
  Eigen::Vector3d at(double s) const;    ///< position at arc length s from the reference
  Eigen::Vector3d dirAt(double s) const; ///< unit 3D direction of travel at s
};

/// Build a helix from a fitted circle and hits ordered ALONG THE TRACK (cm). The reference
/// point is the first hit; rot comes from the hit order; zRef/tanL from a linear z-vs-s fit
/// with s unwrapped hit-to-hit (so tracks turning by more than pi are handled).
/// Returns false if fewer than 2 hits or R <= 0.
bool helixFromHits(double cx, double cy, double R, const std::vector<Eigen::Vector3d>& orderedHitsCm, Helix& out,
                   std::vector<double>* sOfHits = nullptr);

/// Arc length of the point on the helix closest (in angle) to p, unwrapped to lie nearest
/// sGuess. Used to place a hit on the helix.
double arcLengthNear(const Helix& h, const Eigen::Vector3d& p, double sGuess);

struct Crossing {
  bool ok = false;
  double s = 0; ///< arc length at the crossing
  Eigen::Vector3d pos = Eigen::Vector3d::Zero();
};

/// First crossing of the transverse radius Rc, moving from sFrom in direction dir (+1/-1),
/// at least sEps and at most sMax away (in |s|).
Crossing crossCylinder(const Helix& h, double Rc, double sFrom, int dir, double sEps = 1e-6,
                       double sMax = std::numeric_limits<double>::infinity());

/// Crossing of the plane z = zPlane (unique: z is monotonic in s). Same s window.
Crossing crossPlaneZ(const Helix& h, double zPlane, double sFrom, int dir, double sEps = 1e-6,
                     double sMax = std::numeric_limits<double>::infinity());

/// Crossing of a layer within its extents (+- tol): cylinder at rMean for a barrel layer,
/// plane at zMean for a disk. Approximate by construction - a barrel layer is a polygon of
/// flat staves - so use refineToPlane() against the actual sensor for residuals.
Crossing crossLayer(const Helix& h, const LayerDescriptor& L, double sFrom, int dir, double tol_cm,
                    double sMax = std::numeric_limits<double>::infinity());

/// Exact intersection with the plane through `origin` with unit normal `normal`, by Newton
/// iteration starting from sStart. ok=false if it does not converge.
Crossing refineToPlane(const Helix& h, const Eigen::Vector3d& origin, const Eigen::Vector3d& normal, double sStart,
                       int maxIter = 20, double tol_cm = 1e-9);

} // namespace k4rt

#endif // TRACKING_HELIXEXTRAPOLATOR_H
