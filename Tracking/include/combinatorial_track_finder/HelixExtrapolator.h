/**
 * HelixExtrapolator: helix model and layer crossings of CombinatorialTrackFinder.
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

// The class Helix answers the seed-and-extend question: where does a fitted track next cross a
// barrel cylinder, an endcap disk or a sensor plane, going forward or backward?
//
// The helix is GEOMETRIC (no charge or field sign): a transverse circle (centre, radius), the
// rotation sense of travel about the centre (+1: the angle grows along the track, from the hit
// order) and a line z = referenceZ + tanLambda * s, with s the transverse arc length from the
// first hit. The point at arc length s is at the angle referenceAngle + rotation * s / radius.
// Uniform field, units cm; no Gaudi/DD4hep (unit test test_HelixExtrapolator.cpp).
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
  double arcLength = 0; ///< cm
  Eigen::Vector3d position = Eigen::Vector3d::Zero();
};

class Helix {
public:
  double centreX = 0, centreY = 0, radius = 0; ///< transverse circle, cm
  int rotation = +1;                           ///< +1: the angle about the centre grows along the travel
  double referenceAngle = 0;                   ///< angle of the reference point about the centre
  double referenceZ = 0;                       ///< z at the reference point, cm
  double tanLambda = 0;                        ///< dz / ds (s = transverse arc length)

  /// Helix from a fitted circle and hits ordered ALONG THE TRACK (cm): reference = first hit,
  /// rotation from the hit order, z line fitted versus s (s unwrapped hit to hit, so turns beyond
  /// pi work). The arc length of every hit goes to arcLengthsOfHits if given. false with fewer
  /// than 2 hits or radius <= 0.
  static bool fromHits(double circleCentreX, double circleCentreY, double circleRadius,
                       const std::vector<Eigen::Vector3d>& orderedHitsCm, Helix& out,
                       std::vector<double>* arcLengthsOfHits = nullptr);

  double angleAt(double arcLength) const { return referenceAngle + rotation * arcLength / radius; }
  Eigen::Vector3d positionAt(double arcLength) const;
  Eigen::Vector3d directionAt(double arcLength) const; ///< unit 3D direction of travel

  /// Arc length of the helix point closest in angle to the point, the one nearest to guessArcLength.
  double arcLengthNear(const Eigen::Vector3d& point, double guessArcLength) const;

  /// First crossing of the transverse radius cylinderRadius from fromArcLength in the direction
  /// (+1 forward, -1 backward), between minArcLength and maxArcLength away.
  Crossing crossCylinder(double cylinderRadius, double fromArcLength, int direction, double minArcLength = 1e-6,
                         double maxArcLength = std::numeric_limits<double>::infinity()) const;

  /// Crossing of the plane z = planeZ (unique: z is monotonic in s); same arc-length window.
  Crossing crossPlaneZ(double planeZ, double fromArcLength, int direction, double minArcLength = 1e-6,
                       double maxArcLength = std::numeric_limits<double>::infinity()) const;

  /// Crossing of a layer within its extent (+- tolerance): the cylinder at its mean radius (barrel)
  /// or the plane at its mean z (disk). Approximate (staves are flat): refine with refineToPlane().
  Crossing crossLayer(const LayerDescriptor& layer, double fromArcLength, int direction, double toleranceCm,
                      double maxArcLength = std::numeric_limits<double>::infinity()) const;

  /// Exact crossing of a sensor plane (point, unit normal), Newton iteration from startArcLength;
  /// ok = false if it does not converge.
  Crossing refineToPlane(const Eigen::Vector3d& planeOrigin, const Eigen::Vector3d& planeNormal, double startArcLength,
                         int maxIterations = 20, double toleranceCm = 1e-9) const;
};

} // namespace k4rt

#endif // TRACKING_HELIXEXTRAPOLATOR_H
