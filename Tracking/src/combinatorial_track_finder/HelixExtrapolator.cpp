/**
 * HelixExtrapolator: helix model and layer crossings of CombinatorialTrackFinder.
 * See HelixExtrapolator.h for an overview.
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

#include "HelixExtrapolator.h"
#include "TrackingMath.h"

#include <algorithm>
#include <cmath>

namespace k4rt {

namespace {
  constexpr double kTwoPi = 2.0 * M_PI;
  // angle in (0, 2pi]: unlike wrapAngleZeroToTwoPi, 0 becomes a full turn (a crossing at the
  // starting point itself is the next turn)
  double wrapPositive(double angle) {
    angle = std::fmod(angle, kTwoPi);
    if (angle <= 0.0)
      angle += kTwoPi;
    return angle;
  }
} // namespace

Eigen::Vector3d Helix::positionAt(double arcLength) const {
  const double angle = angleAt(arcLength);
  return {centreX + radius * std::cos(angle), centreY + radius * std::sin(angle), referenceZ + tanLambda * arcLength};
}

Eigen::Vector3d Helix::directionAt(double arcLength) const {
  // derivative of positionAt: the unit transverse tangent (turned by the rotation) and tanLambda
  const double angle = angleAt(arcLength);
  Eigen::Vector3d direction(-rotation * std::sin(angle), rotation * std::cos(angle), tanLambda);
  return direction.normalized();
}

bool Helix::fromHits(double circleCentreX, double circleCentreY, double circleRadius,
                     const std::vector<Eigen::Vector3d>& hits, Helix& out, std::vector<double>* arcLengthsOfHits) {
  const int numberOfHits = static_cast<int>(hits.size());
  if (numberOfHits < 2 || !(circleRadius > 0.0))
    return false;
  out.centreX = circleCentreX;
  out.centreY = circleCentreY;
  out.radius = circleRadius;

  // 1. angle of every hit about the centre, unwrapped (each step the smallest turn from the last)
  std::vector<double> angles(numberOfHits);
  angles[0] = std::atan2(hits[0].y() - circleCentreY, hits[0].x() - circleCentreX);
  for (int i = 1; i < numberOfHits; ++i) {
    const double angle = std::atan2(hits[i].y() - circleCentreY, hits[i].x() - circleCentreX);
    angles[i] = angles[i - 1] + wrapAngle(angle - angles[i - 1]);
  }
  // 2. rotation sense: the sign of the total turn
  out.rotation = (angles[numberOfHits - 1] - angles[0] >= 0.0) ? +1 : -1;
  out.referenceAngle = angles[0];

  // 3. arc lengths s_i = rotation * R * (angle_i - angle_0) and the least-squares line z = z0 + t s:
  //    t = (n sum(sz) - sum(s) sum(z)) / det, z0 = (sum(z) sum(s^2) - sum(s) sum(sz)) / det,
  //    det = n sum(s^2) - sum(s)^2
  std::vector<double> arcLengths(numberOfHits);
  double sumS = 0, sumZ = 0, sumSS = 0, sumSZ = 0;
  for (int i = 0; i < numberOfHits; ++i) {
    arcLengths[i] = out.rotation * circleRadius * (angles[i] - angles[0]);
    sumS += arcLengths[i];
    sumZ += hits[i].z();
    sumSS += arcLengths[i] * arcLengths[i];
    sumSZ += arcLengths[i] * hits[i].z();
  }
  const double determinant = numberOfHits * sumSS - sumS * sumS;
  if (std::abs(determinant) > 1e-12) {
    out.referenceZ = (sumZ * sumSS - sumS * sumSZ) / determinant;
    out.tanLambda = (numberOfHits * sumSZ - sumS * sumZ) / determinant;
  } else {
    out.referenceZ = sumZ / numberOfHits; // no lever arm: flat line at the mean z
    out.tanLambda = 0.0;
  }
  if (arcLengthsOfHits)
    *arcLengthsOfHits = arcLengths;
  return true;
}

double Helix::arcLengthNear(const Eigen::Vector3d& point, double guessArcLength) const {
  // the smallest angular move from the guess (wrapped to (-pi, pi]), turned into an arc length
  const double angle = std::atan2(point.y() - centreY, point.x() - centreX);
  const double angleAtGuess = angleAt(guessArcLength);
  return guessArcLength + rotation * radius * wrapAngle(angle - angleAtGuess);
}

Crossing Helix::crossCylinder(double cylinderRadius, double fromArcLength, int direction, double minArcLength,
                              double maxArcLength) const {
  Crossing best;
  // Two circles (track: centre c, radius R; cylinder: radius Rc at the origin), d = |c|: they
  // cross unless d > R + Rc or d < |R - Rc|; the crossings lie on the line along c at
  // a = (Rc^2 - R^2 + d^2) / (2 d) from the origin, +-h = sqrt(Rc^2 - a^2) sideways.
  const double centreDistance = std::hypot(centreX, centreY);
  if (centreDistance < 1e-12)
    return best; // centred on the axis: the radius never changes
  if (centreDistance > radius + cylinderRadius || centreDistance < std::abs(radius - cylinderRadius))
    return best; // never reaches the cylinder
  const double footDistance =
      (cylinderRadius * cylinderRadius - radius * radius + centreDistance * centreDistance) / (2.0 * centreDistance);
  const double halfChordSquared = cylinderRadius * cylinderRadius - footDistance * footDistance;
  const double halfChord = std::sqrt(std::max(0.0, halfChordSquared));
  const double unitX = centreX / centreDistance; // unit vector from the origin to the centre
  const double unitY = centreY / centreDistance;

  // the crossing reached first from fromArcLength in the direction: moving by ds turns the angle
  // by rotation * direction * ds / R, so the turn still needed (in (0, 2pi]) gives ds = R * turn
  const double angleFrom = angleAt(fromArcLength);
  double bestDistance = std::numeric_limits<double>::infinity();
  for (int side = -1; side <= 1; side += 2) {
    const double crossingX = footDistance * unitX - side * halfChord * unitY;
    const double crossingY = footDistance * unitY + side * halfChord * unitX;
    const double crossingAngle = std::atan2(crossingY - centreY, crossingX - centreX);
    const double turn = wrapPositive(rotation * direction * (crossingAngle - angleFrom));
    const double distance = radius * turn;
    if (distance < minArcLength || distance > maxArcLength)
      continue;
    if (distance < bestDistance)
      bestDistance = distance;
  }
  if (!std::isfinite(bestDistance))
    return best;
  best.ok = true;
  best.arcLength = fromArcLength + direction * bestDistance;
  best.position = positionAt(best.arcLength);
  return best;
}

Crossing Helix::crossPlaneZ(double planeZ, double fromArcLength, int direction, double minArcLength,
                            double maxArcLength) const {
  Crossing crossing;
  if (std::abs(tanLambda) < 1e-12)
    return crossing; // parallel to the plane
  const double arcLength = (planeZ - referenceZ) / tanLambda;
  const double distance = direction * (arcLength - fromArcLength); // travelled in the direction
  if (distance < minArcLength || distance > maxArcLength)
    return crossing;
  crossing.ok = true;
  crossing.arcLength = arcLength;
  crossing.position = positionAt(arcLength);
  return crossing;
}

Crossing Helix::crossLayer(const LayerDescriptor& layer, double fromArcLength, int direction, double tolerance,
                           double maxArcLength) const {
  Crossing crossing;
  if (layer.numberOfSurfaces == 0)
    return crossing;
  if (layer.isBarrel()) {
    // the cylinder at the mean radius, inside the z extent (+- tolerance)
    crossing = crossCylinder(layer.rMean_cm, fromArcLength, direction, 1e-6, maxArcLength);
    const double z = crossing.position.z();
    if (crossing.ok && (z < layer.zMin_cm - tolerance || z > layer.zMax_cm + tolerance))
      crossing.ok = false;
  } else {
    // the plane at the mean z, inside the radial extent (+- tolerance)
    crossing = crossPlaneZ(layer.zMean_cm, fromArcLength, direction, 1e-6, maxArcLength);
    if (crossing.ok) {
      const double transverseRadius = std::hypot(crossing.position.x(), crossing.position.y());
      if (transverseRadius < layer.rMin_cm - tolerance || transverseRadius > layer.rMax_cm + tolerance)
        crossing.ok = false;
    }
  }
  return crossing;
}

Crossing Helix::refineToPlane(const Eigen::Vector3d& planeOrigin, const Eigen::Vector3d& planeNormal,
                              double startArcLength, int maxIterations, double tolerance) const {
  Crossing crossing;
  // Newton on f(s) = n . (positionAt(s) - o): s <- s - f / f', f' = n . d(position)/ds
  double arcLength = startArcLength;
  for (int iteration = 0; iteration < maxIterations; ++iteration) {
    const Eigen::Vector3d point = positionAt(arcLength);
    const double distance = planeNormal.dot(point - planeOrigin);
    if (std::abs(distance) < tolerance) {
      crossing.ok = true;
      crossing.arcLength = arcLength;
      crossing.position = point;
      return crossing;
    }
    // d(position)/ds, NOT normalised: unit transverse part, dz/ds = tanLambda
    const double angle = angleAt(arcLength);
    const Eigen::Vector3d derivative(-rotation * std::sin(angle), rotation * std::cos(angle), tanLambda);
    const double slope = planeNormal.dot(derivative);
    if (std::abs(slope) < 1e-12)
      return crossing; // parallel to the plane
    arcLength -= distance / slope;
  }
  return crossing;
}

} // namespace k4rt
