// The kink chi2 and the Highland angle (see TrackingMath.h).
#include "TrackingMath.h"

#include <algorithm>
#include <cmath>

namespace k4rt {

namespace {
  // d/dR of asin(l / 2R), the angle between a chord of length l and the circle:
  // -l / (2 R^2 sqrt(1 - (l/2R)^2)), with l/2R kept below 1
  double chordAngleDerivative(double chordLength, double radius) {
    const double halfChordOverRadius = std::min(0.999999, chordLength / (2.0 * radius));
    return -chordLength / (2.0 * radius * radius * std::sqrt(1.0 - halfChordOverRadius * halfChordOverRadius));
  }
} // namespace

double KinkChi2Calculator::highlandTheta0(double radiationLengths, double momentumGeV) {
  if (!(radiationLengths > 0.0) || !(momentumGeV > 0.0))
    return 0.0;
  // Highland: theta0 = 13.6 MeV / p * sqrt(x/X0) * (1 + 0.038 ln(x/X0)) for beta = 1; the log
  // factor is kept above 0.1, so very thin layers still give a small positive angle
  const double logTerm = 1.0 + 0.038 * std::log(radiationLengths);
  return 0.0136 / momentumGeV * std::sqrt(radiationLengths) * std::max(logTerm, 0.1);
}

double KinkChi2Calculator::triplet(const KinkPoint& first, const KinkPoint& middle, const KinkPoint& last,
                                   double radius, int rotation, double sigmaRadius) const {
  // the two chords first->middle and middle->last in the transverse plane
  const Eigen::Vector2d chord1(middle.position.x() - first.position.x(), middle.position.y() - first.position.y());
  const Eigen::Vector2d chord2(last.position.x() - middle.position.x(), last.position.y() - middle.position.y());
  const double length1 = chord1.norm(), length2 = chord2.norm();
  if (length1 < 1e-6 || length2 < 1e-6 || !(radius > 0.0))
    return 0.0;

  // 1. bending-plane kink: measured turn between the chords minus the helix turn; a chord of
  //    length l makes the angle asin(l / 2R) with the tangent at each end
  const double measuredTurn = std::atan2(chord1.x() * chord2.y() - chord1.y() * chord2.x(), chord1.dot(chord2));
  const double chordAngle1 = std::asin(std::min(1.0, length1 / (2.0 * radius)));
  const double chordAngle2 = std::asin(std::min(1.0, length2 / (2.0 * radius)));
  const double kinkPhi = wrapAngle(measuredTurn - rotation * (chordAngle1 + chordAngle2));

  // 2. dip kink: change of dz/ds between the segments (s = arc length 2 R angle, not the chord)
  const double arc1 = 2.0 * radius * chordAngle1, arc2 = 2.0 * radius * chordAngle2;
  const double dz1 = middle.position.z() - first.position.z(), dz2 = last.position.z() - middle.position.z();
  const double kinkLambda = std::atan2(dz2, arc2) - std::atan2(dz1, arc1);

  // 3. expected scattering at the middle hit: p = pT / cos(lambda); the material grows as
  //    1 / cos(incidence) on the sensor (limited to 1 / 0.05 for grazing tracks)
  const double tanLambda = (dz1 + dz2) / (arc1 + arc2);
  const double cosLambda = 1.0 / std::sqrt(1.0 + tanLambda * tanLambda);
  const double pT = ptFromRadius(m_bFieldTesla, radius); // GeV
  const double momentum = pT / cosLambda;
  const Eigen::Vector3d direction = (last.position - first.position).normalized();
  const double cosIncidence = std::max(0.05, std::abs(middle.normal.normalized().dot(direction)));
  const double theta0 = highlandTheta0(middle.xOverX0 / cosIncidence, momentum);

  // 4. variance of the bending-plane kink: radius error (turn derivative x sigmaR), scattering
  //    projected on the transverse plane (theta0 / cos(lambda)), and measurement errors (a shift
  //    sigma turns a chord by sigma / l; the middle hit moves both chords)
  const double inverseLength1 = 1.0 / length1, inverseLength2 = 1.0 / length2;
  const double sigmaTurnFromRadius =
      (chordAngleDerivative(length1, radius) + chordAngleDerivative(length2, radius)) * sigmaRadius;
  const double variancePhi = sigmaTurnFromRadius * sigmaTurnFromRadius + (theta0 / cosLambda) * (theta0 / cosLambda) +
                             first.sigmaTransverseCm * first.sigmaTransverseCm * inverseLength1 * inverseLength1 +
                             middle.sigmaTransverseCm * middle.sigmaTransverseCm * (inverseLength1 + inverseLength2) *
                                 (inverseLength1 + inverseLength2) +
                             last.sigmaTransverseCm * last.sigmaTransverseCm * inverseLength2 * inverseLength2;
  // 5. variance of the dip kink: scattering plus z errors as slope errors (sigma_z / s), times
  //    cos^4(lambda) because d(lambda) = cos^2(lambda) d(tan lambda)
  const double inverseArc1 = 1.0 / arc1, inverseArc2 = 1.0 / arc2;
  const double cosLambdaToTheFourth = cosLambda * cosLambda * cosLambda * cosLambda;
  const double varianceLambda =
      theta0 * theta0 + cosLambdaToTheFourth * (first.sigmaZCm * first.sigmaZCm * inverseArc1 * inverseArc1 +
                                                middle.sigmaZCm * middle.sigmaZCm * (inverseArc1 + inverseArc2) *
                                                    (inverseArc1 + inverseArc2) +
                                                last.sigmaZCm * last.sigmaZCm * inverseArc2 * inverseArc2);
  // 6. chi2 of the two kinks (2 dof)
  return kinkPhi * kinkPhi / variancePhi + kinkLambda * kinkLambda / varianceLambda;
}

KinkResult KinkChi2Calculator::track(const std::vector<KinkPoint>& points, double radius, int rotation) const {
  KinkResult result;
  const int numberOfPoints = static_cast<int>(points.size());
  result.perTriplet.assign(std::max(numberOfPoints, 0), 0.0);
  if (numberOfPoints < 3)
    return result;
  // one kink at every inner hit; independent kinks, so their chi2 add up
  for (int i = 1; i + 1 < numberOfPoints; ++i) {
    const double tripletChi2 = triplet(points[i - 1], points[i], points[i + 1], radius, rotation);
    result.perTriplet[i] = tripletChi2;
    result.chi2 += tripletChi2;
  }
  // 2 angles per kink, minus 1 for the radius taken from the fit
  result.ndf = 2 * (numberOfPoints - 2) - 1;
  return result;
}

} // namespace k4rt
