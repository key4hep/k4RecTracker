// The hit pool, the circle fit and the EDM4hep track states (see TrackCandidate.h).
#include "TrackCandidate.h"
#include "LayerGeometry.h"
#include "TrackingMath.h"

#include <algorithm>
#include <cmath>

namespace CTF {

namespace {
  /// A hit as the state builder sees it (mm).
  struct StateHitPoint {
    double x, y, z;
    double distance;  ///< from the origin
    double varianceZ; ///< mm^2
  };

  /// The weighted line z = zPca + tanLambda * s of a helix state, with its covariance terms.
  struct ZLineFit {
    double tanLambda = 0.0;
    double zPca = 0.0;
    double sumW = 0.0;     ///< sum of the weights
    double meanS = 0.0;    ///< weighted mean arc length
    double sumWDSDS = 0.0; ///< sum w (s - meanS)^2
  };

  // helpers, defined at the end of the file
  void circleResidualsAndJacobian(const std::vector<Eigen::Vector2d>& points, const std::vector<double>& weights,
                                  const Eigen::Vector3d& circle, Eigen::MatrixXd& jacobian, Eigen::VectorXd& residuals);
  std::vector<StateHitPoint> hitsByDistanceFromOrigin(const std::vector<EventHit>& hits,
                                                      const std::vector<size_t>& trackHits);
  double rotationSense(const std::vector<StateHitPoint>& points, double centreX, double centreY);
  std::vector<double> arcLengthsFromPca(const std::vector<StateHitPoint>& points, double centreX, double centreY,
                                        double radius, double sense, double anglePca, size_t nearestHit);
  bool fitZAlongArc(const std::vector<StateHitPoint>& points, const std::vector<double>& arcLength, ZLineFit& line);
  void setStateCovariance(const CircleFit& circle, double centreToReferenceX, double centreToReferenceY,
                          double centreToReference, double omega, const ZLineFit& line, edm4hep::TrackState& state);
  bool solenoidBarrelCrossing(const OuterHelix& outer, double solenoidRadiusMm, bool tryOtherCrossing,
                              Eigen::Vector3d& crossingCm);
  void circleCylinderCrossings(double centreX, double centreY, double radius, double centreDistance,
                               double solenoidRadius, double crossingX[2], double crossingY[2]);
  double backwardArcAngle(double centreX, double centreY, double angleOfHit, bool clockwise, double pointX,
                          double pointY);
} // namespace

// ---- the free hits of one iteration ----

bool collectHitPool(const std::vector<EventHit>& hits, const std::vector<bool>& used, double maxHitEdepKeV,
                    HitPool& pool) {
  int remainingHits = 0;
  for (size_t i = 0; i < used.size(); ++i) {
    if (!used[i])
      remainingHits++;
  }
  if (remainingHits < 3)
    return false;

  double edepThresholdGeV = -1.0; // MaxHitEdepKeV (EDM4hep edep in GeV); -1 = off
  if (maxHitEdepKeV > 0.0)
    edepThresholdGeV = maxHitEdepKeV * 1e-6;
  pool.hits.clear();
  pool.hitsByLayer.clear();
  for (size_t i = 0; i < hits.size(); ++i) {
    if (used[i] || !hits[i].usable)
      continue;
    if (edepThresholdGeV > 0.0 && hits[i].edepGeV > edepThresholdGeV)
      continue;
    pool.hits.push_back(i);
    pool.hitsByLayer[hits[i].layerKey].push_back(i);
  }
  if (pool.hits.size() < 3)
    return false;

  // layer keys by layer number (innermost first), then by key
  pool.layerOrder.clear();
  for (const auto& [layerKey, layerHits] : pool.hitsByLayer)
    pool.layerOrder.push_back(layerKey);
  std::sort(pool.layerOrder.begin(), pool.layerOrder.end(), [](int left, int right) {
    const int leftLayer = keyWithoutRegion(left);
    const int rightLayer = keyWithoutRegion(right);
    if (leftLayer != rightLayer)
      return leftLayer < rightLayer;
    return left < right;
  });
  return pool.layerOrder.size() >= 3;
}

// Bz (T) at the centroid of the hits
double fieldAtCentroid(const std::vector<EventHit>& hits, const std::vector<size_t>& hitIndices,
                       const FieldFunction& bzAt) {
  const int numberOfHits = static_cast<int>(hitIndices.size());
  double sumX = 0.0, sumY = 0.0, sumZ = 0.0; // cm
  for (size_t hitIndex : hitIndices) {
    const auto& position = hits[hitIndex].positionMm;
    sumX += position[0] / kMmPerCm;
    sumY += position[1] / kMmPerCm;
    sumZ += position[2] / kMmPerCm;
  }
  return bzAt(sumX / numberOfHits, sumY / numberOfHits, sumZ / numberOfHits);
}

// ---- the circle fit ----

std::optional<CircleFit> fitCircle(const std::vector<EventHit>& hits, const std::vector<size_t>& hitIndices,
                                   double defaultSigmaCm) {
  const int numberOfHits = static_cast<int>(hitIndices.size());
  if (numberOfHits < 3)
    return std::nullopt;

  // positions (cm) and weights 1 / du^2 (one radial residual per hit)
  std::vector<Eigen::Vector2d> points(numberOfHits);
  std::vector<double> sigmas(numberOfHits);
  std::vector<double> weights(numberOfHits);
  for (int i = 0; i < numberOfHits; ++i) {
    const EventHit& hit = hits[hitIndices[i]];
    const auto& position = hit.positionMm;
    points[i] = {position[0] / kMmPerCm, position[1] / kMmPerCm};
    double du = hit.duMm; // mm
    if (du <= 0)
      du = defaultSigmaCm * kMmPerCm;
    sigmas[i] = du / kMmPerCm;
    weights[i] = 1.0 / (sigmas[i] * sigmas[i]);
  }

  // Step 1, algebraic first guess: x^2 + y^2 + D x + E y + F = 0 is linear in (D, E, F); the
  // rows (scaled by w, so a hit counts with w^2) are solved by least squares, (A^T A) p = A^T b;
  // then centre (-D/2, -E/2), radius sqrt(D^2 + E^2 - 4F) / 2
  Eigen::MatrixXd designMatrix(numberOfHits, 3);
  Eigen::VectorXd rightHandSide(numberOfHits);
  for (int i = 0; i < numberOfHits; ++i) {
    const double x = points[i].x();
    const double y = points[i].y();
    designMatrix(i, 0) = weights[i] * x;
    designMatrix(i, 1) = weights[i] * y;
    designMatrix(i, 2) = weights[i] * 1.0;
    rightHandSide(i) = -weights[i] * (x * x + y * y);
  }
  const Eigen::Vector3d solution =
      (designMatrix.transpose() * designMatrix).ldlt().solve(designMatrix.transpose() * rightHandSide);
  const double coefficientD = solution(0);
  const double coefficientE = solution(1);
  const double coefficientF = solution(2);
  const double discriminant = coefficientD * coefficientD + coefficientE * coefficientE - 4 * coefficientF;
  if (discriminant <= 0)
    return std::nullopt; // no real circle
  Eigen::Vector3d circle(-coefficientD / 2.0, -coefficientE / 2.0, std::sqrt(discriminant) / 2.0);

  // Step 2, Gauss-Newton on the real hit distances: delta = -(J^T J)^-1 J^T r
  Eigen::MatrixXd jacobian(numberOfHits, 3);
  Eigen::VectorXd residuals(numberOfHits);
  for (int step = 0; step < 10; ++step) { // at most 10 Gauss-Newton steps
    circleResidualsAndJacobian(points, weights, circle, jacobian, residuals);
    const Eigen::Vector3d delta = -(jacobian.transpose() * jacobian).ldlt().solve(jacobian.transpose() * residuals);
    circle += delta;
    if (delta.norm() < 1e-6) // converged: the step moved (cx, cy, R) by less than 1e-6 cm
      break;
  }
  CircleFit fit;
  fit.centreX = circle(0);
  fit.centreY = circle(1);
  fit.radius = circle(2);
  if (fit.radius <= 0)
    return std::nullopt;

  // Step 3, covariance of (cx, cy, R) = (J^T J)^-1 at the final circle
  circleResidualsAndJacobian(points, weights, circle, jacobian, residuals);
  fit.covariance = (jacobian.transpose() * jacobian).ldlt().solve(Eigen::Matrix3d::Identity());

  // Step 4, residuals and chi2/ndf (ndf = N - 3, at least 1)
  fit.residuals.resize(numberOfHits);
  double chi2 = 0.0;
  for (int i = 0; i < numberOfHits; ++i) {
    const double dx = points[i].x() - fit.centreX;
    const double dy = points[i].y() - fit.centreY;
    const double residual = std::abs(std::sqrt(dx * dx + dy * dy) - fit.radius);
    fit.residuals[i] = residual;
    chi2 += (residual / sigmas[i]) * (residual / sigmas[i]);
  }
  int ndf = numberOfHits - 3;
  if (ndf < 1)
    ndf = 1;
  fit.chi2ndf = chi2 / ndf;
  return fit;
}

// ---- EDM4hep track states ----

// The fitted helix (circle + straight line of z versus transverse arc length) as an EDM4hep track
// state. Geometric, so valid for either field sign: the particle moves outward (from the hit
// nearest the origin), omega > 0 = clockwise seen from +z, and d0, phi, z0 are given at the circle
// point closest to the reference point (innermost hit, outermost hit or origin), EDM4hep
// convention: pca = ref + (-d0 sin(phi), d0 cos(phi), z0), centre = pca + (1/omega)(sin(phi), -cos(phi)).
// Circle in cm, state in mm. false (state unchanged) for degenerate input.
bool makeStateFromCircle(const CircleFit& circle, const std::vector<EventHit>& hits,
                         const std::vector<size_t>& trackHits, StatePoint where, edm4hep::TrackState& state) {
  const double centreX = circle.centreX * kMmPerCm; // mm
  const double centreY = circle.centreY * kMmPerCm;
  const double radius = circle.radius * kMmPerCm;
  if (radius <= 0.0 || trackHits.size() < 2)
    return false;

  // 1. hits sorted by distance from the origin (the order an outward-moving particle crosses them)
  const std::vector<StateHitPoint> points = hitsByDistanceFromOrigin(hits, trackHits);

  // 2. sense of rotation seen from +z (+1 = counter-clockwise)
  const double sense = rotationSense(points, centreX, centreY);
  if (sense == 0.0)
    return false;

  // 3. reference point and the closest circle point (pca)
  double referenceX = 0.0, referenceY = 0.0, referenceZ = 0.0;
  size_t nearestHit = 0; // the hit next to the reference point, to count the arc
  if (where == StatePoint::FirstHit) {
    referenceX = points.front().x;
    referenceY = points.front().y;
    referenceZ = points.front().z;
  } else if (where == StatePoint::LastHit) {
    nearestHit = points.size() - 1;
    referenceX = points.back().x;
    referenceY = points.back().y;
    referenceZ = points.back().z;
  }
  const double centreToReferenceX = referenceX - centreX;
  const double centreToReferenceY = referenceY - centreY;
  const double centreToReference =
      std::sqrt(centreToReferenceX * centreToReferenceX + centreToReferenceY * centreToReferenceY);
  if (centreToReference < 1e-6)
    return false; // reference at the centre: no closest point
  const double anglePca = std::atan2(centreToReferenceY, centreToReferenceX);
  const double pcaX = centreX + radius * centreToReferenceX / centreToReference;
  const double pcaY = centreY + radius * centreToReferenceY / centreToReference;

  // 4. phi = direction of travel at the pca (radius vector turned by +-90 degrees), omega > 0 for
  //    clockwise motion, d0 from pca - ref = (-d0 sin(phi), d0 cos(phi))
  const double phi = wrapAngle(anglePca + sense * M_PI / 2.0);
  const double omega = -sense / radius; // 1/mm
  const double d0 = -(pcaX - referenceX) * std::sin(phi) + (pcaY - referenceY) * std::cos(phi);

  // 5. the arc length of every hit from the pca
  const std::vector<double> arcLength =
      arcLengthsFromPca(points, centreX, centreY, radius, sense, anglePca, nearestHit);

  // 6. the line z = zPca + tanLambda * s
  ZLineFit line;
  if (!fitZAlongArc(points, arcLength, line))
    return false;

  // 7. covariance, then the parameters
  setStateCovariance(circle, centreToReferenceX, centreToReferenceY, centreToReference, omega, line, state);
  state.D0 = static_cast<float>(d0);
  state.phi = static_cast<float>(phi);
  state.omega = static_cast<float>(omega);
  state.Z0 = static_cast<float>(line.zPca - referenceZ);
  state.tanLambda = static_cast<float>(line.tanLambda);
  state.referencePoint = edm4hep::Vector3f(referenceX, referenceY, referenceZ);
  return true;
}

// track state with these parameters and a default diagonal covariance (D0, Z0 in mm^2,
// phi, tanLambda in rad^2, omega in (1/mm)^2)
edm4hep::TrackState stateWithDefaultCovariance(double d0, double phi, double omega, double z0, double tanLambda,
                                               int location) {
  edm4hep::TrackState state;
  state.D0 = d0;
  state.phi = phi;
  state.omega = omega;
  state.Z0 = z0;
  state.tanLambda = tanLambda;
  state.location = location;
  for (int row = 0; row < 6; ++row)
    for (int column = 0; column <= row; ++column)
      state.setCovMatrix(0.0, static_cast<edm4hep::TrackParams>(row), static_cast<edm4hep::TrackParams>(column));
  // the default diagonal: d0 and z0 1 mm^2, phi and tanLambda 0.01, omega 1e-10 (1/mm)^2
  state.setCovMatrix(1.0, edm4hep::TrackParams::d0, edm4hep::TrackParams::d0);
  state.setCovMatrix(0.01, edm4hep::TrackParams::phi, edm4hep::TrackParams::phi);
  state.setCovMatrix(1e-10, edm4hep::TrackParams::omega, edm4hep::TrackParams::omega);
  state.setCovMatrix(1.0, edm4hep::TrackParams::z0, edm4hep::TrackParams::z0);
  state.setCovMatrix(0.01, edm4hep::TrackParams::tanLambda, edm4hep::TrackParams::tanLambda);
  return state;
}

// ---- the SolenoidCrossing AtIP state ----

bool Solenoid::contains(const Eigen::Vector3d& positionCm) const {
  const double transverseRadius = std::sqrt(positionCm.x() * positionCm.x() + positionCm.y() * positionCm.y());
  const double absoluteZ = std::abs(positionCm.z());
  return (transverseRadius < radiusMm / kMmPerCm) && (absoluteZ < halfLengthMm / kMmPerCm);
}

// The solenoid boundary = a cylinder at R_sol closed by disks at |z| = Z_sol; the outer circle is
// followed back from the innermost hit to it, then the inner helix (inner field, same pT and
// tanLambda) is attached and its perigee returned.
std::optional<edm4hep::TrackState> solenoidCrossingState(const OuterHelix& outer, const Solenoid& solenoid) {
  const double solenoidHalfLengthCm = solenoid.halfLengthMm / kMmPerCm;
  const double centreX = outer.centreXCm;
  const double centreY = outer.centreYCm;
  const double radius = outer.radiusCm;
  // backward on the outer circle (+1: angle grows, clockwise track); also the reversed curvature of
  // the inner helix
  const double innerSign = outer.clockwise ? 1.0 : -1.0;

  // 1. the crossing of the boundary (cm)
  Eigen::Vector3d crossing;
  if (!outer.innermostHitInEndcap) {
    if (!solenoidBarrelCrossing(outer, solenoid.radiusMm, true, crossing))
      return std::nullopt;
  } else {
    // endcap hit: the end disk, at the arc s back from the hit with z(s) = +-Z_sol
    if (std::abs(outer.tanLambda) < 1e-6)
      return std::nullopt;
    const double hitZ = outer.innermostHitCm.z();
    const double diskZ = (hitZ < 0.0) ? -solenoidHalfLengthCm : +solenoidHalfLengthCm;
    const double arcToDisk = (hitZ - diskZ) / outer.tanLambda; // hitZ - tanLambda * s = diskZ
    if (arcToDisk < 0.0) {
      // the hit is already at |z| < Z_sol: the barrel crossing instead
      if (!solenoidBarrelCrossing(outer, solenoid.radiusMm, false, crossing))
        return std::nullopt;
    } else {
      if (arcToDisk > 2.0 * M_PI * radius)
        return std::nullopt;
      const double angleOfHit = std::atan2(outer.innermostHitCm.y() - centreY, outer.innermostHitCm.x() - centreX);
      const double angleAtDisk = angleOfHit + innerSign * (arcToDisk / radius);
      crossing =
          Eigen::Vector3d(centreX + radius * std::cos(angleAtDisk), centreY + radius * std::sin(angleAtDisk), diskZ);
    }
  }

  // 2. the inner field, at the origin
  const double innerField = solenoid.innerFieldTesla;
  if (std::abs(innerField) < 1e-6)
    return std::nullopt;

  // 3. inner circle: same pT and tangent at the crossing, bent the other way, so its centre lies
  //    on the outward radial line of the outer circle, R_inner = pT / (0.3 B_inner) beyond it
  const double radialX = (crossing.x() - centreX) / radius;
  const double radialY = (crossing.y() - centreY) / radius;
  const double innerRadius = radiusFromPt(innerField, outer.ptGeV); // cm
  const double innerCentreX = crossing.x() + innerRadius * radialX;
  const double innerCentreY = crossing.y() + innerRadius * radialY;

  // 4. perigee: d0 = |inner centre| - R_inner, direction perpendicular to the line to the centre
  const double d0Cm = std::sqrt(innerCentreX * innerCentreX + innerCentreY * innerCentreY) - innerRadius;
  const double phi0 = wrapAngle(std::atan2(innerCentreY, innerCentreX) - innerSign * M_PI / 2.0);
  const double omega = innerSign / (innerRadius * kMmPerCm); // 1/mm, signed

  // 5. z at the perigee: back from the crossing, the shorter way (< pi)
  const double angleOfCrossing = std::atan2(crossing.y() - innerCentreY, crossing.x() - innerCentreX);
  const double angleOfPerigee = std::atan2(-innerCentreY, -innerCentreX);
  double backwardAngle = (omega > 0) ? wrapAngleZeroToTwoPi(angleOfCrossing - angleOfPerigee)
                                     : wrapAngleZeroToTwoPi(angleOfPerigee - angleOfCrossing);
  if (backwardAngle > M_PI)
    backwardAngle = 2.0 * M_PI - backwardAngle;
  const double z0Cm = crossing.z() - outer.tanLambda * (innerRadius * backwardAngle);

  return stateWithDefaultCovariance(d0Cm * kMmPerCm, phi0, omega, z0Cm * kMmPerCm, outer.tanLambda,
                                    edm4hep::TrackState::AtIP);
}

// ============================================================================================
// helpers (declared at the top)
// ============================================================================================

namespace {

  // weighted residuals sqrt(w) (|p - c| - R) of the hits from the circle (cx, cy, R) and their
  // Jacobian d residual / d(cx, cy, R) = -sqrt(w) (dx / distance, dy / distance, 1)
  void circleResidualsAndJacobian(const std::vector<Eigen::Vector2d>& points, const std::vector<double>& weights,
                                  const Eigen::Vector3d& circle, Eigen::MatrixXd& jacobian,
                                  Eigen::VectorXd& residuals) {
    const double centreX = circle(0);
    const double centreY = circle(1);
    const double radius = circle(2);
    for (size_t i = 0; i < points.size(); ++i) {
      const double dx = points[i].x() - centreX;
      const double dy = points[i].y() - centreY;
      const double distance = std::sqrt(dx * dx + dy * dy);
      const double sqrtWeight = std::sqrt(weights[i]);
      residuals(i) = sqrtWeight * (distance - radius);
      if (distance > 1e-12) {
        jacobian(i, 0) = -sqrtWeight * dx / distance;
        jacobian(i, 1) = -sqrtWeight * dy / distance;
        jacobian(i, 2) = -sqrtWeight;
      } else {
        jacobian(i, 0) = jacobian(i, 1) = jacobian(i, 2) = 0.0; // hit at the centre: no direction
      }
    }
  }

  // step 1 of makeStateFromCircle: the hits sorted by distance from the origin
  std::vector<StateHitPoint> hitsByDistanceFromOrigin(const std::vector<EventHit>& hits,
                                                      const std::vector<size_t>& trackHits) {
    std::vector<StateHitPoint> points;
    for (size_t hitIndex : trackHits) {
      const auto& position = hits[hitIndex].positionMm;
      StateHitPoint point;
      point.x = position[0];
      point.y = position[1];
      point.z = position[2];
      point.distance = std::sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
      point.varianceZ = hits[hitIndex].varianceZMm2;
      points.push_back(point);
    }
    std::sort(points.begin(), points.end(),
              [](const StateHitPoint& left, const StateHitPoint& right) { return left.distance < right.distance; });
    return points;
  }

  // step 2: sense of rotation seen from +z (+1 = counter-clockwise) from the sum of the cross
  // products of consecutive hits about the centre; 0 when the sum is 0 (no sense)
  double rotationSense(const std::vector<StateHitPoint>& points, double centreX, double centreY) {
    double crossSum = 0.0;
    for (size_t i = 0; i + 1 < points.size(); ++i) {
      const double fromCentreX = points[i].x - centreX;
      const double fromCentreY = points[i].y - centreY;
      const double nextFromCentreX = points[i + 1].x - centreX;
      const double nextFromCentreY = points[i + 1].y - centreY;
      crossSum += fromCentreX * nextFromCentreY - fromCentreY * nextFromCentreX;
    }
    if (crossSum == 0.0)
      return 0.0;
    return (crossSum > 0.0) ? 1.0 : -1.0;
  }

  // step 5: angle turned along the travel, hit by hit (no 2 pi jumps), and from it the arc length
  // of every hit from the pca
  std::vector<double> arcLengthsFromPca(const std::vector<StateHitPoint>& points, double centreX, double centreY,
                                        double radius, double sense, double anglePca, size_t nearestHit) {
    std::vector<double> turned(points.size());
    double previousAngle = std::atan2(points[0].y - centreY, points[0].x - centreX);
    turned[0] = 0.0;
    for (size_t i = 1; i < points.size(); ++i) {
      const double angle = std::atan2(points[i].y - centreY, points[i].x - centreX);
      const double step = wrapAngle(sense * (angle - previousAngle));
      turned[i] = turned[i - 1] + step;
      previousAngle = angle;
    }
    const double angleNearest = std::atan2(points[nearestHit].y - centreY, points[nearestHit].x - centreX);
    const double pcaStep = wrapAngle(sense * (anglePca - angleNearest));
    const double turnedPca = turned[nearestHit] + pcaStep;
    std::vector<double> arcLength(points.size());
    for (size_t i = 0; i < points.size(); ++i)
      arcLength[i] = radius * (turned[i] - turnedPca);
    return arcLength;
  }

  // step 6: weighted line fit z = zPca + tanLambda * s (weights 1 / sigma_z^2; a disk hit measures
  // z exactly, its sigma_z^2 is 0 or rounding noise: floor (1 um)^2), about the weighted means (no
  // cancellation when one weight is much larger); false if the hits do not spread along the arc
  bool fitZAlongArc(const std::vector<StateHitPoint>& points, const std::vector<double>& arcLength, ZLineFit& line) {
    std::vector<double> weight(points.size());
    double sumW = 0.0, sumWS = 0.0, sumWZ = 0.0;
    for (size_t i = 0; i < points.size(); ++i) {
      weight[i] = 1.0 / std::max(points[i].varianceZ, 1e-6); // floor (1 um)^2, in mm^2
      sumW += weight[i];
      sumWS += weight[i] * arcLength[i];
      sumWZ += weight[i] * points[i].z;
    }
    const double meanS = sumWS / sumW;
    const double meanZ = sumWZ / sumW;
    double sumWDSDS = 0.0, sumWDSDZ = 0.0;
    for (size_t i = 0; i < points.size(); ++i) {
      const double deviationS = arcLength[i] - meanS;
      sumWDSDS += weight[i] * deviationS * deviationS;
      sumWDSDZ += weight[i] * deviationS * (points[i].z - meanZ);
    }
    if (!(sumWDSDS / sumW > 1e-6))
      return false; // the hits do not spread along the arc (weighted spread below 1 um)
    line.tanLambda = sumWDSDZ / sumWDSDS;
    line.zPca = meanZ - line.tanLambda * meanS;
    line.sumW = sumW;
    line.meanS = meanS;
    line.sumWDSDS = sumWDSDS;
    return true;
  }

  // step 7: covariance: d0, phi, omega from the circle covariance (gradients in x0, y0, R, cm),
  // with lower limits; z0, tanLambda from the line fit
  void setStateCovariance(const CircleFit& circle, double centreToReferenceX, double centreToReferenceY,
                          double centreToReference, double omega, const ZLineFit& line, edm4hep::TrackState& state) {
    const double centreToReferenceSquaredCm2 = centreToReference * centreToReference / 100.0;
    const Eigen::Vector3d gradientD0(-centreToReferenceX / centreToReference, -centreToReferenceY / centreToReference,
                                     -1.0); // cm/cm
    const Eigen::Vector3d gradientPhi(centreToReferenceY / 10.0 / centreToReferenceSquaredCm2,
                                      -centreToReferenceX / 10.0 / centreToReferenceSquaredCm2, 0.0); // rad/cm
    const Eigen::Vector3d gradientOmega(0.0, 0.0, 1.0 / (circle.radius * circle.radius * 10.0));      // (1/mm)/cm
    const Eigen::Matrix3d& covariance = circle.covariance;
    // lower limits: 0.01 mm^2 (d0), 0.01 rad^2 (phi), (10% of omega)^2 (omega)
    const double varianceD0 = std::max(gradientD0.dot(covariance * gradientD0) * 100.0, 0.01);
    const double variancePhi = std::max(gradientPhi.dot(covariance * gradientPhi), 0.01);
    const double varianceOmega = std::max(gradientOmega.dot(covariance * gradientOmega), std::pow(0.1 * omega, 2));

    using Parameter = edm4hep::TrackParams;
    for (int row = 0; row < 6; ++row)
      for (int column = 0; column <= row; ++column)
        state.setCovMatrix(0.0f, static_cast<Parameter>(row), static_cast<Parameter>(column));
    state.setCovMatrix(static_cast<float>(varianceD0), Parameter::d0, Parameter::d0);
    state.setCovMatrix(static_cast<float>(variancePhi), Parameter::phi, Parameter::phi);
    state.setCovMatrix(static_cast<float>(varianceOmega), Parameter::omega, Parameter::omega);
    // line fit, with S = sum w (s - meanS)^2: var(tanLambda) = 1/S, var(z0) = 1/sumW + meanS^2/S,
    // cov(z0, tanLambda) = -meanS/S
    state.setCovMatrix(static_cast<float>(1.0 / line.sumW + line.meanS * line.meanS / line.sumWDSDS), Parameter::z0,
                       Parameter::z0);
    state.setCovMatrix(static_cast<float>(1.0 / line.sumWDSDS), Parameter::tanLambda, Parameter::tanLambda);
    state.setCovMatrix(static_cast<float>(-line.meanS / line.sumWDSDS), Parameter::z0, Parameter::tanLambda);
  }

  // crossing (cm) of the outer circle with the solenoid barrel, first reached going back from the
  // hit and on its IP side in z (else the other one, with tryOtherCrossing); false if none
  bool solenoidBarrelCrossing(const OuterHelix& outer, double solenoidRadiusMm, bool tryOtherCrossing,
                              Eigen::Vector3d& crossingCm) {
    const double solenoidRadiusCm = solenoidRadiusMm / kMmPerCm;
    const double centreX = outer.centreXCm;
    const double centreY = outer.centreYCm;
    const double radius = outer.radiusCm;
    const double centreDistance = std::sqrt(centreX * centreX + centreY * centreY); // from the IP
    if (centreDistance < 1e-6)
      return false;
    if (radius <= std::abs(centreDistance - solenoidRadiusCm) || radius >= centreDistance + solenoidRadiusCm)
      return false; // the circles miss each other
    double crossingX[2];
    double crossingY[2];
    circleCylinderCrossings(centreX, centreY, radius, centreDistance, solenoidRadiusCm, crossingX, crossingY);

    // the two crossings in the order reached going back from the hit
    const double hitX = outer.innermostHitCm.x();
    const double hitY = outer.innermostHitCm.y();
    const double hitZ = outer.innermostHitCm.z();
    const double angleOfHit = std::atan2(hitY - centreY, hitX - centreX);
    double backwardAngle[2];
    for (int i = 0; i < 2; ++i)
      backwardAngle[i] = backwardArcAngle(centreX, centreY, angleOfHit, outer.clockwise, crossingX[i], crossingY[i]);
    const int first = (backwardAngle[0] <= backwardAngle[1]) ? 0 : 1;
    const int order[2] = {first, 1 - first};

    // the first (with tryOtherCrossing also the other) crossing on the IP side of the hit in z
    const int crossingsToTry = tryOtherCrossing ? 2 : 1;
    for (int attempt = 0; attempt < crossingsToTry; ++attempt) {
      const int crossingIndex = order[attempt];
      const double z = hitZ - outer.tanLambda * (radius * backwardAngle[crossingIndex]);
      if (!(std::abs(z) > std::abs(hitZ) + 1e-3)) { // 10 um (cm) tolerance
        crossingCm = Eigen::Vector3d(crossingX[crossingIndex], crossingY[crossingIndex], z);
        return true;
      }
    }
    return false;
  }

  // The two crossings of a circle (centre at the distance d from the origin, radius R) with the
  // solenoid cylinder (radius Rs): on the line from the origin to the centre at
  // h = (Rs^2 + d^2 - R^2) / (2 d), and +-l = sqrt(Rs^2 - h^2) sideways.
  void circleCylinderCrossings(double centreX, double centreY, double radius, double centreDistance,
                               double solenoidRadius, double crossingX[2], double crossingY[2]) {
    const double along =
        (solenoidRadius * solenoidRadius + centreDistance * centreDistance - radius * radius) / (2.0 * centreDistance);
    const double sideways = std::sqrt(std::max(0.0, solenoidRadius * solenoidRadius - along * along));
    const double unitX = centreX / centreDistance, unitY = centreY / centreDistance;
    const double perpendicularX = -unitY, perpendicularY = unitX;
    crossingX[0] = along * unitX + sideways * perpendicularX;
    crossingY[0] = along * unitY + sideways * perpendicularY;
    crossingX[1] = along * unitX - sideways * perpendicularX;
    crossingY[1] = along * unitY - sideways * perpendicularY;
  }

  // angle on the circle from the hit (at angleOfHit about the centre) back to the point, in [0, 2pi),
  // against the direction of travel (the track came from the IP side)
  double backwardArcAngle(double centreX, double centreY, double angleOfHit, bool clockwise, double pointX,
                          double pointY) {
    const double angleOfPoint = std::atan2(pointY - centreY, pointX - centreX);
    return clockwise ? wrapAngleZeroToTwoPi(angleOfPoint - angleOfHit)
                     : wrapAngleZeroToTwoPi(angleOfHit - angleOfPoint);
  }

} // namespace

} // namespace CTF
