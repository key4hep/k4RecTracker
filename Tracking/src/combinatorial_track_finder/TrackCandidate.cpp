// The circle fit (see TrackCandidate.h).
#include "TrackCandidate.h"
#include "TrackingMath.h"

#include <cmath>

namespace k4rt {

namespace {
  constexpr int kMaxGaussNewtonSteps = 10;
  constexpr double kGaussNewtonTolerance = 1e-6; // cm: stop when a step moves (cx, cy, R) less

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
} // namespace

std::optional<CircleFit> fitCircle(const edm4hep::TrackerHitPlaneCollection& hits,
                                   const std::vector<size_t>& hitIndices, double defaultSigmaCm) {
  const int numberOfHits = static_cast<int>(hitIndices.size());
  if (numberOfHits < 3)
    return std::nullopt;

  // positions (cm) and weights 1 / du^2 (one radial residual per hit)
  std::vector<Eigen::Vector2d> points(numberOfHits);
  std::vector<double> sigmas(numberOfHits);
  std::vector<double> weights(numberOfHits);
  for (int i = 0; i < numberOfHits; ++i) {
    const auto hit = hits[hitIndices[i]];
    const auto& position = hit.getPosition();
    points[i] = {position[0] / kMmPerCm, position[1] / kMmPerCm};
    double du = hit.getDu(); // mm
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
  for (int step = 0; step < kMaxGaussNewtonSteps; ++step) {
    circleResidualsAndJacobian(points, weights, circle, jacobian, residuals);
    const Eigen::Vector3d delta = -(jacobian.transpose() * jacobian).ldlt().solve(jacobian.transpose() * residuals);
    circle += delta;
    if (delta.norm() < kGaussNewtonTolerance)
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

} // namespace k4rt
