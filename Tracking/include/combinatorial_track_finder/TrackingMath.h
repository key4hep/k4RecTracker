// TrackingMath - the formulas shared by the CombinatorialTrackFinder code (one version for all):
// angle wrapping, pT <-> radius <-> curvature, the Highland angle and the kink chi2, the mm <-> cm
// factors, and the settings shared by both search strategies. Needs only Eigen (unit-testable).
//
// The kink chi2: multiple scattering enters as KINKS. For three consecutive hits a-b-c the change
// of direction at b (in the bending plane and in dip) is compared with the turn of the fitted helix,
// each with the expected spread
//     sigma^2 = theta0(b)^2 (Highland, x/X0 of b's layer / cos(incidence)) + measurement errors.
// Kinks at different layers are independent, so
//     chi2 = sum over triplets of (dphi/sigma_phi)^2 + (dlambda/sigma_lambda)^2,  ndf = 2 (n - 2) - 1
// gives chi2/ndf ~ 1 for a real track at any momentum (the idea of General Broken Lines, used here
// only as a gate and a ranking). Inflating every hit error by the upstream scattering instead
// over-counts by ~3x, because those displacements are correlated.
//
// Units: radius in cm, curvature omega in 1/mm (EDM4hep), field in tesla, momentum in GeV.
#ifndef TRACKING_TRACKINGMATH_H
#define TRACKING_TRACKINGMATH_H

#include <Eigen/Dense>

#include <cmath>
#include <vector>

namespace k4rt {

// EDM4hep positions are in mm, DD4hep and the circle fit in cm. x / kMmPerCm and x * kCmPerMm differ
// in the last bit: each place keeps the form it always used (unchanged outputs).
constexpr double kMmPerCm = 10.0;
constexpr double kCmPerMm = 0.1;

/// Angle wrapped into (-pi, pi].
inline double wrapAngle(double angle) {
  while (angle > M_PI)
    angle -= 2.0 * M_PI;
  while (angle <= -M_PI)
    angle += 2.0 * M_PI;
  return angle;
}

/// Angle wrapped into [0, 2 pi).
inline double wrapAngleZeroToTwoPi(double angle) {
  angle = std::fmod(angle, 2.0 * M_PI);
  if (angle < 0.0)
    angle += 2.0 * M_PI;
  return angle;
}

// pT [GeV] = 0.299792458 * |B| [T] * R [m] for charge 1 (0.299792458 = c / 1e9)
constexpr double kPtPerTeslaMetre = 0.299792458;

/// pT (GeV) of a circle of radius radiusCm (cm) in the field bFieldT (T).
inline double ptFromRadius(double bFieldT, double radiusCm) {
  return kPtPerTeslaMetre * std::abs(bFieldT) * radiusCm / 100.0;
}

/// Bending radius (cm) for pT ptGeV (GeV) in the field bFieldT (T).
inline double radiusFromPt(double bFieldT, double ptGeV) {
  return ptGeV * 100.0 / (kPtPerTeslaMetre * std::abs(bFieldT));
}

/// pT (GeV) from the EDM4hep curvature omega (1/mm, either sign).
inline double ptFromOmega(double bFieldT, double omegaPerMm) {
  return kPtPerTeslaMetre * std::abs(bFieldT) / (std::abs(omegaPerMm) * 1000.0);
}

/// |omega| (1/mm) for pT ptGeV (GeV); the caller gives it the sign of the charge.
inline double omegaFromPt(double bFieldT, double ptGeV) {
  return kPtPerTeslaMetre * std::abs(bFieldT) / (ptGeV * 1000.0);
}

/// Settings used by both search strategies, filled once from the options (so they cannot disagree).
struct CommonSearchConfig {
  int minTrackHits = 3;                ///< MinTrackHits
  double minPtGeV = 0.0;               ///< MinComboPT (<= 0 = off)
  double maxPairHitDistanceCm = 300.0; ///< MaxPairHitDistance
  double bFieldTesla = 0.0;            ///< Bz at the IP, from the field map
  double sigmaHitDefaultCm = 0.04;     ///< SigmaHitDefault: hit resolution when a hit has no du
};

/// One hit as the kink chi2 sees it.
struct KinkPoint {
  Eigen::Vector3d position = Eigen::Vector3d::Zero(); ///< cm
  double sigmaTransverseCm = 1e-3;                    ///< measurement error, transverse
  double sigmaZCm = 1e-3;                             ///< measurement error, along z / dip
  double xOverX0 = 0.0;                               ///< material of the hit's layer at normal incidence
  Eigen::Vector3d normal = Eigen::Vector3d::UnitX();  ///< sensor normal, for the incidence angle
};

/// The kink chi2 of a track, its ndf, and the chi2 of every triplet.
struct KinkResult {
  double chi2 = 0.0;
  int ndf = 0;
  std::vector<double> perTriplet; ///< chi2 of the triplet centred at hit i (0 at both ends)
  double chi2ndf() const { return ndf > 0 ? chi2 / ndf : 0.0; }
};

/// The kink chi2 described above, for a given field (the momentum, and so the expected
/// scattering, follows from the radius and the field).
class KinkChi2Calculator {
public:
  explicit KinkChi2Calculator(double bFieldTesla) : m_bFieldTesla(bFieldTesla) {}

  /// Highland multiple-scattering angle (one plane, with the log term) for beta = 1 (the particle
  /// type is unknown): within 1% of a muon or pion above 1 GeV, 6-10% smaller at 0.3 GeV.
  static double highlandTheta0(double radiationLengths, double momentumGeV);

  /// Kink chi2 (2 dof) of the triplet first-middle-last (order of travel) on a helix of radius
  /// radiusCm and rotation sense (+1 = anticlockwise seen from +z). sigmaRadiusCm: error of the
  /// radius, when it comes from other hits than these three (0 when fitted to them).
  double triplet(const KinkPoint& first, const KinkPoint& middle, const KinkPoint& last, double radiusCm, int rotation,
                 double sigmaRadiusCm = 0.0) const;

  /// Kink chi2 of hits ordered along the track: the sum over all triplets.
  KinkResult track(const std::vector<KinkPoint>& points, double radiusCm, int rotation) const;

private:
  double m_bFieldTesla = 0.0; ///< tesla, from the field map
};

} // namespace k4rt

#endif // TRACKING_TRACKINGMATH_H
