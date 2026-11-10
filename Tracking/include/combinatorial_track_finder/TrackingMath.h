// TrackingMath - the formulas shared by the CombinatorialTrackFinder code, kept in one place
// so every part of the code uses exactly the same version:
//   - wrapping an angle into a fixed range
//   - converting between transverse momentum, circle radius and curvature in a field
//   - the multiple-scattering angle (Highland) and the kink chi2 of SeedAndExtend
//   - the settings shared by both search strategies (CommonSearchConfig)
//
// It needs only Eigen, so the standalone unit tests can use it too.
//
// The kink chi2: a multiple-scattering-aware chi2 for track finding, built from kinks.
//
// The circle fit weights hits by their measurement error only (micrometres). At pT ~ 0.5 GeV
// multiple scattering moves hits by 0.1-1 mm, so real low-pT tracks fail any chi2 cut tuned on
// it and extension stops part-way. Inflating every hit's error with the scattering of all
// upstream material over-counts by ~3x, because those displacements are correlated.
//
// Here the scattering enters where it is independent: as KINKS. For three consecutive hits
// a-b-c the change of direction at b, in the bending plane and in dip, is compared with what
// the fitted helix predicts. Each kink has an expected spread
//     sigma^2 = theta0(b)^2  (Highland, x/X0 of b's layer / cos(incidence), at the track's p)
//             + the measurement errors of a, b, c propagated to the angle.
// Kinks at different layers are statistically independent, so
//     chi2 = sum over triplets of (dphi/sigma_phi)^2 + (dlambda/sigma_lambda)^2,
//     ndf  = 2 (n - 2) - 1          (the radius is taken from the fit)
// has chi2/ndf ~ 1 for a real track at any momentum. This is the idea behind General Broken
// Lines, used here only as a gate and a ranking, not for the final track parameters.
//
// Units: radius in cm, curvature omega in 1/mm (EDM4hep), field in tesla, momentum in GeV.
#ifndef TRACKING_TRACKINGMATH_H
#define TRACKING_TRACKINGMATH_H

#include <Eigen/Dense>

#include <cmath>
#include <vector>

namespace k4rt {

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

// ---- transverse momentum, bending radius and curvature in a solenoid field ----
// A particle of charge 1 and transverse momentum pT bends on a circle of radius R in a
// field B with  pT [GeV] = 0.299792458 * |B| [T] * R [m]  (0.299792458 = c / 1e9).

/// pT per tesla and per metre of bending radius (GeV / (T m))
constexpr double kPtPerTeslaMetre = 0.299792458;

/// Transverse momentum (GeV) of a circle of radius radiusCm (cm) in the field bFieldT (T).
inline double ptFromRadius(double bFieldT, double radiusCm) {
  return kPtPerTeslaMetre * std::abs(bFieldT) * radiusCm / 100.0;
}

/// Bending radius (cm) of a particle of transverse momentum ptGeV (GeV) in the field bFieldT (T).
inline double radiusFromPt(double bFieldT, double ptGeV) {
  return ptGeV * 100.0 / (kPtPerTeslaMetre * std::abs(bFieldT));
}

/// Transverse momentum (GeV) from the EDM4hep curvature omega (1/mm, either sign).
inline double ptFromOmega(double bFieldT, double omegaPerMm) {
  return kPtPerTeslaMetre * std::abs(bFieldT) / (std::abs(omegaPerMm) * 1000.0);
}

/// Size of the EDM4hep curvature |omega| (1/mm) for the transverse momentum ptGeV (GeV);
/// the caller gives it the sign of the charge.
inline double omegaFromPt(double bFieldT, double ptGeV) {
  return kPtPerTeslaMetre * std::abs(bFieldT) / (ptGeV * 1000.0);
}

// ---- kink chi2 (multiple scattering included) ----

/// Settings used by both search strategies (Exhaustive and SeedAndExtend), set once from the
/// CombinatorialTrackFinder options so that the two searches cannot disagree.
struct CommonSearchConfig {
  int minTrackHits = 3;                ///< MinTrackHits: fewest hits of a track
  double minPtGeV = 0.0;               ///< MinComboPT: minimum pT (GeV); <= 0 = off
  double maxPairHitDistanceCm = 300.0; ///< MaxPairHitDistance: largest 3D distance between hits (cm)
  double bFieldTesla = 0.0;            ///< Bz at the IP, from the field map
};

/// One hit as the kink chi2 sees it.
struct KinkPoint {
  Eigen::Vector3d pos = Eigen::Vector3d::Zero();     ///< cm
  double sigT_cm = 1e-3;                             ///< measurement error, transverse
  double sigZ_cm = 1e-3;                             ///< measurement error, along z / dip
  double xOverX0 = 0.0;                              ///< material of the hit's layer at normal incidence
  Eigen::Vector3d normal = Eigen::Vector3d::UnitX(); ///< sensor normal, for the incidence angle
};

/// The kink chi2 of a track, its ndf, and the chi2 of every triplet.
struct KinkResult {
  double chi2 = 0.0;
  int ndf = 0;
  std::vector<double> perTriplet; ///< chi2 of the triplet centred at hit i (0 at both ends)
  double chi2ndf() const { return ndf > 0 ? chi2 / ndf : 0.0; }
};

/// Computes the kink chi2 described above, for a given field (the momentum, and so the
/// expected scattering, follows from the radius and the field).
class KinkChi2Calculator {
public:
  explicit KinkChi2Calculator(double bFieldTesla) : m_bFieldTesla(bFieldTesla) {}

  /// Highland multiple-scattering angle (in one plane), with the log term, for a particle with
  /// beta = 1: the track finder does not know the particle type. Above p ~ 1 GeV this is
  /// within 1% of the angle for a muon or a pion (6-10% smaller at 0.3 GeV).
  static double highlandTheta0(double xOverX0, double p_GeV);

  /// Kink chi2 (2 dof) of the triplet a-b-c for a helix of radius R_cm and rotation sense rot
  /// (+1 = anticlockwise seen from +z, the order of travel being a -> b -> c).
  /// sigmaR_cm: uncertainty of R, propagated to the expected turn. Pass it when R comes from
  /// hits other than a, b, c (e.g. testing a new hit against the current track's curvature);
  /// leave it 0 when R was fitted to these hits.
  double triplet(const KinkPoint& a, const KinkPoint& b, const KinkPoint& c, double R_cm, int rot,
                 double sigmaR_cm = 0.0) const;

  /// Kink chi2 of hits ordered along the track: the sum over all triplets.
  KinkResult track(const std::vector<KinkPoint>& pts, double R_cm, int rot) const;

private:
  double m_bFieldTesla = 0.0; ///< tesla, from the field map
};

} // namespace k4rt

#endif // TRACKING_TRACKINGMATH_H
