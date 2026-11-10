// The kink chi2 and the Highland angle (see TrackingMath.h).
#include "TrackingMath.h"

#include <algorithm>
#include <cmath>

namespace k4rt {

double KinkChi2Calculator::highlandTheta0(double x, double p) {
  if (!(x > 0.0) || !(p > 0.0))
    return 0.0;
  // Highland formula: theta0 = 13.6 MeV / (beta c p) * sqrt(x/X0) * (1 + 0.038 ln(x/X0)),
  // x = material thickness in radiation lengths, p in GeV, beta = 1. The log factor is kept
  // above 0.1 so that very thin layers still give a small positive angle.
  const double logTerm = 1.0 + 0.038 * std::log(x);
  return 0.0136 / p * std::sqrt(x) * std::max(logTerm, 0.1);
}

double KinkChi2Calculator::triplet(const KinkPoint& a, const KinkPoint& b, const KinkPoint& c, double R, int rot,
                                   double sigmaR) const {
  // the two chords a->b and b->c in the transverse plane, and their lengths
  const Eigen::Vector2d u1(b.pos.x() - a.pos.x(), b.pos.y() - a.pos.y());
  const Eigen::Vector2d u2(c.pos.x() - b.pos.x(), c.pos.y() - b.pos.y());
  const double l1 = u1.norm(), l2 = u2.norm();
  if (l1 < 1e-6 || l2 < 1e-6 || !(R > 0.0))
    return 0.0;

  // 1. kink in the bending plane: the measured turn between the two chords (signed angle from
  //    the cross and dot products) minus the turn the helix itself makes. A chord of length l on
  //    a circle of radius R makes the angle asin(l / 2R) with the tangent at each end, so the
  //    helix turns by h1 + h2 between the chords.
  const double turn = std::atan2(u1.x() * u2.y() - u1.y() * u2.x(), u1.dot(u2));
  const double h1 = std::asin(std::min(1.0, l1 / (2.0 * R))), h2 = std::asin(std::min(1.0, l2 / (2.0 * R)));
  const double dphi = wrapAngle(turn - rot * (h1 + h2));

  // 2. kink in the dip: change of the slope dz/ds between the two segments, with s the arc
  //    length on the circle (2 R h for a chord), not the chord length
  const double s1 = 2.0 * R * h1, s2 = 2.0 * R * h2;
  const double dz1 = b.pos.z() - a.pos.z(), dz2 = c.pos.z() - b.pos.z();
  const double dlam = std::atan2(dz2, s2) - std::atan2(dz1, s1);

  // 3. expected multiple scattering at b: momentum from the radius (pT) and the mean dip
  //    (p = pT / cos(lambda)); the material seen grows as 1 / cos(incidence angle) on the sensor
  //    (limited to 1 / 0.05 for grazing tracks)
  const double tanl = (dz1 + dz2) / (s1 + s2);
  const double cosl = 1.0 / std::sqrt(1.0 + tanl * tanl);
  const double pT = ptFromRadius(m_bFieldTesla, R); // GeV, R in cm
  const double p = pT / cosl;
  const Eigen::Vector3d dir = (c.pos - a.pos).normalized();
  const double cosInc = std::max(0.05, std::abs(b.normal.normalized().dot(dir)));
  const double th0 = highlandTheta0(b.xOverX0 / cosInc, p);

  // 4. expected spread of the bending-plane kink, the sum of:
  //    - the uncertainty of the helix turn from sigmaR: d/dR asin(l/2R) = -l / (2 R^2 sqrt(1 - (l/2R)^2)),
  //    - the scattering angle, projected on the transverse plane (theta0 / cos(lambda)),
  //    - the measurement errors: a transverse shift sigma of a hit turns a chord of length l by
  //      sigma / l; the middle hit b moves both chords (1/l1 + 1/l2)
  const double e1 = 1.0 / l1, e2 = 1.0 / l2;
  const auto dh = [R](double l) {
    const double x = std::min(0.999999, l / (2.0 * R));
    return -l / (2.0 * R * R * std::sqrt(1.0 - x * x));
  };
  const double sTurnR = (dh(l1) + dh(l2)) * sigmaR;
  const double sPhi2 = sTurnR * sTurnR + (th0 / cosl) * (th0 / cosl) + a.sigT_cm * a.sigT_cm * e1 * e1 +
                       b.sigT_cm * b.sigT_cm * (e1 + e2) * (e1 + e2) + c.sigT_cm * c.sigT_cm * e2 * e2;
  // 5. expected spread of the dip kink: the scattering angle, plus the z measurement errors
  //    turned into slope errors (sigma_z / s, the middle hit in both segments); d(lambda) =
  //    cos^2(lambda) d(tan lambda), hence the cos^4 on the variance
  const double f1 = 1.0 / s1, f2 = 1.0 / s2, c4 = cosl * cosl * cosl * cosl;
  const double sLam2 =
      th0 * th0 + c4 * (a.sigZ_cm * a.sigZ_cm * f1 * f1 + b.sigZ_cm * b.sigZ_cm * (f1 + f2) * (f1 + f2) +
                        c.sigZ_cm * c.sigZ_cm * f2 * f2);
  // 6. chi2 of the two kinks, each divided by its expected variance (2 dof)
  return dphi * dphi / sPhi2 + dlam * dlam / sLam2;
}

KinkResult KinkChi2Calculator::track(const std::vector<KinkPoint>& pts, double R, int rot) const {
  KinkResult r;
  const int n = static_cast<int>(pts.size());
  r.perTriplet.assign(std::max(n, 0), 0.0);
  if (n < 3)
    return r;
  // one kink at every inner hit; the kinks are independent, so their chi2 add up
  for (int i = 1; i + 1 < n; ++i) {
    const double c = triplet(pts[i - 1], pts[i], pts[i + 1], R, rot);
    r.perTriplet[i] = c;
    r.chi2 += c;
  }
  // 2 measured angles per kink, minus 1 for the radius taken from the fit
  r.ndf = 2 * (n - 2) - 1;
  return r;
}

} // namespace k4rt
