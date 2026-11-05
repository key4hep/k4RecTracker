// The kink chi2 and the Highland angle (see TrackingMath.h).
#include "TrackingMath.h"

#include <algorithm>
#include <cmath>

namespace k4rt {

double highlandTheta0(double x, double p) {
  if (!(x > 0.0) || !(p > 0.0))
    return 0.0;
  // x = material thickness in radiation lengths, p in GeV, beta = 1
  const double logTerm = 1.0 + 0.038 * std::log(x);
  return 0.0136 / p * std::sqrt(x) * std::max(logTerm, 0.1);
}

double tripletKinkChi2(const KinkPoint& a, const KinkPoint& b, const KinkPoint& c, double R, int rot,
                       const KinkConfig& cfg, double sigmaR) {
  const Eigen::Vector2d u1(b.pos.x() - a.pos.x(), b.pos.y() - a.pos.y());
  const Eigen::Vector2d u2(c.pos.x() - b.pos.x(), c.pos.y() - b.pos.y());
  const double l1 = u1.norm(), l2 = u2.norm();
  if (l1 < 1e-6 || l2 < 1e-6 || !(R > 0.0))
    return 0.0;

  // bending plane: measured turn between the chords vs the helix's own turn
  const double turn = std::atan2(u1.x() * u2.y() - u1.y() * u2.x(), u1.dot(u2));
  const double h1 = std::asin(std::min(1.0, l1 / (2.0 * R))), h2 = std::asin(std::min(1.0, l2 / (2.0 * R)));
  const double dphi = wrapAngle(turn - rot * (h1 + h2));

  // dip: change of slope between the two segments (arc lengths, not chords)
  const double s1 = 2.0 * R * h1, s2 = 2.0 * R * h2;
  const double dz1 = b.pos.z() - a.pos.z(), dz2 = c.pos.z() - b.pos.z();
  const double dlam = std::atan2(dz2, s2) - std::atan2(dz1, s1);

  // expected scattering at b
  const double tanl = (dz1 + dz2) / (s1 + s2);
  const double cosl = 1.0 / std::sqrt(1.0 + tanl * tanl);
  const double pT = ptFromRadius(cfg.bField_T, R); // GeV, R in cm
  const double p = pT / cosl;
  const Eigen::Vector3d dir = (c.pos - a.pos).normalized();
  const double cosInc = std::max(0.05, std::abs(b.normal.normalized().dot(dir)));
  const double th0 = highlandTheta0(b.xOverX0 / cosInc, p);

  const double e1 = 1.0 / l1, e2 = 1.0 / l2;
  // error of the expected turn from the curvature: d/dR asin(l/2R) = -l / (2 R^2 sqrt(1 - (l/2R)^2))
  const auto dh = [R](double l) {
    const double x = std::min(0.999999, l / (2.0 * R));
    return -l / (2.0 * R * R * std::sqrt(1.0 - x * x));
  };
  const double sTurnR = (dh(l1) + dh(l2)) * sigmaR;
  const double sPhi2 = sTurnR * sTurnR + (th0 / cosl) * (th0 / cosl) + a.sigT_cm * a.sigT_cm * e1 * e1 +
                       b.sigT_cm * b.sigT_cm * (e1 + e2) * (e1 + e2) + c.sigT_cm * c.sigT_cm * e2 * e2;
  const double f1 = 1.0 / s1, f2 = 1.0 / s2, c4 = cosl * cosl * cosl * cosl;
  const double sLam2 =
      th0 * th0 + c4 * (a.sigZ_cm * a.sigZ_cm * f1 * f1 + b.sigZ_cm * b.sigZ_cm * (f1 + f2) * (f1 + f2) +
                        c.sigZ_cm * c.sigZ_cm * f2 * f2);
  return dphi * dphi / sPhi2 + dlam * dlam / sLam2;
}

KinkResult kinkChi2(const std::vector<KinkPoint>& pts, double R, int rot, const KinkConfig& cfg) {
  KinkResult r;
  const int n = static_cast<int>(pts.size());
  r.perTriplet.assign(std::max(n, 0), 0.0);
  if (n < 3)
    return r;
  for (int i = 1; i + 1 < n; ++i) {
    const double c = tripletKinkChi2(pts[i - 1], pts[i], pts[i + 1], R, rot, cfg);
    r.perTriplet[i] = c;
    r.chi2 += c;
  }
  r.ndf = 2 * (n - 2) - 1;
  return r;
}

} // namespace k4rt
