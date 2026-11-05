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
  // to (0, 2pi]. Unlike wrapAngleZeroToTwoPi (TrackingMath.h), an angle of exactly 0
  // becomes a full turn: a crossing at the starting point itself is the next turn.
  double wrapPositive(double a) {
    a = std::fmod(a, kTwoPi);
    if (a <= 0.0)
      a += kTwoPi;
    return a;
  }
} // namespace

Eigen::Vector3d Helix::at(double s) const {
  const double a = angleAt(s);
  return {cx + R * std::cos(a), cy + R * std::sin(a), zRef + tanL * s};
}

Eigen::Vector3d Helix::dirAt(double s) const {
  const double a = angleAt(s);
  Eigen::Vector3d d(-rot * std::sin(a), rot * std::cos(a), tanL);
  return d.normalized();
}

bool helixFromHits(double cx, double cy, double R, const std::vector<Eigen::Vector3d>& hits, Helix& out,
                   std::vector<double>* sOfHits) {
  const int n = static_cast<int>(hits.size());
  if (n < 2 || !(R > 0.0))
    return false;
  out.cx = cx;
  out.cy = cy;
  out.R = R;

  // Unwrapped angle of every hit about the centre, accumulated hit-to-hit.
  std::vector<double> ang(n);
  ang[0] = std::atan2(hits[0].y() - cy, hits[0].x() - cx);
  for (int i = 1; i < n; ++i) {
    const double a = std::atan2(hits[i].y() - cy, hits[i].x() - cx);
    ang[i] = ang[i - 1] + wrapAngle(a - ang[i - 1]);
  }
  out.rot = (ang[n - 1] - ang[0] >= 0.0) ? +1 : -1;
  out.phiRef = ang[0];

  // s_i = rot * R * (ang_i - ang_0) >= 0 along the track; linear fit z = zRef + tanL * s.
  std::vector<double> s(n);
  double sS = 0, sZ = 0, sSS = 0, sSZ = 0;
  for (int i = 0; i < n; ++i) {
    s[i] = out.rot * R * (ang[i] - ang[0]);
    sS += s[i];
    sZ += hits[i].z();
    sSS += s[i] * s[i];
    sSZ += s[i] * hits[i].z();
  }
  const double det = n * sSS - sS * sS;
  if (std::abs(det) > 1e-12) {
    out.zRef = (sZ * sSS - sS * sSZ) / det;
    out.tanL = (n * sSZ - sS * sZ) / det;
  } else {
    out.zRef = sZ / n;
    out.tanL = 0.0;
  }
  if (sOfHits)
    *sOfHits = s;
  return true;
}

double arcLengthNear(const Helix& h, const Eigen::Vector3d& p, double sGuess) {
  const double a = std::atan2(p.y() - h.cy, p.x() - h.cx);
  const double a0 = h.angleAt(sGuess);
  // smallest angular move from the guess, converted back to arc length
  return sGuess + h.rot * h.R * wrapAngle(a - a0);
}

Crossing crossCylinder(const Helix& h, double Rc, double sFrom, int dir, double sEps, double sMax) {
  Crossing best;
  const double d = std::hypot(h.cx, h.cy);
  if (d < 1e-12)
    return best; // circle centred on the axis
  if (d > h.R + Rc || d < std::abs(h.R - Rc))
    return best;                                              // never reaches Rc
  const double a = (Rc * Rc - h.R * h.R + d * d) / (2.0 * d); // along the centre direction
  const double h2 = Rc * Rc - a * a;
  const double hh = std::sqrt(std::max(0.0, h2));
  const double ux = h.cx / d, uy = h.cy / d;

  const double a0 = h.angleAt(sFrom);
  double bestDs = std::numeric_limits<double>::infinity();
  for (int sgn = -1; sgn <= 1; sgn += 2) {
    const double px = a * ux - sgn * hh * uy;
    const double py = a * uy + sgn * hh * ux;
    const double ap = std::atan2(py - h.cy, px - h.cx);
    // travelling dir along s turns the angle by rot*dir per unit s
    const double dAng = wrapPositive(h.rot * dir * (ap - a0));
    const double ds = h.R * dAng;
    if (ds < sEps || ds > sMax)
      continue;
    if (ds < bestDs)
      bestDs = ds;
  }
  if (!std::isfinite(bestDs))
    return best;
  best.ok = true;
  best.s = sFrom + dir * bestDs;
  best.pos = h.at(best.s);
  return best;
}

Crossing crossPlaneZ(const Helix& h, double zPlane, double sFrom, int dir, double sEps, double sMax) {
  Crossing c;
  if (std::abs(h.tanL) < 1e-12)
    return c;
  const double s = (zPlane - h.zRef) / h.tanL;
  const double ds = dir * (s - sFrom);
  if (ds < sEps || ds > sMax)
    return c;
  c.ok = true;
  c.s = s;
  c.pos = h.at(s);
  return c;
}

Crossing crossLayer(const Helix& h, const LayerDescriptor& L, double sFrom, int dir, double tol, double sMax) {
  Crossing c;
  if (L.nSurfaces == 0)
    return c;
  if (L.isBarrel()) {
    c = crossCylinder(h, L.rMean_cm, sFrom, dir, 1e-6, sMax);
    if (c.ok && (c.pos.z() < L.zMin_cm - tol || c.pos.z() > L.zMax_cm + tol))
      c.ok = false;
  } else {
    c = crossPlaneZ(h, L.zMean_cm, sFrom, dir, 1e-6, sMax);
    if (c.ok) {
      const double r = std::hypot(c.pos.x(), c.pos.y());
      if (r < L.rMin_cm - tol || r > L.rMax_cm + tol)
        c.ok = false;
    }
  }
  return c;
}

Crossing refineToPlane(const Helix& h, const Eigen::Vector3d& o, const Eigen::Vector3d& n, double sStart, int maxIter,
                       double tol) {
  Crossing c;
  double s = sStart;
  for (int it = 0; it < maxIter; ++it) {
    const Eigen::Vector3d p = h.at(s);
    const double f = n.dot(p - o);
    if (std::abs(f) < tol) {
      c.ok = true;
      c.s = s;
      c.pos = p;
      return c;
    }
    // d/ds of the position, NOT normalised: (dx/ds, dy/ds) has unit length, dz/ds = tanL
    const double a = h.angleAt(s);
    const Eigen::Vector3d dp(-h.rot * std::sin(a), h.rot * std::cos(a), h.tanL);
    const double fp = n.dot(dp);
    if (std::abs(fp) < 1e-12)
      return c; // travelling parallel to the plane
    s -= f / fp;
  }
  return c;
}

} // namespace k4rt
