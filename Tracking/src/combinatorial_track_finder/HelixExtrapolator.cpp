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
  // An angle brought into (0, 2pi]. Unlike wrapAngleZeroToTwoPi (TrackingMath.h), an angle of
  // exactly 0 becomes a full turn: a crossing at the starting point itself is the next turn.
  double wrapPositive(double a) {
    a = std::fmod(a, kTwoPi);
    if (a <= 0.0)
      a += kTwoPi;
    return a;
  }
} // namespace

Eigen::Vector3d Helix::at(double s) const {
  // a point on the circle at the angle reached after s, and the straight line in z
  const double a = angleAt(s);
  return {cx + R * std::cos(a), cy + R * std::sin(a), zRef + tanL * s};
}

Eigen::Vector3d Helix::dirAt(double s) const {
  // derivative of at(s): the transverse tangent of the circle (unit length, turned by rot) and
  // dz/ds = tanL; normalised to a unit 3D vector
  const double a = angleAt(s);
  Eigen::Vector3d d(-rot * std::sin(a), rot * std::cos(a), tanL);
  return d.normalized();
}

bool Helix::fromHits(double centreX, double centreY, double radius, const std::vector<Eigen::Vector3d>& hits,
                     Helix& out, std::vector<double>* sOfHits) {
  const int n = static_cast<int>(hits.size());
  if (n < 2 || !(radius > 0.0))
    return false;
  out.cx = centreX;
  out.cy = centreY;
  out.R = radius;

  // 1. angle of every hit about the centre, "unwrapped": each step is the smallest turn from the
  //    previous hit, so the accumulated angle keeps growing (or shrinking) beyond +-pi
  std::vector<double> ang(n);
  ang[0] = std::atan2(hits[0].y() - centreY, hits[0].x() - centreX);
  for (int i = 1; i < n; ++i) {
    const double a = std::atan2(hits[i].y() - centreY, hits[i].x() - centreX);
    ang[i] = ang[i - 1] + wrapAngle(a - ang[i - 1]);
  }
  // 2. the rotation sense: the sign of the total turn from the first to the last hit
  out.rot = (ang[n - 1] - ang[0] >= 0.0) ? +1 : -1;
  out.phiRef = ang[0];

  // 3. arc length of each hit, s_i = rot * R * (ang_i - ang_0) >= 0 along the track, and the
  //    least-squares straight line z = zRef + tanL * s through the points (s_i, z_i):
  //      tanL = (n * sum(s z) - sum(s) sum(z)) / det,  zRef = (sum(z) sum(s^2) - sum(s) sum(s z)) / det,
  //      det = n * sum(s^2) - sum(s)^2
  std::vector<double> s(n);
  double sS = 0, sZ = 0, sSS = 0, sSZ = 0;
  for (int i = 0; i < n; ++i) {
    s[i] = out.rot * radius * (ang[i] - ang[0]);
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
    // all hits at the same s (no lever arm): a flat line at the mean z
    out.zRef = sZ / n;
    out.tanL = 0.0;
  }
  if (sOfHits)
    *sOfHits = s;
  return true;
}

double Helix::arcLengthNear(const Eigen::Vector3d& p, double sGuess) const {
  // angle of p about the centre, compared with the angle at the guess; the smallest angular move
  // (wrapped to (-pi, pi]) is converted back to an arc length (angle * R, with the rotation sense)
  const double a = std::atan2(p.y() - cy, p.x() - cx);
  const double a0 = angleAt(sGuess);
  return sGuess + rot * R * wrapAngle(a - a0);
}

Crossing Helix::crossCylinder(double rCylinder, double sFrom, int dir, double sEps, double sMax) const {
  Crossing best;
  // Intersection of two circles in the transverse plane: the track circle (centre c, radius R)
  // and the cylinder (centre at the origin, radius Rc). With d = |c|:
  //   - no crossing if d > R + Rc (too far) or d < |R - Rc| (one circle inside the other);
  //   - otherwise the two crossing points lie on the line through the origin along c at the
  //     distance a = (Rc^2 - R^2 + d^2) / (2 d) from the origin, and +-hh perpendicular to it,
  //     with hh = sqrt(Rc^2 - a^2) (Pythagoras in the triangle origin - crossing - foot point).
  const double d = std::hypot(cx, cy);
  if (d < 1e-12)
    return best; // circle centred on the axis: its radius never changes
  if (d > R + rCylinder || d < std::abs(R - rCylinder))
    return best; // never reaches the radius rCylinder
  const double a = (rCylinder * rCylinder - R * R + d * d) / (2.0 * d);
  const double h2 = rCylinder * rCylinder - a * a;
  const double hh = std::sqrt(std::max(0.0, h2));
  const double ux = cx / d, uy = cy / d; // unit vector from the origin towards the centre

  // of the two crossing points, the one reached first when moving from sFrom in the direction dir
  const double a0 = angleAt(sFrom);
  double bestDs = std::numeric_limits<double>::infinity();
  for (int sgn = -1; sgn <= 1; sgn += 2) {
    const double px = a * ux - sgn * hh * uy;
    const double py = a * uy + sgn * hh * ux;
    const double ap = std::atan2(py - cy, px - cx);
    // moving by ds along dir turns the angle about the centre by rot * dir * ds / R; the turn
    // still needed to reach the crossing, in (0, 2pi], gives the arc length ds = R * turn
    const double dAng = wrapPositive(rot * dir * (ap - a0));
    const double ds = R * dAng;
    if (ds < sEps || ds > sMax)
      continue;
    if (ds < bestDs)
      bestDs = ds;
  }
  if (!std::isfinite(bestDs))
    return best;
  best.ok = true;
  best.s = sFrom + dir * bestDs;
  best.pos = at(best.s);
  return best;
}

Crossing Helix::crossPlaneZ(double zPlane, double sFrom, int dir, double sEps, double sMax) const {
  Crossing c;
  // z = zRef + tanL * s, solved for s; impossible for a track flying parallel to the plane
  if (std::abs(tanL) < 1e-12)
    return c;
  const double s = (zPlane - zRef) / tanL;
  const double ds = dir * (s - sFrom); // distance travelled in the direction dir
  if (ds < sEps || ds > sMax)
    return c;
  c.ok = true;
  c.s = s;
  c.pos = at(s);
  return c;
}

Crossing Helix::crossLayer(const LayerDescriptor& layer, double sFrom, int dir, double tol, double sMax) const {
  Crossing c;
  if (layer.nSurfaces == 0)
    return c;
  if (layer.isBarrel()) {
    // barrel layer: the cylinder at its mean radius, accepted inside its z extent (+- tol)
    c = crossCylinder(layer.rMean_cm, sFrom, dir, 1e-6, sMax);
    if (c.ok && (c.pos.z() < layer.zMin_cm - tol || c.pos.z() > layer.zMax_cm + tol))
      c.ok = false;
  } else {
    // disk: the plane at its mean z, accepted inside its radial extent (+- tol)
    c = crossPlaneZ(layer.zMean_cm, sFrom, dir, 1e-6, sMax);
    if (c.ok) {
      const double r = std::hypot(c.pos.x(), c.pos.y());
      if (r < layer.rMin_cm - tol || r > layer.rMax_cm + tol)
        c.ok = false;
    }
  }
  return c;
}

Crossing Helix::refineToPlane(const Eigen::Vector3d& o, const Eigen::Vector3d& n, double sStart, int maxIter,
                              double tol) const {
  Crossing c;
  // Newton's method on f(s) = n . (at(s) - o), the signed distance of the helix point from the
  // plane: s <- s - f(s) / f'(s), with f'(s) = n . d(at)/ds.
  double s = sStart;
  for (int it = 0; it < maxIter; ++it) {
    const Eigen::Vector3d p = at(s);
    const double f = n.dot(p - o);
    if (std::abs(f) < tol) {
      c.ok = true;
      c.s = s;
      c.pos = p;
      return c;
    }
    // d/ds of the position, NOT normalised: (dx/ds, dy/ds) has unit length, dz/ds = tanL
    const double a = angleAt(s);
    const Eigen::Vector3d dp(-rot * std::sin(a), rot * std::cos(a), tanL);
    const double fp = n.dot(dp);
    if (std::abs(fp) < 1e-12)
      return c; // travelling parallel to the plane
    s -= f / fp;
  }
  return c;
}

} // namespace k4rt
