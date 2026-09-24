#pragma once

/** @file ClusterCounting.h
 *
 * Utilities for a parametrised cluster counting (dN/dx) simulation in gaseous tracking detectors:
 *
 *  - ClusterCounting::Parametrisation: number of primary ionisation clusters per unit length as a function of the
 *    beta*gamma of the particle, for a given gas mixture. The built-in tables are the ones of the Delphes
 *    TrackCovariance module (TrkUtil::Nclusters), interpolated with the same cubic spline, so that the numbers are
 *    identical to what was obtained when linking against Delphes.
 *  - ClusterCounting::trackLengthInCylinder: path length of a helix inside a cylindrical gas volume, for the first
 *    traversal of the volume by the track.
 *
 * All lengths are in mm (curvature in 1/mm) and all cluster densities are in clusters/mm.
 */

// ROOT
#include "TSpline.h"

// STL
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace k4RecTracker::ClusterCounting {

/// Gas mixtures for which built-in cluster density tables are available. The integer values match the "GasType" /
/// "GasSel" properties of the algorithms (and the options of the Delphes TrkUtil class).
enum class GasMixture : int {
  HeIsobutane_90_10 = 0, ///< He(90%)-Isobutane(10%)
  He = 1,                ///< pure He
  ArEthane_50_50 = 2,    ///< Ar(50%)-Ethane(50%)
  Ar = 3,                ///< pure Ar
};

/// Convert an integer gas selection to a GasMixture, returns std::nullopt if the value is not a known gas mixture
inline std::optional<GasMixture> toGasMixture(int gas) {
  if (gas < static_cast<int>(GasMixture::HeIsobutane_90_10) || gas > static_cast<int>(GasMixture::Ar)) {
    return std::nullopt;
  }
  return static_cast<GasMixture>(gas);
}

/** Number of primary ionisation clusters per mm as a function of beta*gamma.
 *
 * The data points are interpolated with a cubic spline (ROOT TSpline3). The spline is built once at construction, so
 * evaluating the parametrisation is cheap and does not allocate. Outside of the tabulated beta*gamma range the value
 * at the closest edge of the table is returned (the high beta*gamma end is the Fermi plateau); use inRange() to check
 * whether a value is covered by the table.
 *
 * Known limitation of the built-in tables (inherited from Delphes): the spline undershoots between the last two data
 * points (beta*gamma = 1000 and 10000), e.g. by ~27% at beta*gamma = 8000 for He-Isobutane, although both points are
 * on the plateau.
 */
class Parametrisation {
public:
  /// Build from a table of beta*gamma values (strictly increasing) and the corresponding number of clusters per mm
  Parametrisation(std::vector<double> betaGamma, std::vector<double> clustersPerMM)
      : m_betaGamma(std::move(betaGamma)), m_clustersPerMM(std::move(clustersPerMM)) {
    if (m_betaGamma.size() != m_clustersPerMM.size()) {
      throw std::invalid_argument("ClusterCounting::Parametrisation: beta*gamma and cluster density tables have "
                                  "different sizes");
    }
    if (m_betaGamma.size() < 3) {
      throw std::invalid_argument("ClusterCounting::Parametrisation: at least 3 data points are needed");
    }
    if (!std::is_sorted(m_betaGamma.begin(), m_betaGamma.end(), std::less_equal<>{})) {
      throw std::invalid_argument("ClusterCounting::Parametrisation: beta*gamma values must be strictly increasing");
    }
    m_spline = TSpline3("ClusterCountingParametrisation", m_betaGamma.data(), m_clustersPerMM.data(),
                        static_cast<int>(m_betaGamma.size()));
  }

  /// Built-in table for one of the supported gas mixtures
  static Parametrisation forGas(GasMixture gas) {
    // beta*gamma values of the data points
    std::vector<double> betaGamma = {0.5,  0.8,  1.,   2.,   3.,   4.,    5.,    8.,     10.,
                                     12.0, 15.0, 20.0, 50.0, 100., 200.0, 500.0, 1000.0, 10000.0};
    // Number of clusters per cm, from Delphes TrackCovariance TrkUtil::Nclusters (the last point is repeated to
    // smoothen the spline on the plateau)
    std::vector<double> clustersPerCM;
    switch (gas) {
    case GasMixture::HeIsobutane_90_10:
      clustersPerCM = {42.94, 23.6,  18.97, 12.98, 12.2,  12.13, 12.24, 12.73, 13.03,
                       13.29, 13.63, 14.08, 15.56, 16.43, 16.8,  16.95, 16.98, 16.98};
      break;
    case GasMixture::He:
      clustersPerCM = {11.79, 6.5, 5.23, 3.59, 3.38, 3.37, 3.4,  3.54, 3.63,
                       3.7,   3.8, 3.92, 4.33, 4.61, 4.78, 4.87, 4.89, 4.89};
      break;
    case GasMixture::ArEthane_50_50:
      clustersPerCM = {130.04, 71.55, 57.56, 39.44, 37.08, 36.9,  37.25, 38.76, 39.68,
                       40.49,  41.53, 42.91, 46.8,  48.09, 48.59, 48.85, 48.93, 48.93};
      break;
    case GasMixture::Ar:
      clustersPerCM = {88.69, 48.93, 39.41, 27.09, 25.51, 25.43, 25.69, 26.78, 27.44,
                       28.02, 28.77, 29.78, 32.67, 33.75, 34.24, 34.57, 34.68, 34.68};
      break;
    default:
      throw std::invalid_argument("ClusterCounting::Parametrisation: unknown gas mixture " +
                                  std::to_string(static_cast<int>(gas)));
    }
    for (auto& value : clustersPerCM) {
      value /= 10.0; // clusters/cm -> clusters/mm
    }
    return Parametrisation(std::move(betaGamma), std::move(clustersPerCM));
  }

  /// Number of clusters per mm for the given beta*gamma, clamped to the tabulated range
  double clustersPerMM(double betaGamma) const {
    return m_spline.Eval(std::clamp(betaGamma, minBetaGamma(), maxBetaGamma()));
  }

  double minBetaGamma() const { return m_betaGamma.front(); }
  double maxBetaGamma() const { return m_betaGamma.back(); }
  bool inRange(double betaGamma) const { return betaGamma >= minBetaGamma() && betaGamma <= maxBetaGamma(); }

private:
  std::vector<double> m_betaGamma;
  std::vector<double> m_clustersPerMM;
  TSpline3 m_spline;
};

/// Cylindrical volume, coaxial with the z axis, with inner radius rMin, outer radius rMax, and extending in z between
/// zMin and zMax (all in mm)
struct Cylinder {
  double rMin;
  double rMax;
  double zMin;
  double zMax;
};

/** Path length (mm) of a helix inside a cylindrical volume, for the first traversal of the volume.
 *
 * The helix is given by its (EDM4hep / LCIO convention) track parameters with respect to the origin: transverse impact
 * parameter d0 (mm), signed curvature omega (1/mm), longitudinal impact parameter z0 (mm) and tanLambda. The azimuthal
 * angle phi is not needed since the volume is rotationally symmetric. The track is followed from the point of closest
 * approach in the direction of motion, and the length of the first segment inside the volume is returned (so for
 * loopers only the first passage is counted).
 *
 * Returns 0 if the track never enters the volume, if the parameters do not describe a valid point of closest approach
 * (omega * d0 >= 1), or if the track never leaves the volume (e.g. a looper with tanLambda = 0 confined inside it).
 */
inline double trackLengthInCylinder(double d0, double omega, double z0, double tanLambda, const Cylinder& volume) {
  constexpr double infinity = std::numeric_limits<double>::infinity();
  constexpr double pi = std::numbers::pi;
  // Below this curvature (i.e. radius of curvature > 1000 km) the track is treated as a straight line
  constexpr double minAbsOmega = 1e-9;
  const double absOmega = std::abs(omega);
  const bool isStraight = absOmega < minAbsOmega;

  // With s the arc length in the transverse plane (s = 0 at the point of closest approach):
  //   r^2(s) = d0^2 + (1 - omega d0) * (2 sin(omega s / 2) / omega)^2,   z(s) = z0 + tanLambda * s
  const double b = 1.0 - omega * d0;
  if (!(b > 0.0)) {
    return 0.0;
  }
  auto r2 = [&](double s) {
    const double chord = isStraight ? s : 2.0 * std::sin(0.5 * absOmega * s) / absOmega;
    return d0 * d0 + b * chord * chord;
  };
  auto z = [&](double s) { return z0 + tanLambda * s; };
  auto isInside = [&](double s) {
    const double r2s = r2(s);
    const double zs = z(s);
    return r2s > volume.rMin * volume.rMin && r2s < volume.rMax * volume.rMax && zs > volume.zMin && zs < volume.zMax;
  };

  // Crossings with the cylinders at s > 0 within one turn (period) of the helix: for a circle there are up to two
  // (outgoing and incoming), repeating every turn; for a straight line there is at most one.
  const double period = isStraight ? infinity : 2.0 * pi / absOmega;
  std::vector<double> radialCrossings;
  for (const double radius : {volume.rMin, volume.rMax}) {
    const double a = radius * radius - d0 * d0;
    if (!(a > 0.0)) {
      continue; // radius <= |d0|: the track is always outside this cylinder
    }
    if (isStraight) {
      radialCrossings.push_back(std::sqrt(a / b));
      continue;
    }
    const double sinHalfPhase2 = omega * omega * a / (4.0 * b);
    if (sinHalfPhase2 > 1.0) {
      continue; // the track never reaches this radius
    }
    const double halfPhase = std::asin(std::sqrt(sinHalfPhase2));
    radialCrossings.push_back(2.0 * halfPhase / absOmega);
    radialCrossings.push_back(2.0 * (pi - halfPhase) / absOmega);
  }

  // First crossing of any boundary strictly after sFrom
  auto nextCrossing = [&](double sFrom) {
    double next = infinity;
    for (const double crossing : radialCrossings) {
      double s = crossing;
      if (s <= sFrom && std::isfinite(period)) {
        s += (std::floor((sFrom - crossing) / period) + 1.0) * period;
        if (s <= sFrom) {
          s += period;
        }
      }
      if (s > sFrom) {
        next = std::min(next, s);
      }
    }
    if (tanLambda != 0.0) {
      for (const double zWall : {volume.zMin, volume.zMax}) {
        const double s = (zWall - z0) / tanLambda;
        if (s > sFrom) {
          next = std::min(next, s);
        }
      }
    }
    return next;
  };

  // Walk along the track through the intervals between consecutive boundary crossings until the first one inside the
  // volume. A track crosses at most 6 boundaries (4 radial, 2 walls) before completing its first traversal of the
  // volume, the limit below only protects against numerical corner cases.
  constexpr int maxSteps = 16;
  double start = 0.0;
  for (int step = 0; step < maxSteps; ++step) {
    const double end = nextCrossing(start);
    if (!std::isfinite(end)) {
      return 0.0;
    }
    if (isInside(0.5 * (start + end))) {
      return (end - start) * std::sqrt(1.0 + tanLambda * tanLambda);
    }
    start = end;
  }
  return 0.0;
}

} // namespace k4RecTracker::ClusterCounting
