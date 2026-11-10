/**
 * Unit test of the helix extrapolation of CombinatorialTrackFinder (HelixExtrapolator).
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

// No Gaudi, no DD4hep, no input files.
//
// Every crossing computed by the helix code is compared with an independent reference: a
// 4th-order Runge-Kutta integration of the Lorentz force in a uniform field, started from the
// same position, momentum, charge and field. The two must agree to better than 10 um.
#include "HelixExtrapolator.h"
#include "TrackingMath.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

// pT [GeV] = kPtPerTeslaCm * B [T] * R [cm]
const double kPtPerTeslaCm = k4rt::kPtPerTeslaMetre / 100.0;

// largest allowed distance between the helix crossing and the Runge-Kutta one: 1e-3 cm = 10 um
const double kToleranceCm = 1e-3;

int numberOfFailures = 0;

// Print the result of one check, and count it if it failed.
void printCheck(bool passed, const std::string& description, double value, double limit) {
  const char* verdict = "ok";
  if (!passed) {
    verdict = "FAIL";
    numberOfFailures += 1;
  }
  std::printf("  %-58s %12.3e  (limit %8.1e)  %s\n", description.c_str(), value, limit, verdict);
}

// One simulated particle of the test.
struct TestParticle {
  std::string name;
  Eigen::Vector3d startPosition; // cm
  Eigen::Vector3d startMomentum; // GeV
  double charge;                 // +1 or -1
  double bFieldTesla;            // field along z
};

// ---------------------------------------------------------------------------------------------
// The reference: Runge-Kutta integration of the Lorentz force
// ---------------------------------------------------------------------------------------------

// A point of the particle's path: its position (cm) and momentum (GeV).
struct ParticleState {
  Eigen::Vector3d position;
  Eigen::Vector3d momentum;
};

// How the state changes per cm of 3D path length l, in a field B along +z:
//   d(position) / dl = the unit direction t
//   d(momentum) / dl = kPtPerTeslaCm * charge * B * (t x z)
ParticleState pathDerivative(const ParticleState& state, double charge, double bFieldTesla) {
  const Eigen::Vector3d direction = state.momentum.normalized();
  const Eigen::Vector3d zAxis(0, 0, 1);
  ParticleState derivative;
  derivative.position = direction;
  derivative.momentum = kPtPerTeslaCm * charge * bFieldTesla * direction.cross(zAxis);
  return derivative;
}

// state + factor * derivative (one Runge-Kutta sub-step)
ParticleState addStep(const ParticleState& state, const ParticleState& derivative, double factor) {
  ParticleState result;
  result.position = state.position + factor * derivative.position;
  result.momentum = state.momentum + factor * derivative.momentum;
  return result;
}

// One 4th-order Runge-Kutta step of stepLength cm along the path.
ParticleState rungeKuttaStep(const ParticleState& state, double stepLength, double charge, double bFieldTesla) {
  const ParticleState k1 = pathDerivative(state, charge, bFieldTesla);
  const ParticleState k2 = pathDerivative(addStep(state, k1, stepLength / 2), charge, bFieldTesla);
  const ParticleState k3 = pathDerivative(addStep(state, k2, stepLength / 2), charge, bFieldTesla);
  const ParticleState k4 = pathDerivative(addStep(state, k3, stepLength), charge, bFieldTesla);

  // weighted average of the four slopes
  ParticleState next;
  next.position = state.position + stepLength / 6 * (k1.position + 2 * k2.position + 2 * k3.position + k4.position);
  next.momentum = state.momentum + stepLength / 6 * (k1.momentum + 2 * k2.momentum + 2 * k3.momentum + k4.momentum);
  return next;
}

// A surface the particle crosses: a cylinder around z, a disk (plane at fixed z) or a tilted plane.
struct Surface {
  std::string type;            // "cylinder", "disk" or "plane"
  double radius = 0;           // cylinder: radius (cm)
  double z = 0;                // disk: z position (cm)
  Eigen::Vector3d planePoint;  // plane: one point on it
  Eigen::Vector3d planeNormal; // plane: its unit normal
};

// Signed distance of a point from the surface: negative on one side, positive on the other.
double signedDistance(const Surface& surface, const Eigen::Vector3d& point) {
  if (surface.type == "cylinder") {
    return std::hypot(point.x(), point.y()) - surface.radius;
  }
  if (surface.type == "disk") {
    return point.z() - surface.z;
  }
  return surface.planeNormal.dot(point - surface.planePoint);
}

// Follow the particle with Runge-Kutta steps of 0.05 cm until it crosses the surface, then find
// the crossing point by halving the last step. direction = +1 follows the particle forward and
// -1 backward (going backward = going forward with the momentum and the charge reversed).
// Returns NaN if the surface is not reached within 200 m of path.
Eigen::Vector3d rungeKuttaCrossing(ParticleState state, double charge, double bFieldTesla, int direction,
                                   const Surface& surface) {
  if (direction < 0) {
    state.momentum = -state.momentum;
    charge = -charge;
  }
  const double stepLength = 0.05; // cm
  double distanceBefore = signedDistance(surface, state.position);

  for (int step = 0; step < 400000; ++step) {
    const ParticleState next = rungeKuttaStep(state, stepLength, charge, bFieldTesla);
    const double distanceAfter = signedDistance(surface, next.position);
    const bool insideBefore = distanceBefore < 0;
    const bool insideAfter = distanceAfter < 0;

    if (insideBefore != insideAfter) {
      // the surface is crossed during this step: halve the step 80 times to find where
      double shortLength = 0;
      double longLength = stepLength;
      for (int halving = 0; halving < 80; ++halving) {
        const double middleLength = 0.5 * (shortLength + longLength);
        const ParticleState middle = rungeKuttaStep(state, middleLength, charge, bFieldTesla);
        const bool insideMiddle = signedDistance(surface, middle.position) < 0;
        if (insideMiddle == insideBefore) {
          shortLength = middleLength; // not crossed yet
        } else {
          longLength = middleLength; // already crossed
        }
      }
      return rungeKuttaStep(state, 0.5 * (shortLength + longLength), charge, bFieldTesla).position;
    }

    state = next;
    distanceBefore = distanceAfter;
  }
  return Eigen::Vector3d(NAN, NAN, NAN);
}

// ---------------------------------------------------------------------------------------------
// The helix of a particle, and the checks
// ---------------------------------------------------------------------------------------------

// The helix of a particle (charge +-1) at startPosition with startMomentum, in a field along z.
k4rt::Helix helixFromParticle(const TestParticle& particle) {
  const Eigen::Vector3d& position = particle.startPosition;
  const Eigen::Vector3d& momentum = particle.startMomentum;
  const double pT = std::hypot(momentum.x(), momentum.y());
  // transverse direction of flight (unit vector)
  const double directionX = momentum.x() / pT;
  const double directionY = momentum.y() / pT;

  k4rt::Helix helix;
  helix.R = k4rt::radiusFromPt(particle.bFieldTesla, pT);
  // a positive charge in a field along +z turns clockwise (rot = -1)
  if (particle.charge * particle.bFieldTesla > 0) {
    helix.rot = -1;
  } else {
    helix.rot = +1;
  }
  // the circle centre is one radius away from the start, perpendicular to the direction
  helix.cx = position.x() + helix.rot * helix.R * (-directionY);
  helix.cy = position.y() + helix.rot * helix.R * directionX;
  helix.phiRef = std::atan2(position.y() - helix.cy, position.x() - helix.cx);
  helix.zRef = position.z();
  helix.tanL = momentum.z() / pT;
  return helix;
}

// 1. The helix starts at the particle position and points along its momentum.
void checkStartPoint(const k4rt::Helix& helix, const TestParticle& particle) {
  const double positionError = (helix.at(0) - particle.startPosition).norm();
  printCheck(positionError < 1e-12, "helix passes through x0", positionError, 1e-12);

  const double directionError = (helix.dirAt(0) - particle.startMomentum.normalized()).norm();
  printCheck(directionError < 1e-12, "direction at x0 = p0 direction", directionError, 1e-12);
}

// 2. Forward crossing of a cylinder (a barrel layer) vs Runge-Kutta. The cylinder radius is half
// way between the start radius and the largest radius the circle reaches, so it is always
// crossed. Returns the crossing, which check 6 uses.
k4rt::Crossing checkCylinderCrossing(const k4rt::Helix& helix, const TestParticle& particle,
                                     const ParticleState& startState) {
  const double startRadius = std::hypot(particle.startPosition.x(), particle.startPosition.y());
  const double largestRadius = std::hypot(helix.cx, helix.cy) + helix.R;
  const double cylinderRadius = startRadius + 0.5 * (largestRadius - startRadius);

  const k4rt::Crossing crossing = helix.crossCylinder(cylinderRadius, 0.0, +1);
  double notFound = 1;
  if (crossing.ok) {
    notFound = 0;
  }
  printCheck(crossing.ok, "forward cylinder crossing found", notFound, 0.5);
  if (!crossing.ok) {
    return crossing;
  }

  const double distanceFromCylinder = std::abs(std::hypot(crossing.pos.x(), crossing.pos.y()) - cylinderRadius);
  printCheck(distanceFromCylinder < 1e-9, "crossing lies on the cylinder", distanceFromCylinder, 1e-9);

  Surface cylinder;
  cylinder.type = "cylinder";
  cylinder.radius = cylinderRadius;
  const Eigen::Vector3d reference = rungeKuttaCrossing(startState, particle.charge, particle.bFieldTesla, +1, cylinder);
  const double difference = (crossing.pos - reference).norm();
  printCheck(difference < kToleranceCm, "forward cylinder vs RK4", difference, kToleranceCm);
  return crossing;
}

// 3. Forward crossing of a disk (an endcap layer) 20 cm away in z vs Runge-Kutta.
void checkDiskCrossing(const k4rt::Helix& helix, const TestParticle& particle, const ParticleState& startState) {
  double diskZ = particle.startPosition.z() - 20.0;
  if (particle.startMomentum.z() > 0) {
    diskZ = particle.startPosition.z() + 20.0;
  }

  const k4rt::Crossing crossing = helix.crossPlaneZ(diskZ, 0.0, +1);
  double notFound = 1;
  if (crossing.ok) {
    notFound = 0;
  }
  printCheck(crossing.ok, "forward disk crossing found", notFound, 0.5);
  if (!crossing.ok) {
    return;
  }

  Surface disk;
  disk.type = "disk";
  disk.z = diskZ;
  const Eigen::Vector3d reference = rungeKuttaCrossing(startState, particle.charge, particle.bFieldTesla, +1, disk);
  const double difference = (crossing.pos - reference).norm();
  printCheck(difference < kToleranceCm, "disk vs RK4", difference, kToleranceCm);
}

// 4. Backward crossing of a cylinder at half the start radius (only for a displaced start).
void checkBackwardCrossing(const k4rt::Helix& helix, const TestParticle& particle, const ParticleState& startState) {
  const double startRadius = std::hypot(particle.startPosition.x(), particle.startPosition.y());
  if (startRadius <= 1.0) {
    return; // starts (almost) on the axis: nothing to cross backward
  }
  const double cylinderRadius = 0.5 * startRadius;

  const k4rt::Crossing crossing = helix.crossCylinder(cylinderRadius, 0.0, -1);
  if (!crossing.ok) {
    std::printf("  (backward cylinder r=%.2f not reached - skipped)\n", cylinderRadius);
    return;
  }

  Surface cylinder;
  cylinder.type = "cylinder";
  cylinder.radius = cylinderRadius;
  const Eigen::Vector3d reference = rungeKuttaCrossing(startState, particle.charge, particle.bFieldTesla, -1, cylinder);
  const double difference = (crossing.pos - reference).norm();
  printCheck(difference < kToleranceCm, "backward cylinder vs RK4", difference, kToleranceCm);
  printCheck(crossing.s < 0, "backward crossing has s < 0", crossing.s, 0);
}

// 5. Rebuild the helix from 7 points taken on it (more than half a turn), knowing its circle:
// the rotation sense, tan(lambda), the z reference and the unwrapped arc lengths must come back.
void checkHelixFromHits(const k4rt::Helix& helix) {
  const double pointSpacing = 0.45 * helix.R; // arc length between two points (cm)
  std::vector<Eigen::Vector3d> points;
  for (int pointIndex = 0; pointIndex < 7; ++pointIndex) {
    points.push_back(helix.at(pointSpacing * pointIndex));
  }

  k4rt::Helix fittedHelix;
  std::vector<double> arcLengths;
  const bool fitWorked = k4rt::Helix::fromHits(helix.cx, helix.cy, helix.R, points, fittedHelix, &arcLengths);

  double rotationDifference = 9; // shown when the fit failed
  if (fitWorked) {
    rotationDifference = std::abs(fittedHelix.rot - helix.rot);
  }
  printCheck(fitWorked && fittedHelix.rot == helix.rot, "helixFromHits: rotation sense", rotationDifference, 0.5);

  const double tanLambdaDifference = std::abs(fittedHelix.tanL - helix.tanL);
  printCheck(tanLambdaDifference < 1e-9, "helixFromHits: tanLambda", tanLambdaDifference, 1e-9);

  const double zReferenceDifference = std::abs(fittedHelix.zRef - helix.zRef);
  printCheck(zReferenceDifference < 1e-9, "helixFromHits: zRef", zReferenceDifference, 1e-9);

  const double lastArcLengthDifference = std::abs(arcLengths.back() - pointSpacing * 6);
  printCheck(lastArcLengthDifference < 1e-9, "helixFromHits: unwrapped arc length > pi R", lastArcLengthDifference,
             1e-9);
}

// 6. Exact crossing of a tilted sensor plane near the cylinder crossing of check 2, vs Runge-Kutta.
void checkPlaneCrossing(const k4rt::Helix& helix, const TestParticle& particle, const ParticleState& startState,
                        const k4rt::Crossing& cylinderCrossing) {
  // the plane normal: the radial direction at the cylinder crossing, tilted a little
  const Eigen::Vector3d radialDirection =
      Eigen::Vector3d(cylinderCrossing.pos.x(), cylinderCrossing.pos.y(), 0).normalized();
  const Eigen::Vector3d tilt(0.2, -0.1, 0.15);
  const Eigen::Vector3d planeNormal = (radialDirection + tilt).normalized();
  // the plane goes through a point 0.3 cm sideways from the cylinder crossing
  const Eigen::Vector3d sideways(-planeNormal.y(), planeNormal.x(), 0);
  const Eigen::Vector3d planePoint = cylinderCrossing.pos + 0.3 * sideways;

  const k4rt::Crossing crossing = helix.refineToPlane(planePoint, planeNormal, cylinderCrossing.s);
  double notFound = 1;
  if (crossing.ok) {
    notFound = 0;
  }
  printCheck(crossing.ok, "refineToPlane converged", notFound, 0.5);
  if (!crossing.ok) {
    return;
  }

  const double distanceFromPlane = std::abs(planeNormal.dot(crossing.pos - planePoint));
  printCheck(distanceFromPlane < 1e-9, "point lies on the plane", distanceFromPlane, 1e-9);

  Surface plane;
  plane.type = "plane";
  plane.planePoint = planePoint;
  plane.planeNormal = planeNormal;
  const Eigen::Vector3d reference = rungeKuttaCrossing(startState, particle.charge, particle.bFieldTesla, +1, plane);
  const double difference = (crossing.pos - reference).norm();
  printCheck(difference < kToleranceCm, "plane intersection vs RK4", difference, kToleranceCm);
}

} // namespace

int main() {
  // prompt and displaced particles, both charges, both field directions, low and higher pT
  std::vector<TestParticle> particles;
  particles.push_back(
      {"prompt mu+  pT=1 GeV  B=+2T", Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0.8, 0.6, 0.5), +1, 2.0});
  particles.push_back(
      {"prompt mu-  pT=1 GeV  B=+2T", Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0.6, -0.8, -0.3), -1, 2.0});
  particles.push_back(
      {"displaced mu+ pT=3   B=-1.7", Eigen::Vector3d(10, 5, 3), Eigen::Vector3d(2.0, 2.2, 1.5), +1, -1.7});
  particles.push_back(
      {"low-pT mu-  pT=0.4  B=+2T", Eigen::Vector3d(1, -2, 0.5), Eigen::Vector3d(-0.3, 0.26, 0.2), -1, 2.0});

  for (const TestParticle& particle : particles) {
    std::printf("%s\n", particle.name.c_str());
    const k4rt::Helix helix = helixFromParticle(particle);
    ParticleState startState;
    startState.position = particle.startPosition;
    startState.momentum = particle.startMomentum;

    checkStartPoint(helix, particle);
    const k4rt::Crossing cylinderCrossing = checkCylinderCrossing(helix, particle, startState);
    checkDiskCrossing(helix, particle, startState);
    checkBackwardCrossing(helix, particle, startState);
    checkHelixFromHits(helix);
    if (cylinderCrossing.ok) {
      checkPlaneCrossing(helix, particle, startState, cylinderCrossing);
    }
  }

  if (numberOfFailures > 0) {
    std::printf("\nFAILED: %d failure(s)\n", numberOfFailures);
    return 1;
  }
  std::printf("\nPASSED: %d failure(s)\n", numberOfFailures);
  return 0;
}
