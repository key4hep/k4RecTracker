/**
 * Unit test of the kink chi2 of CombinatorialTrackFinder (TrackingMath).
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

// No Gaudi, no DD4hep, no input files.
//
// Toy calibration of the kink chi2: muons cross 9 barrel layers, with Highland multiple
// scattering in each layer and Gaussian measurement noise on each hit. When the kink chi2 is
// computed with the true circle radius, each kink (one triplet of hits) has 2 degrees of freedom,
// so the mean of chi2 / (2 (n - 2)) must be 1 (within 5%) at every momentum.
#include "HelixExtrapolator.h"
#include "TrackingMath.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

namespace {

// detector of the toy: 9 barrel layers (radius in cm, as in the IDEA vertex + inner tracker)
const std::vector<double> kLayerRadiiCm = {1.4, 2.4, 3.5, 13.0, 40.0, 67.0, 100.0, 157.0, 204.0};
const double kLayerThicknessX0 = 0.01;          // material of each layer, in radiation lengths
const double kHitResolutionCm = 5e-4;           // measurement error of each hit, both directions (5 um)
const double kBFieldTesla = 2.0;                // uniform field along z
const double kPolarAngle = 70.0 * M_PI / 180.0; // all muons fly at 70 degrees to the beam
const int kTracksPerMomentum = 2000;

// The random numbers of the toy. One generator with a fixed seed, so the test is reproducible.
struct RandomNumbers {
  std::mt19937 generator{7};
  std::normal_distribution<double> gaussian{0, 1};
  std::uniform_real_distribution<double> azimuth{0, 2 * M_PI};
};

// The helix of a particle at startPosition flying along direction with transverse momentum pT.
k4rt::Helix helixFromDirection(const Eigen::Vector3d& startPosition, const Eigen::Vector3d& direction, double pT,
                               double charge) {
  const double transverseLength = std::hypot(direction.x(), direction.y());
  // transverse direction of flight (unit vector)
  const double directionX = direction.x() / transverseLength;
  const double directionY = direction.y() / transverseLength;

  k4rt::Helix helix;
  helix.radius = k4rt::radiusFromPt(kBFieldTesla, pT);
  // a positive charge in a field along +z turns clockwise (rot = -1)
  if (charge * kBFieldTesla > 0) {
    helix.rotation = -1;
  } else {
    helix.rotation = +1;
  }
  // the circle centre is one radius away from the start, perpendicular to the direction
  helix.centreX = startPosition.x() + helix.rotation * helix.radius * (-directionY);
  helix.centreY = startPosition.y() + helix.rotation * helix.radius * directionX;
  helix.referenceAngle = std::atan2(startPosition.y() - helix.centreY, startPosition.x() - helix.centreX);
  helix.referenceZ = startPosition.z();
  helix.tanLambda = direction.z() / transverseLength;
  return helix;
}

// Simulate one muon from the origin through the layers. In each layer: move along the helix to
// the layer, store a hit with measurement noise, then turn the direction by a random scattering
// angle (Highland). Returns the hits as the kink chi2 sees them.
std::vector<k4rt::KinkPoint> simulateTrack(const Eigen::Vector3d& startDirection, double momentum, double charge,
                                           RandomNumbers& random) {
  const double pT = momentum * std::sin(kPolarAngle);
  Eigen::Vector3d position(0, 0, 0);
  Eigen::Vector3d direction = startDirection;
  std::vector<k4rt::KinkPoint> hits;

  for (double layerRadius : kLayerRadiiCm) {
    // 1. move to the next layer
    const k4rt::Helix helix = helixFromDirection(position, direction, pT, charge);
    const k4rt::Crossing crossing = helix.crossCylinder(layerRadius, 0.0, +1);
    if (!crossing.ok) {
      break; // the muon curls up before this layer
    }
    position = crossing.position;
    direction = helix.directionAt(crossing.arcLength);

    // 2. the measured hit: the true position moved by Gaussian noise along the two measurement
    // directions of a barrel sensor (along phi and along z)
    k4rt::KinkPoint hit;
    hit.normal = Eigen::Vector3d(position.x(), position.y(), 0).normalized();
    hit.xOverX0 = kLayerThicknessX0;
    hit.sigmaTransverseCm = kHitResolutionCm;
    hit.sigmaZCm = kHitResolutionCm;
    const Eigen::Vector3d alongPhi(-hit.normal.y(), hit.normal.x(), 0);
    const double noiseAlongPhi = kHitResolutionCm * random.gaussian(random.generator);
    const double noiseAlongZ = kHitResolutionCm * random.gaussian(random.generator);
    hit.position = position + noiseAlongPhi * alongPhi + noiseAlongZ * Eigen::Vector3d::UnitZ();
    hits.push_back(hit);

    // 3. multiple scattering in the layer: a random kick in two planes perpendicular to the
    // direction; the material seen grows as 1 / cos(incidence angle)
    const double cosIncidence = std::abs(hit.normal.dot(direction));
    const double scatteringAngle = k4rt::KinkChi2Calculator::highlandTheta0(kLayerThicknessX0 / cosIncidence, momentum);
    const Eigen::Vector3d kickDirection1 = direction.cross(Eigen::Vector3d::UnitZ()).normalized();
    const Eigen::Vector3d kickDirection2 = direction.cross(kickDirection1).normalized();
    const double kick1 = scatteringAngle * random.gaussian(random.generator);
    const double kick2 = scatteringAngle * random.gaussian(random.generator);
    direction = (direction + kick1 * kickDirection1 + kick2 * kickDirection2).normalized();
  }
  return hits;
}

} // namespace

int main() {
  RandomNumbers random;
  const k4rt::KinkChi2Calculator kinkCalculator(kBFieldTesla);
  int numberOfFailures = 0;

  std::printf("  p [GeV]   tracks   <chi2>/(2(n-2))   expected 1.00 +- 0.05\n");
  const std::vector<double> momenta = {0.3, 0.5, 1.0, 2.0, 5.0, 20.0};
  for (double momentum : momenta) {
    double sumOfNormalisedChi2 = 0;
    int numberOfTracks = 0;

    for (int trackIndex = 0; trackIndex < kTracksPerMomentum; ++trackIndex) {
      // a muon from the origin at a random azimuth; the charge alternates
      const double azimuth = random.azimuth(random.generator);
      double charge = -1.0;
      if (trackIndex % 2 == 1) {
        charge = 1.0;
      }
      const Eigen::Vector3d startDirection(std::sin(kPolarAngle) * std::cos(azimuth),
                                           std::sin(kPolarAngle) * std::sin(azimuth), std::cos(kPolarAngle));

      const std::vector<k4rt::KinkPoint> hits = simulateTrack(startDirection, momentum, charge, random);
      if (hits.size() < 4) {
        continue; // too few layers reached
      }

      // the kink chi2 with the true circle (the one at the start, before any scattering)
      const double pT = momentum * std::sin(kPolarAngle);
      const k4rt::Helix trueHelix = helixFromDirection(Eigen::Vector3d(0, 0, 0), startDirection, pT, charge);
      const k4rt::KinkResult result = kinkCalculator.track(hits, trueHelix.radius, trueHelix.rotation);
      const double numberOfDegreesOfFreedom = 2.0 * (hits.size() - 2);
      sumOfNormalisedChi2 += result.chi2 / numberOfDegreesOfFreedom;
      numberOfTracks += 1;
    }

    const double meanNormalisedChi2 = sumOfNormalisedChi2 / numberOfTracks;
    const bool passed = std::abs(meanNormalisedChi2 - 1.0) < 0.05;
    const char* verdict = "ok";
    if (!passed) {
      verdict = "FAIL";
      numberOfFailures += 1;
    }
    std::printf("  %6.1f   %6d   %8.3f            %s\n", momentum, numberOfTracks, meanNormalisedChi2, verdict);
  }

  if (numberOfFailures > 0) {
    std::printf("FAILED: %d failure(s)\n", numberOfFailures);
    return 1;
  }
  std::printf("PASSED: %d failure(s)\n", numberOfFailures);
  return 0;
}
