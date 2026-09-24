// Unit tests for the cluster counting utilities in k4RecTracker/ClusterCounting.h

#include "k4RecTracker/ClusterCounting.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>

using namespace k4RecTracker::ClusterCounting;

namespace {

int n_failures = 0;

void check_close(double value, double expected, const std::string& what, double tolerance = 1e-9) {
  if (std::abs(value - expected) > tolerance * std::max(1.0, std::abs(expected))) {
    std::cerr << "FAILED: " << what << ": got " << value << ", expected " << expected << std::endl;
    ++n_failures;
  }
}

void check(bool condition, const std::string& what) {
  if (!condition) {
    std::cerr << "FAILED: " << what << std::endl;
    ++n_failures;
  }
}

void test_parametrisation() {
  // The spline has to go through the data points (values from Delphes TrkUtil::Nclusters, in clusters/cm)
  const auto he_iso = Parametrisation::forGas(GasMixture::HeIsobutane_90_10);
  check_close(he_iso.clustersPerMM(0.5), 4.294, "He-Isobutane at bg = 0.5");
  check_close(he_iso.clustersPerMM(1.0), 1.897, "He-Isobutane at bg = 1");
  check_close(he_iso.clustersPerMM(100.0), 1.643, "He-Isobutane at bg = 100");
  check_close(Parametrisation::forGas(GasMixture::He).clustersPerMM(4.0), 0.337, "He at bg = 4");
  check_close(Parametrisation::forGas(GasMixture::ArEthane_50_50).clustersPerMM(20.0), 4.291, "Ar-Ethane at bg = 20");
  check_close(Parametrisation::forGas(GasMixture::Ar).clustersPerMM(1000.0), 3.468, "Ar at bg = 1000");

  // Outside of the table the value at the edge is returned
  check(he_iso.inRange(0.5) && he_iso.inRange(10000.0), "table edges are in range");
  check(!he_iso.inRange(0.49) && !he_iso.inRange(10001.0), "values outside the table are out of range");
  check_close(he_iso.clustersPerMM(0.1), he_iso.clustersPerMM(0.5), "clamping below the table");
  check_close(he_iso.clustersPerMM(20000.0), 1.698, "clamping above the table (plateau)");

  // Reference values of the Delphes implementation between data points (including the undershoot of the spline
  // between beta*gamma = 1000 and 10000), to make sure the parametrisation is unchanged
  check_close(he_iso.clustersPerMM(1.5), 1.406035318, "He-Isobutane at bg = 1.5", 1e-9);
  check_close(he_iso.clustersPerMM(8000.0), 1.238399472, "He-Isobutane at bg = 8000", 1e-9);

  // Gas selection
  check(toGasMixture(0) == GasMixture::HeIsobutane_90_10, "gas 0");
  check(toGasMixture(3) == GasMixture::Ar, "gas 3");
  check(!toGasMixture(-1).has_value() && !toGasMixture(4).has_value(), "invalid gas selections");

  // Custom tables
  bool thrown = false;
  try {
    Parametrisation({1.0, 0.5, 2.0}, {1.0, 1.0, 1.0});
  } catch (const std::invalid_argument&) {
    thrown = true;
  }
  check(thrown, "non increasing beta*gamma table is rejected");
  const Parametrisation flat({1.0, 2.0, 3.0}, {0.5, 0.5, 0.5});
  check_close(flat.clustersPerMM(2.5), 0.5, "custom flat table");
}

void test_track_length() {
  const Cylinder chamber{350.0, 2000.0, -2000.0, 2000.0};
  const double length_radial = chamber.rMax - chamber.rMin;

  // Straight track from the origin, perpendicular to the beam
  check_close(trackLengthInCylinder(0.0, 0.0, 0.0, 0.0, chamber), length_radial, "straight radial track");
  check_close(trackLengthInCylinder(0.0, 1e-12, 0.0, 0.0, chamber), length_radial, "almost straight radial track");

  // Straight track with polar angle: path length scales with sqrt(1 + tanLambda^2)
  check_close(trackLengthInCylinder(0.0, 0.0, 0.0, 0.5, chamber), length_radial * std::sqrt(1.25),
              "straight track with tanLambda = 0.5");

  // Curved track from the origin: s(r) = 2R asin(r / 2R) for radius of curvature R, independent of the charge
  const double radius = 3000.0;
  const double expected_curved =
      2.0 * radius * (std::asin(chamber.rMax / (2.0 * radius)) - std::asin(chamber.rMin / (2.0 * radius)));
  check_close(trackLengthInCylinder(0.0, 1.0 / radius, 0.0, 0.0, chamber), expected_curved, "positive curvature");
  check_close(trackLengthInCylinder(0.0, -1.0 / radius, 0.0, 0.0, chamber), expected_curved, "negative curvature");
  check_close(trackLengthInCylinder(0.0, 1.0 / radius, 0.0, -0.3, chamber), expected_curved * std::sqrt(1.09),
              "curved track with tanLambda = -0.3");

  // Forward track leaving through the endcap: enters at r = rMin, leaves at z = zMax
  const double tan_lambda = 3.0;
  check_close(trackLengthInCylinder(0.0, 0.0, 0.0, tan_lambda, chamber),
              (chamber.zMax / tan_lambda - chamber.rMin) * std::sqrt(1.0 + tan_lambda * tan_lambda),
              "track leaving through the endcap");
  check_close(trackLengthInCylinder(0.0, 0.0, 0.0, -tan_lambda, chamber),
              (-chamber.zMin / tan_lambda - chamber.rMin) * std::sqrt(1.0 + tan_lambda * tan_lambda),
              "track leaving through the negative endcap");

  // Very forward track never entering the chamber (leaves through the endcap plane before reaching rMin)
  check_close(trackLengthInCylinder(0.0, 0.0, 0.0, 10.0, chamber), 0.0, "track passing inside the inner cylinder");

  // Low momentum looper that does not reach rMax: enters and leaves through the inner cylinder
  const double looper_radius = 500.0;
  check_close(trackLengthInCylinder(0.0, 1.0 / looper_radius, 0.0, 0.0, chamber),
              2.0 * looper_radius * (std::numbers::pi - 2.0 * std::asin(chamber.rMin / (2.0 * looper_radius))),
              "looper");

  // Looper that never reaches the chamber
  check_close(trackLengthInCylinder(0.0, 1.0 / 150.0, 0.0, 0.0, chamber), 0.0, "looper inside the inner cylinder");

  // Track starting inside the chamber (large impact parameter), straight line
  const double d0 = 1000.0;
  check_close(trackLengthInCylinder(d0, 0.0, 0.0, 0.0, chamber), std::sqrt(chamber.rMax * chamber.rMax - d0 * d0),
              "track with large impact parameter");

  // Track starting outside of the chamber in z, moving away from it
  check_close(trackLengthInCylinder(0.0, 0.0, 2500.0, 1.0, chamber), 0.0, "track outside in z moving away");

  // Invalid point of closest approach
  check_close(trackLengthInCylinder(10.0, 0.2, 0.0, 0.0, chamber), 0.0, "invalid parameters (omega * d0 > 1)");
}

} // namespace

int main() {
  test_parametrisation();
  test_track_length();
  if (n_failures > 0) {
    std::cerr << n_failures << " check(s) failed" << std::endl;
    return 1;
  }
  std::cout << "All checks passed" << std::endl;
  return 0;
}
