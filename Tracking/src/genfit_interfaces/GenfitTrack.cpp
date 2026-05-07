/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "GenfitTrack.h"

namespace GenfitInterface {

GenfitTrack::GenfitTrack(const edm4hep::Track& track, const bool skipTrackOrdering,
                         const dd4hep::rec::DCH_info* dch_info, const dd4hep::DDSegmentation::BitFieldCoder* decoder,
                         const GenfitInterface::GenfitField* fieldMap)
    : m_originalTrack(track), m_posInit(0., 0., 0.), m_momInit(0., 0., 0.), m_covInit(6), m_genfitTrackRep(nullptr),
      m_genfitTrack(nullptr), m_edm4hepTrack(),

      m_dch_info(dch_info), m_dc_decoder(decoder), m_fieldMap(fieldMap)

{

  CheckInitialization();
  OrderHits(track, skipTrackOrdering);
}

GenfitTrack::~GenfitTrack() {}

/**
 * @brief Check if required Genfit components are properly initialized.
 *
 * This method verifies whether the singleton instances of `genfit::FieldManager`
 * and `genfit::MaterialEffects` have been initialized. These components are essential
 * for tracking operations in the Genfit framework.
 *
 * @note This method should be called before performing any tracking-related operations
 * that depend on magnetic field or material effects.
 */
void GenfitTrack::CheckInitialization() {

  if (!genfit::FieldManager::getInstance()->isInitialized()) {
    throw std::runtime_error("Error: FieldManager is not initialized!");
  }

  if (!genfit::MaterialEffects::getInstance()->isInitialized()) {
    throw std::runtime_error("Error: MaterialEffects is not initialized!");
  }
}

/**
 * @brief Orders the tracker hits of a track according to their spatial progression.
 *
 * This method determines a suitable starting point among the track hits,
 * typically the one closest to the interaction region (|r| minimum).
 *
 * All hits are then sorted according to their 3D distance from the chosen
 * starting point, ensuring a consistent ordering along the track trajectory.
 * The ordered hits are stored in the internal MutableTrack object for
 * subsequent processing (e.g. track fitting).
 *
 * @param track The input track whose tracker hits must be reordered.
 * @param skipTrackOrdering If true, skip the track ordering step.
 */
void GenfitTrack::OrderHits(const edm4hep::Track& track, bool skipTrackOrdering) {

  if (skipTrackOrdering) {

    for (const auto& hit : track.getTrackerHits()) {

      m_edm4hepTrack.addToTrackerHits(hit);
    }
    return;
  }

  const auto hits = track.getTrackerHits();
  std::vector<std::pair<float, std::size_t>> distIndex;

  // Track endpoints in cylindrical radius
  double rMin = std::numeric_limits<double>::max();
  double firstX = 0., firstY = 0., firstZ = 0.;

  // Loop over hits to find the hit with smallest r = sqrt(x^2 + y^2)
  for (const auto& hit : hits) {
    const auto p = hit.getPosition();
    double r = std::hypot(p.x, p.y);

    if (r < rMin) {
      rMin = r;
      firstX = p.x;
      firstY = p.y;
      firstZ = p.z;
    }
  }

  // Compute distances from the chosen starting point
  distIndex.reserve(hits.size());

  for (std::size_t i = 0; i < hits.size(); ++i) {
    const auto p = hits[i].getPosition();
    const float d = std::hypot(p.x - firstX, p.y - firstY, p.z - firstZ);

    distIndex.emplace_back(d, i);
  }

  // Sort hits along the track
  std::ranges::sort(distIndex, {}, &std::pair<float, std::size_t>::first);

  m_edm4hepTrack = edm4hep::MutableTrack();
  for (const auto& [_, idx] : distIndex) {
    m_edm4hepTrack.addToTrackerHits(hits[idx]);
  }
}

/**
 * @brief Initialize the GenfitTrack object with configurable initialization strategies.
 *
 * This method initializes the internal state of the `GenfitTrack` object by defining:
 * - the initial position (`m_posInit`)
 * - the initial momentum (`m_momInit`)
 * - the charge hypothesis (`m_charge_hypothesis`)
 * - the initial covariance matrix (`m_covInit`)
 *
 * The initialization procedure depends on the selected strategy and may rely on:
 * - tracker hits,
 * - pre-existing track states,
 * - or user-provided parameters.
 *
 * In all cases, the magnetic field at the reference position is evaluated and used
 * to ensure a consistent initialization of both kinematics and covariance.
 *
 * ---
 * Initialization workflow:
 *
 * 1. Optional hit reduction:
 *    - If `LimitHits` is true, the tracker hits are preprocessed using
 *      `LimitNumberHits(Epsilon, Window)`.
 *    - This step is intended to reduce the number of hits, retaining only those associated with the first loop of
 * looping tracks.
 *
 * 2. Definition of initial track parameters (controlled by InitializationType):
 *
 *    - InitializationType == 0 (default: geometric seed from hits):
 *        - Requires at least two tracker hits.
 *        - The initial position is set to the first hit.
 *        - The initial momentum direction is defined as the normalized vector
 *          from the first to the second hit.
 *        - The charge hypothesis is estimated via `ComputeInitialParameters(Bz)`.
 *        - The covariance matrix is computed using `ComputeInitialCovarianceMatrix`.
 *
 *    - InitializationType == 1 (refined: parameter-based estimation):
 *        - Uses `ComputeInitialParameters(Bz)` to obtain a more accurate estimate
 *          of position and momentum, typically exploiting multiple hits.
 *        - Optionally updates the reference point (`m_VP_referencePoint`):
 *            - if `UseFirstHitAsReference` is true, or
 *            - if displaced tracking is enabled and the first hit lies outside
 *              a sphere of radius `RadiusForDisplacedTracking`.
 *        - The covariance matrix is computed consistently with the estimated parameters.
 *
 *    - InitializationType == 2 (track state based):
 *        - Extracts initialization from an existing `edm4hep::TrackState`.
 *        - The reference point is taken from the track state.
 *        - The momentum is reconstructed from helix parameters:
 *            - transverse momentum is derived from curvature (`omega`) and Bz
 *            - direction is defined using `phi` and `tanLambda`
 *        - The covariance matrix is:
 *            - first reconstructed in helix parameter space
 *            - then transformed to Cartesian coordinates via
 *              `CovarianceMatrixHelixToCartesian`.
 *
 *    - InitializationType == 3 (custom user-defined):
 *        - Uses user-provided `Init_position` and `Init_momentum`.
 *        - Both parameters are mandatory.
 *        - The charge hypothesis is still inferred via `ComputeInitialParameters(Bz)`.
 *        - The covariance matrix is computed accordingly.
 *
 * ---
 *
 * Error handling:
 * - The method throws `std::runtime_error` in the following cases:
 *     - missing required parameters (e.g. TrackStateLocation, Init_position, ...)
 *     - insufficient number of hits for the selected strategy
 *     - requested track state not found
 *     - unknown InitializationType
 *
 * ---
 *
 * @param RadiusForDisplacedTracking If > 0, enables displaced tracking logic:
 *        tracks with first hit outside this radius may use a different reference point.
 * @param UseFirstHitAsReference If true, the first hit is always used as reference point.
 * @param LimitHits If true, reduces the number of hits before initialization.
 * @param InitializationType Strategy selector:
 *        0 = default (first two hits),
 *        1 = refined (computed parameters),
 *        2 = track state (from edm4hep::TrackState),
 *        3 = custom (user-defined).
 * @param TrackStateLocation Required if InitializationType == 2:
 *        0 = edm4hep::TrackState::AtOther
 *        1 = edm4hep::TrackState::AtFirstHit
 *        2 = edm4hep::TrackState::AtLastHit
 *        3 = edm4hep::TrackState::AtCalorimeter
 * @param Init_position Custom initial position (required if InitializationType == 3).
 * @param Init_momentum Custom initial momentum (required if InitializationType == 3).
 * @param Epsilon Parameter for hit reduction (required if LimitHits == true).
 * @param Window Parameter for hit reduction (required if LimitHits == true).
 * @param sigma_d0        Optional standard deviation for d0 (used in covariance construction).
 * @param sigma_phi       Optional standard deviation for phi (used in covariance construction).
 * @param omega_factor    Optional factor for omega (used in covariance construction).
 * @param z0_factor       Optional factor for z0 (used in covariance construction).
 * @param sigma_tanLambda Optional standard deviation for tanLambda (used in covariance construction).
 *
 *
 * @note The optional sigma parameters (sigma_d0, sigma_phi, omega_factor, z0_factor, sigma_tanLambda) are used to
 * construct the initial covariance matrix in the absence of a track state (InitializationType == 0,1,3).
 *
 * They are applied only in:
 * - InitializationType == 0 (default)
 * - InitializationType == 1 (refined initialization)
 * - InitializationType == 3 (custom initialization)
 *
 * If not provided, default physically motivated values are used internally.
 * These parameters are not used in InitializationType == 2, where the
 * covariance is directly taken from the input TrackState and transformed
 * to Cartesian coordinates.
 */
void GenfitTrack::InitializeTrack(double RadiusForDisplacedTracking, bool UseFirstHitAsReference, bool LimitHits,
                                  int InitializationType, std::optional<int> TrackStateLocation,
                                  std::optional<TVector3> Init_position, std::optional<TVector3> Init_momentum,
                                  std::optional<double> Epsilon, std::optional<int> Window,
                                  std::optional<double> sigma_d0, std::optional<double> sigma_phi,
                                  std::optional<double> omega_factor, std::optional<double> z0_factor,
                                  std::optional<double> sigma_tanLambda) {

  // -------------------------
  // Optional hit preprocessing
  // -------------------------
  if (LimitHits) {
    if (Epsilon.has_value() && Window.has_value()) {
      LimitNumberHits(Epsilon.value(), Window.value());
    } else {
      throw std::runtime_error("InitializeTrack: Epsilon and Window are required when LimitHits is true");
    }
  }

  // =========================
  // Initialization strategies
  // =========================

  if (InitializationType == 0) // --- Default: first two hits
  {
    auto hits = m_edm4hepTrack.getTrackerHits();

    if (hits.size() < 2) {
      throw std::runtime_error("InitializeTrack: Not enough hits for InitializationType == 0");
    }

    auto pos1 = hits[0].getPosition();
    auto pos2 = hits[1].getPosition();

    TVector3 p1(pos1.x * dd4hep::mm, pos1.y * dd4hep::mm, pos1.z * dd4hep::mm);

    TVector3 p2(pos2.x * dd4hep::mm, pos2.y * dd4hep::mm, pos2.z * dd4hep::mm);

    // Seed definition
    m_posInit = p1;
    m_momInit = (p2 - p1).Unit();

    double Bz = m_fieldMap->getBz(m_posInit) / (dd4hep::tesla / dd4hep::kilogauss); // From kilogauss to Tesla

    HelperInitialization initInfo = ComputeInitialParameters(Bz);
    m_charge_hypothesis = initInfo.Charge;

    m_covInit = ComputeInitialCovarianceMatrix(Bz, m_charge_hypothesis, sigma_d0.value(), sigma_phi.value(),
                                               omega_factor.value(), z0_factor.value(), sigma_tanLambda.value());
  }

  else if (InitializationType == 1) // --- Refined
  {
    auto hits = m_edm4hepTrack.getTrackerHits();

    if (hits.empty()) {
      throw std::runtime_error("InitializeTrack: No hits available for InitializationType == 1");
    }

    auto firstHit = hits[0].getPosition();

    TVector3 firstHitVec(firstHit.x * dd4hep::mm, firstHit.y * dd4hep::mm, firstHit.z * dd4hep::mm);

    double Bz = m_fieldMap->getBz(firstHitVec) / (dd4hep::tesla / dd4hep::kilogauss); // From kilogauss to Tesla

    // Optional reference point update
    if (UseFirstHitAsReference || RadiusForDisplacedTracking > 0.) {
      if (UseFirstHitAsReference || firstHitVec.Mag() > RadiusForDisplacedTracking) {
        m_VP_referencePoint = firstHitVec;
      }
    }

    HelperInitialization initInfo = ComputeInitialParameters(Bz);

    m_charge_hypothesis = initInfo.Charge;
    m_posInit = initInfo.Position;
    m_momInit = initInfo.Momentum;

    m_covInit = ComputeInitialCovarianceMatrix(Bz, m_charge_hypothesis, sigma_d0.value(), sigma_phi.value(),
                                               omega_factor.value(), z0_factor.value(), sigma_tanLambda.value());
  }

  else if (InitializationType == 2) // --- From TrackState
  {
    if (!TrackStateLocation.has_value()) {
      throw std::runtime_error("InitializeTrack: TrackStateLocation is required for InitializationType == 2");
    }

    auto trackStates = m_originalTrack.getTrackStates();
    bool found = false;

    for (auto ts : trackStates) {
      if (ts.location == TrackStateLocation.value()) {
        found = true;

        m_VP_referencePoint = TVector3(ts.referencePoint.x * dd4hep::mm, ts.referencePoint.y * dd4hep::mm,
                                       ts.referencePoint.z * dd4hep::mm); // in cm

        double Bz =
            m_fieldMap->getBz(m_VP_referencePoint) / (dd4hep::tesla / dd4hep::kilogauss); // From kilogauss to Tesla

        HelperInitialization initInfo = ComputeInitialParameters(Bz);
        m_posInit = initInfo.Position;
        m_charge_hypothesis = initInfo.Charge;

        // Helix -> momentum conversion
        double pT = ConversionUnits::a_lcio * std::abs(Bz) / std::abs(ts.omega);
        double px = pT * std::cos(ts.phi);
        double py = pT * std::sin(ts.phi);
        double pz = pT * ts.tanLambda;

        m_momInit = TVector3(px, py, pz);

        // Build symmetric covariance matrix (helix parameters)
        TMatrixDSym C_helix(5);
        C_helix.Zero();

        int index = 0;
        for (int i = 0; i < 5; ++i) {
          for (int j = 0; j <= i; ++j) {
            C_helix(i, j) = ts.covMatrix[index++];
          }
        }

        for (int i = 0; i < 5; ++i) {
          for (int j = i + 1; j < 5; ++j) {
            C_helix(i, j) = C_helix(j, i);
          }
        }

        // Conversion factors:
        // d0 : mm -> cm
        // phi0 : unchanged
        // omega : 1/mm -> 1/cm => x10
        // z0 : mm -> cm        => x0.1
        // tanLambda : unchanged

        double scale[5] = {dd4hep::mm, dd4hep::cm, 1. / dd4hep::mm, dd4hep::mm, dd4hep::cm};
        for (int i = 0; i < 5; ++i) {
          for (int j = 0; j < 5; ++j) {
            C_helix(i, j) *= scale[i] * scale[j];
          }
        }

        // Convert to Cartesian covariance
        m_covInit = CovarianceMatrixHelixToCartesian(C_helix, m_posInit, m_momInit, m_VP_referencePoint,
                                                     m_charge_hypothesis, Bz);
      }
    }

    if (!found) {
      throw std::runtime_error("InitializeTrack: Requested TrackStateLocation not found");
    }
  }

  else if (InitializationType == 3) // --- Custom
  {
    if (!(Init_position.has_value() && Init_momentum.has_value())) {
      throw std::runtime_error(
          "InitializeTrack: Init_position and Init_momentum are required for InitializationType == 3");
    }

    m_posInit = Init_position.value();
    m_momInit = Init_momentum.value();

    double Bz = m_fieldMap->getBz(m_posInit) / (dd4hep::tesla / dd4hep::kilogauss); // From kilogauss to Tesla

    HelperInitialization initInfo = ComputeInitialParameters(Bz);
    m_charge_hypothesis = initInfo.Charge;

    m_covInit = ComputeInitialCovarianceMatrix(Bz, m_charge_hypothesis, sigma_d0.value(), sigma_phi.value(),
                                               omega_factor.value(), z0_factor.value(), sigma_tanLambda.value());
  }

  else {
    throw std::runtime_error("InitializeTrack: Unknown InitializationType");
  }
}

/**
 * @brief Limits the number of tracker hits in the track based on a local extremum
 *        along the y-coordinate after smoothing.
 *
 * This method identifies the first significant local extremum (maximum or minimum)
 * in the y-coordinate of the track hits, after applying a moving average smoothing.
 * It then truncates the track to include only hits up to that extremum.
 *
 * Procedure:
 * 1. Extracts the z and y positions of all tracker hits.
 * 2. Applies a moving average smoothing to the y-coordinates using the specified window.
 * 3. Computes the first derivative of the smoothed y with respect to z.
 * 4. Finds the first index where the derivative changes sign and exceeds the
 *    specified threshold epsilon, indicating a local extremum.
 * 5. If such an extremum is found, rebuilds the track keeping only hits up to that index.
 *
 * @param epsilon       Threshold for detecting significant changes in the derivative.
 * @param smoothWindow  Size of the smoothing window (number of hits) for the moving average.
 */
void GenfitTrack::LimitNumberHits(double epsilon, int smoothWindow) {

  auto hits = m_edm4hepTrack.getTrackerHits();
  int maxHit = hits.size();
  int n = maxHit;

  if (maxHit == 0) {
    throw std::runtime_error("Internal edm4hep::Track is empty.");
  }

  if (n < smoothWindow || n < 3)
    maxHit = -1;
  int W = smoothWindow / 2;

  std::vector<double> z_val;
  std::vector<double> y_val_raw;
  z_val.reserve(hits.size());
  y_val_raw.reserve(hits.size());

  // Extract positions
  for (const auto& hit : hits) {
    auto pos = hit.getPosition();
    z_val.push_back(pos.z);
    y_val_raw.push_back(pos.y);
  }

  // Apply moving average smoothing
  std::vector<double> y_val;
  y_val.reserve(hits.size());
  for (std::size_t i = 0; i < hits.size(); ++i) {
    int start = (i < static_cast<std::size_t>(W)) ? 0 : i - W;
    int end = std::min(i + W + 1, hits.size());

    double sum = 0.0;
    for (int j = start; j < end; ++j)
      sum += y_val_raw[j];
    y_val.push_back(sum / (end - start));
  }

  // Compute first derivative
  std::vector<double> d(n, 0.0);

  double dz0 = z_val[1] - z_val[0];
  d[0] = (dz0 != 0.0) ? (y_val[1] - y_val[0]) / dz0 : 0.0;

  for (std::size_t i = 1; i < static_cast<std::size_t>(W - 1); ++i) {
    double dz = z_val[i + 1] - z_val[i - 1];
    d[i] = (dz != 0.0) ? (y_val[i + 1] - y_val[i - 1]) / dz : 0.0;
  }

  double dzn = z_val[n - 1] - z_val[n - 2];
  d[n - 1] = (dzn != 0.0) ? (y_val[n - 1] - y_val[n - 2]) / dzn : 0.0;

  // Find first significant extremum
  bool MaxFound = false;
  for (int i = 1; i < n - 1; ++i) {
    if ((d[i - 1] > epsilon && d[i] < -epsilon) || (d[i - 1] < -epsilon && d[i] > epsilon)) {
      maxHit = i;
      MaxFound = true;
      break;
    }
  }

  if (!MaxFound)
    maxHit = -1;

  // Rebuild track keeping only hits up to maxHit
  if (maxHit != -1) {
    auto temp_track = m_edm4hepTrack;
    m_edm4hepTrack = edm4hep::MutableTrack();

    int idx_fill_track = 0;
    auto temp_hits = temp_track.getTrackerHits();
    for (const auto& hit : temp_hits) {
      m_edm4hepTrack.addToTrackerHits(hit);
      ++idx_fill_track;
      if (idx_fill_track >= maxHit)
        break;
    }
  }
}

/**
 * @brief Computes the initial covariance matrix of the track parameters
 *        in Cartesian coordinates starting from helix parameter uncertainties.
 *
 * This method builds an initial estimate of the track parameter covariance
 * in helix parameter space and propagates it to Cartesian coordinates using
 * a transformation that depends on the magnetic field and the initial state
 * (position and momentum).
 *
 * The procedure is as follows:
 *
 * - Initializes a 5x5 diagonal covariance matrix in helix parameter space
 *   (d0, phi, omega, z0, tanLambda), assuming no initial correlations.
 *
 * - d0 (transverse impact parameter):
 *     Uses sigma_d0 if provided, otherwise a default value of 0.05.
 *     Sets the variance as sigma_d0 squared.
 *
 * - phi (azimuthal angle):
 *     Uses sigma_phi if provided, otherwise a default value of 0.1.
 *     Sets the variance as sigma_phi squared.
 *
 * - omega (curvature):
 *     Computes curvature from the initial momentum magnitude in the transverse plane
 *     and the magnetic field Bz, using a conversion constant.
 *     The uncertainty is scaled by omega_factor (default 0.5).
 *     The variance is (omega_factor * omega)^2.
 *
 * - z0 (longitudinal impact parameter):
 *     Uses the absolute value of the initial z position as a scale.
 *     Scales it using z0_factor (default 0.1).
 *     The variance is (z0_factor * |z_init_z|)^2.
 *
 * - tanLambda (dip angle):
 *     Uses sigma_tanLambda if provided, otherwise a default value of 0.1.
 *     Sets the variance as sigma_tanLambda squared.
 *
 * - Assumes no correlations between helix parameters (diagonal covariance only).
 *
 * - Converts the helix covariance matrix to Cartesian coordinates:
 *     The transformation depends on:
 *       - Initial position
 *       - Initial momentum
 *       - Reference point
 *       - Magnetic field Bz
 *
 * - Returns the covariance matrix in Cartesian state representation.
 *
 * @param Bz Magnetic field component along the z axis.
 * @param Charge Particle charge hypothesis (+1 or -1).
 * @param sigma_d0 Optional uncertainty on d0 [cm] (default 0.05).
 * @param sigma_phi Optional uncertainty on phi [rad] (default 0.1).
 * @param omega_factor Scaling factor for curvature uncertainty (default 0.5).
 * @param z0_factor Scaling factor for z0 uncertainty (default 0.1).
 * @param sigma_tanLambda Optional uncertainty on tanLambda (default 0.1).
 *
 * @return Covariance matrix in Cartesian coordinates (6x6).
 */
TMatrixDSym GenfitTrack::ComputeInitialCovarianceMatrix(double Bz, int Charge, std::optional<double> sigma_d0,
                                                        std::optional<double> sigma_phi,
                                                        std::optional<double> omega_factor,
                                                        std::optional<double> z0_factor,
                                                        std::optional<double> sigma_tanLambda) {
  TMatrixDSym C_helix(5);
  C_helix.Zero();

  // d0 (in cm)
  double sd0 = sigma_d0.value_or(0.05);
  C_helix(0, 0) = sd0 * sd0;

  // phi (in rad)
  double sphi = sigma_phi.value_or(0.1);
  C_helix(1, 1) = sphi * sphi;

  // omega (in 1/cm)
  double pt = m_momInit.Perp();
  double omega = std::abs(ConversionUnits::a_lcio * Bz / pt * dd4hep::mm); // in 1/cm
  C_helix(2, 2) = std::pow(omega_factor.value_or(0.5) * omega, 2);

  // z0 (in cm)
  double z0_scale = std::abs(m_posInit.Z());
  C_helix(3, 3) = std::pow(z0_factor.value_or(0.1) * z0_scale, 2);

  // tanLambda
  double stl = sigma_tanLambda.value_or(0.1);
  C_helix(4, 4) = stl * stl;

  // Convert to Cartesian covariance
  TMatrixDSym covState =
      CovarianceMatrixHelixToCartesian(C_helix, m_posInit, m_momInit, m_VP_referencePoint, Charge, Bz);

  return covState;
}

/**
 * @brief Computes initial track parameters (position, momentum, charge)
 *        from tracker hits using geometric fits.
 *
 * The method estimates a first approximation of the particle trajectory
 * starting from a collection of tracker hits. It reconstructs the transverse
 * motion via a circular fit in the XY plane and the longitudinal motion
 * via a linear fit in the R–Z plane.
 *
 * The procedure is as follows:
 *
 * - Extracts 3D hit positions (x, y, z) from the input track.
 *
 * - Projects the hits onto the XY plane and performs a circular fit:
 *     - Uses FastCircleFit to determine the best-fit circle.
 *     - Computes the point of closest approach (PCA) to a reference point.
 *     - Evaluates the tangent direction at the PCA.
 *
 * - Estimates the transverse momentum (pT):
 *     - Uses the relation:
 *           pT = |rho * 0.3 * Bz| / 1000
 *       where rho is the circle radius and Bz is the magnetic field.
 *     - Builds the initial momentum vector in the transverse plane
 *       from the tangent direction.
 *
 * - Reconstructs the longitudinal momentum (pZ):
 *     - Converts each hit into (R, z), where R is the radial distance
 *       from the reference point.
 *     - Performs a linear regression:
 *           z = a * R + b
 *     - Computes:
 *           pZ = a * pT
 *     - Handles degenerate cases (small denominator) by setting pZ = 0.
 *
 * - Determines the z-coordinate at the PCA:
 *     - Uses the intercept b from the linear fit.
 *
 * - Estimates the particle charge:
 *     - Computes the curvature direction from the circle center.
 *     - Computes the Lorentz force direction: p × B.
 *     - Compares the two directions via a dot product:
 *           charge = +1 if aligned, -1 otherwise.
 *
 * - Fills and returns a HelperInitialization structure containing:
 *     - Position  : PCA in 3D (converted to cm)
 *     - Momentum  : (px, py, pz)
 *     - Charge    : +-1
 *
 * @param Bz Magnetic field component along the z-axis.
 *
 * @return HelperInitialization structure with initial track parameters.
 */
GenfitTrack::HelperInitialization GenfitTrack::ComputeInitialParameters(double Bz) {

  Point2D_xy referencePoint_xy(m_VP_referencePoint.X() / dd4hep::mm, m_VP_referencePoint.Y() / dd4hep::mm); // in mm

  auto hits = m_edm4hepTrack.getTrackerHits();
  const size_t N = hits.size();

  if (N < 2) {
    throw std::runtime_error("Not enough hits to compute initial parameters");
  }

  // XY + RZ projection
  std::vector<Point2D_xy> points_xy;
  points_xy.reserve(N);

  double sumR = 0.0, sumZ = 0.0, sumRZ = 0.0, sumR2 = 0.0;

  for (const auto& hit : hits) {

    const auto p = hit.getPosition();

    // XY for circle fit
    points_xy.emplace_back(p.x, p.y);

    // RZ for linear fit
    const double dx = p.x - referencePoint_xy.x;
    const double dy = p.y - referencePoint_xy.y;
    const double R = std::sqrt(dx * dx + dy * dy);
    const double Z = p.z;

    sumR += R;
    sumZ += Z;
    sumRZ += R * Z;
    sumR2 += R * R;
  }

  // Circle fit (XY)
  FastCircleFit circle(points_xy);

  Point2D_xy closestPoint = circle.closestPointTo(referencePoint_xy);
  Point2D_xy tangent_xy = circle.tangentAtPCA(closestPoint, points_xy[1]);

  double rho = circle.rho();
  double init_pT = std::abs(rho * 0.3 * Bz) / 1000.0; // 1000 is used to convert from m to mm
  TVector3 init_mom(tangent_xy.x * init_pT, tangent_xy.y * init_pT, 0.0);

  //  Z fit (linear R-Z)
  double denominator = N * sumR2 - sumR * sumR;

  double a = 0.0;
  double b = 0.0;

  if (std::abs(denominator) > 1e-12 * (N * sumR2)) {
    a = (N * sumRZ - sumR * sumZ) / denominator;
    b = (sumZ - a * sumR) / N;
  } else {
    a = 0.0;
    b = sumZ / N;
  }

  double z_PCA = b;

  init_mom.SetZ(a * init_pT);

  // Charge
  TVector3 pos(closestPoint.x, closestPoint.y, z_PCA);
  TVector3 center(circle.x0(), circle.y0(), 0.0);

  const double rx = center.X() - pos.X();
  const double ry = center.Y() - pos.Y();

  const int sign = (init_mom.X() * ry - init_mom.Y() * rx > 0) ? -1 : 1;

  const int bzSign = (Bz > 0) ? 1 : -1;
  int charge = sign * bzSign;

  // Output
  HelperInitialization helper;
  helper.Position = TVector3(closestPoint.x * dd4hep::mm, closestPoint.y * dd4hep::mm, z_PCA * dd4hep::mm); // in cm

  helper.Momentum = init_mom;
  helper.Charge = charge;

  return helper;
}

/**
 * @brief Builds the genfit::Track from the initialized state and tracker hits.
 *
 * The method performs the following steps:
 *
 * - Sets the signed particle hypothesis (including charge sign handling
 *   for leptons).
 * - Deletes any existing genfit::TrackRep and genfit::Track objects.
 * - Constructs the 6D state vector (x, y, z, px, py, pz) from m_posInit
 *   and m_momInit.
 * - Creates a genfit::RKTrackRep using the signed particle hypothesis.
 * - Instantiates a genfit::Track with the track representation,
 *   state vector, and m_covInit.
 * - Iterates over tracker hits stored in m_edm4hepTrack.
 * - Converts each hit into a Genfit measurement:
 *     - edm4hep::TrackerHitPlane: planar measurement (detID = 0)
 *     - edm4hep::SenseWireHit: wire measurement   (detID = 1)
 * - Wraps each measurement into a genfit::TrackPoint
 *   and inserts it into the genfit::Track.
 * - Terminates execution if an unsupported hit type is found.
 *
 * @param particle_hypotesis PDG code hypothesis (sign adjusted internally).
 * @param debug_lvl Debug verbosity level passed to measurement constructors.
 *
 * @note InitializeTrack() must be called before this method.
 */
void GenfitTrack::CreateGenFitTrack(int particle_hypotesis, int debug_lvl) {

  delete m_genfitTrack;
  delete m_genfitTrackRep;
  m_genfitTrack = nullptr;
  m_genfitTrackRep = nullptr;

  m_signed_particle_hypothesis = particle_hypotesis;
  if (particle_hypotesis == 11 || particle_hypotesis == 13) {
    m_signed_particle_hypothesis = -m_charge_hypothesis * particle_hypotesis;

  } else {
    m_signed_particle_hypothesis = m_charge_hypothesis * particle_hypotesis;
  }

  // Create stateVec
  TVectorD stateVec(6);

  stateVec[0] = m_posInit.X();
  stateVec[1] = m_posInit.Y();
  stateVec[2] = m_posInit.Z();

  stateVec[3] = m_momInit.X();
  stateVec[4] = m_momInit.Y();
  stateVec[5] = m_momInit.Z();

  m_genfitTrackRep = new genfit::RKTrackRep(m_signed_particle_hypothesis);
  m_genfitTrack = new genfit::Track(m_genfitTrackRep, stateVec, m_covInit);

  auto hits_for_genfit = m_edm4hepTrack.getTrackerHits();

  int hit_idx(0);
  int detID(-1);

  // The following code creates genfit::TrackPoints from EDM4hep hits and adds them to the genfit::Track.
  //
  // GENFIT manages ownership of the entire object chain:
  // - A genfit::Track owns its TrackPoints
  // - Deleting a Track automatically deletes all associated TrackPoints
  // - Deleting a TrackPoint also releases the underlying measurement object used to create it
  for (auto hit : hits_for_genfit) {

    if (hit.isA<edm4hep::TrackerHitPlane>()) {
      detID = 0;
      auto planar_hit = hit.as<edm4hep::TrackerHitPlane>();
      GenfitInterface::PlanarMeasurement measurement(planar_hit, detID, ++hit_idx, debug_lvl);
      m_genfitTrack->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), m_genfitTrack));
    } else if (hit.isA<edm4hep::SenseWireHit>()) {
      detID = 1;
      auto wire_hit = hit.as<edm4hep::SenseWireHit>();
      GenfitInterface::WireMeasurement measurement(wire_hit, m_dch_info, m_dc_decoder, detID, ++hit_idx, debug_lvl);
      m_genfitTrack->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), m_genfitTrack));
    } else {
      throw std::runtime_error("InitializeTrack: Unknown hit type encountered - Hit will be skipped.");
    }
  }
}

/**
 * @brief Fit the Genfit track using a chosen track fitting algorithm.
 *
 * This method fits the `GenfitTrack` object using the specified fitting strategy
 * (Deterministic Annealing Filter, standard Kalman, or Kalman with reference track),
 * and updates the associated `edm4hep::Track` track states at key locations:
 * the Interaction Point (IP), the first hit, and the last hit.
 *
 * If DAF is used and `FilterHits` is enabled, track hits are filtered after fitting
 * based on their measurement weights, retaining only accepted hits in `m_fittedHits`.
 *
 * The method also computes helix parameters (d0, z0, phi, omega, tanLambda) for each
 * track state based on the fitted position and momentum, taking into account the
 * assumed charge hypothesis and magnetic field.
 *
 * @param FitterType Fitting strategy to use (https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf):
 *        - "DAF"        : Deterministic Annealing Filter
 *        - "KALMAN"     : Standard Kalman filter
 *        - "KALMAN_REF" : Kalman filter with reference track
 * @param debug_lvl   Verbosity level for debug output
 *                       - 0 = silent,
 *                       - 1 = fitter printout
 *                       - 2 = fitter printout + fit results
 * @param Beta_init   Initial annealing parameter (only for DAF, default = 100)
 * @param Beta_final  Final annealing parameter (only for DAF, default = 0.1)
 * @param Beta_steps  Number of annealing steps (only for DAF, default = 10)
 * @param FilterHits  If true and DAF is used, applies hit filtering after fitting (default = true)
 *
 * @return true if the fit was successful, false otherwise.
 *
 * @note If any exception occurs during fitting or state extrapolation, the function
 *       returns false and does not update the track.
 */
bool GenfitTrack::Fit(std::string FitterType = "DAF", int debug_lvl = 0, std::optional<double> Beta_init = 100.,
                      std::optional<double> Beta_final = 0.1, std::optional<int> Beta_steps = 10,
                      std::optional<bool> FilterHits = true) {

  edm4hep::Track Track_temp = m_edm4hepTrack;
  for (size_t i = 0; i < Track_temp.trackStates_size(); ++i) {

    Track_temp.getTrackStates(i) = edm4hep::TrackState();
  }

  // Initialize the genfit fitter
  std::unique_ptr<genfit::AbsKalmanFitter> genfitFitter = nullptr;

  if (FitterType == "DAF") {

    auto daf = std::make_unique<genfit::DAF>(true, 1e-3, 1e-3);
    daf->setAnnealingScheme(Beta_init.value(), Beta_final.value(), Beta_steps.value());
    daf->setProbCut(1e-5);
    daf->setConvergenceDeltaWeight(1e-2);

    genfitFitter = std::move(daf);

  } else if (FitterType == "KALMAN") {

    genfitFitter = std::make_unique<genfit::KalmanFitter>();

  } else if (FitterType == "KALMAN_REF") {

    genfitFitter = std::make_unique<genfit::KalmanFitterRefTrack>();

  } else {

    throw std::invalid_argument("Unknown fit method: " + FitterType);
  }

  int debug_lvl_fit = debug_lvl;
  if (debug_lvl > 1)
    debug_lvl_fit = 0;
  genfitFitter->setDebugLvl(debug_lvl_fit);

  // Process track
  genfit::Track genfitTrack = *m_genfitTrack;
  genfit::AbsTrackRep* trackRep = genfitTrack.getTrackRep(0);

  try {

    genfitFitter->processTrackWithRep(&genfitTrack, trackRep);

  } catch (const std::exception& e) {

    std::cerr << "Exception during track fitting: " << e.what() << std::endl;
    m_edm4hepTrack.setChi2(-1);
    m_edm4hepTrack.setNdf(-1);

    m_trackWithFit.setChi2(-1);
    m_trackWithFit.setNdf(-1);

    return false;
  }

  if (genfitFitter->isTrackFitted(&genfitTrack, trackRep)) {

    // Hit Filtering based on measurement weights (if FitterType == "DAF" and FilterHits == true)
    if (FilterHits.value() && FitterType == "DAF") {

      int numPoints = genfitTrack.getNumPoints();
      for (int idx = 0; idx < numPoints; idx++) {

        // Retrieve the genfit::TrackPoint for this index and get its associated fitter information.
        // Each TrackPoint can have one or more assigned measurements:
        // - Planar measurements: typically 1 measurement
        // - Drift chamber hits: typically 2 measurements (left and right positions)
        //
        // The loop iterates over all measurements to check their assigned weights and determine
        // whether the hit should be accepted or discarded.
        auto point_track = genfitTrack.getPoint(idx);
        auto fitterInfo = point_track->getFitterInfo(trackRep);
        genfit::KalmanFitterInfo* kfi = static_cast<genfit::KalmanFitterInfo*>(fitterInfo);

        unsigned int nMeas = kfi->getNumMeasurements();
        bool isAccepted = false;

        for (unsigned int j = 0; j < nMeas; j++) {
          genfit::MeasurementOnPlane* mop = kfi->getMeasurementOnPlane(j);
          double weight = mop->getWeight();

          // Weight threshold: 1e-5
          // - Genfit assigns 1e-10 to measurements it discards.
          // - Exception: when left and right hits are very close (within detector resolution),
          //   both may get weight 0.5. In this case, we resolve the ambiguity by selecting the left hit.
          if (weight > 1e-5) {
            isAccepted = true;
            break;
          }
        }

        if (isAccepted) {
          // Retrieve the fitted state from the KalmanFitterInfo
          genfit::StateOnPlane state = kfi->getFittedState();
          genfit::MeasuredStateOnPlane measState = kfi->getFittedState();

          // Extract the 3D position of the fitted state
          TVector3 pos = state.getPos();

          // Extract the measurement plane's orientation vectors
          auto planeMeas = state.getPlane();
          auto U = planeMeas->getU(); // Plane u-direction
          auto V = planeMeas->getV(); // Plane v-direction
          auto O = planeMeas->getO(); // Plane origin

          // Extract covariance matrix in Genfit's 6D (pos+momentum) space
          TMatrixDSym covMatrix_cm_gev = measState.get6DCov();

          // Scale the covariance matrix from cm/GeV to desired mm/GeV
          TMatrixDSym covMatrix = covMatrix_cm_gev;
          for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
              double scale = 1.0;
              bool pos_i = (i < 3);
              bool pos_j = (j < 3);
              if (pos_i && pos_j)
                scale = 100.0; // position-position
              else if (pos_i || pos_j)
                scale = 10.0; // position-momentum
              covMatrix(i, j) *= scale;
            }
          }

          // Extract only the 3x3 position covariance submatrix
          TMatrixDSym covXYZ(3);
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              covXYZ(i, j) = covMatrix(i, j);
            }
          }

          // Convert plane axes to TVector3 for error computation
          TVector3 u(U.X(), U.Y(), U.Z());
          TVector3 v(V.X(), V.Y(), V.Z());

          // Convert plane axes to TVectorD for similarity calculation with covariance
          TVectorD uVec(3);
          uVec[0] = U.X();
          uVec[1] = U.Y();
          uVec[2] = U.Z();
          TVectorD vVec(3);
          vVec[0] = V.X();
          vVec[1] = V.Y();
          vVec[2] = V.Z();

          // Compute variances along plane directions using covariance similarity
          double var_u = covXYZ.Similarity(uVec);
          double var_v = covXYZ.Similarity(vVec);

          // Convert variance to standard deviation (errors)
          double err_u = std::sqrt(var_u);
          double err_v = std::sqrt(var_v);

          // Create a new fitted hit object and set its position
          auto hit3D = m_fittedHits.create();
          hit3D.setPosition(edm4hep::Vector3d(pos.X() / dd4hep::mm, pos.Y() / dd4hep::mm, pos.Z() / dd4hep::mm));

          // Set the 3x3 position covariance matrix in EDM4hep format
          hit3D.setCovMatrix({
              static_cast<float>(covXYZ(0, 0)), // xx
              static_cast<float>(covXYZ(1, 0)), // yx
              static_cast<float>(covXYZ(1, 1)), // yy
              static_cast<float>(covXYZ(2, 0)), // zx
              static_cast<float>(covXYZ(2, 1)), // zy
              static_cast<float>(covXYZ(2, 2))  // zz
          });

          // Set the errors along the local plane axes
          hit3D.setDu(err_u);
          hit3D.setDv(err_v);

          // Set plane orientation in EDM4hep
          edm4hep::Vector2f edm4hepU = {static_cast<float>(u.X()), static_cast<float>(u.Y())};
          edm4hep::Vector2f edm4hepV = {static_cast<float>(v.X()), static_cast<float>(v.Y())};
          hit3D.setU(edm4hepU);
          hit3D.setV(edm4hepV);

          hit3D.setType(1); // Mark as accepted hit
          m_trackWithFit.addToTrackerHits(hit3D);
        } else {
          // Create a placeholder for the rejected hit
          auto hit3D = m_fittedHits.create();
          hit3D.setType(0); // Mark as rejected hit
        }
      }
    }

    genfit::MeasuredStateOnPlane fittedState;

    // trackState First Hit
    fittedState = genfitTrack.getFittedState();
    edm4hep::TrackState trackStateFirstHit = UpdateTrackState(fittedState, edm4hep::TrackState::AtFirstHit);

    // trackState lastHit
    fittedState = genfitTrack.getFittedState(genfitTrack.getNumPoints() - 1);
    edm4hep::TrackState trackStateLastHit = UpdateTrackState(fittedState, edm4hep::TrackState::AtLastHit);

    // Extrapolation to IP
    genfit::TrackPoint* tp = genfitTrack.getPointWithFitterInfo(0);
    auto* fi = static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(trackRep));
    fittedState = fi->getFittedState(true);

    try {
      trackRep->extrapolateToLine(fittedState, TVector3(0, 0, 0), TVector3(0, 0, 1));
    } catch (const std::exception& e) {
      std::cerr << "Exception during extrapolation to IP: " << e.what() << std::endl;

      m_edm4hepTrack.setChi2(-1);
      m_edm4hepTrack.setNdf(-1);
      m_trackWithFit.setChi2(-1);
      m_trackWithFit.setNdf(-1);

      return false;
    }

    edm4hep::TrackState trackStateIP = UpdateTrackState(fittedState, edm4hep::TrackState::AtIP);

    if (debug_lvl == 2) {

      std::cout << "GenfitTrack    DEBUG : TrackState at IP: " << std::endl;
      std::cout << "GenfitTrack    DEBUG :  D0: " << trackStateIP.D0 << " mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  Z0: " << trackStateIP.Z0 << " mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  phi: " << trackStateIP.phi << " rad" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  omega: " << trackStateIP.omega << " 1/mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  tanLambda: " << trackStateIP.tanLambda << std::endl;
      std::cout << "GenfitTrack    DEBUG :  location: " << trackStateIP.location << std::endl;
      std::cout << "GenfitTrack    DEBUG :  reference point: (" << trackStateIP.referencePoint.x << ", "
                << trackStateIP.referencePoint.y << ", " << trackStateIP.referencePoint.z << ") mm" << std::endl;

      std::cout << "\nGenfitTrack   DEBUG : TrackState at First Hit: " << std::endl;
      std::cout << "GenfitTrack    DEBUG :  D0: " << trackStateFirstHit.D0 << " mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  Z0: " << trackStateFirstHit.Z0 << " mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  phi: " << trackStateFirstHit.phi << " rad" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  omega: " << trackStateFirstHit.omega << " 1/mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  tanLambda: " << trackStateFirstHit.tanLambda << std::endl;
      std::cout << "GenfitTrack    DEBUG :  location: " << trackStateFirstHit.location << std::endl;
      std::cout << "GenfitTrack    DEBUG :  reference point: (" << trackStateFirstHit.referencePoint.x << ", "
                << trackStateFirstHit.referencePoint.y << ", " << trackStateFirstHit.referencePoint.z << ") mm"
                << std::endl;

      std::cout << "\nGenfitTrack   DEBUG : TrackState at Last Hit: " << std::endl;
      std::cout << "GenfitTrack    DEBUG :  D0: " << trackStateLastHit.D0 << " mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  Z0: " << trackStateLastHit.Z0 << " mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  phi: " << trackStateLastHit.phi << " rad" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  omega: " << trackStateLastHit.omega << " 1/mm" << std::endl;
      std::cout << "GenfitTrack    DEBUG :  tanLambda: " << trackStateLastHit.tanLambda << std::endl;
      std::cout << "GenfitTrack    DEBUG :  location: " << trackStateLastHit.location << std::endl;
      std::cout << "GenfitTrack    DEBUG :  reference point: (" << trackStateLastHit.referencePoint.x << ", "
                << trackStateLastHit.referencePoint.y << ", " << trackStateLastHit.referencePoint.z << ") mm\n"
                << std::endl;
    }

    m_edm4hepTrack.addToTrackStates(trackStateIP);
    m_edm4hepTrack.addToTrackStates(trackStateFirstHit);
    m_edm4hepTrack.addToTrackStates(trackStateLastHit);

    m_trackWithFit.addToTrackStates(trackStateIP);
    m_trackWithFit.addToTrackStates(trackStateFirstHit);
    m_trackWithFit.addToTrackStates(trackStateLastHit);

    m_edm4hepTrack.setChi2(genfitTrack.getFitStatus()->getChi2());
    m_edm4hepTrack.setNdf(genfitTrack.getFitStatus()->getNdf());

    m_trackWithFit.setChi2(genfitTrack.getFitStatus()->getChi2());
    m_trackWithFit.setNdf(genfitTrack.getFitStatus()->getNdf());

    return true;
  } else {
    m_edm4hepTrack.setChi2(-1);
    m_edm4hepTrack.setNdf(-1);

    m_trackWithFit.setChi2(-1);
    m_trackWithFit.setNdf(-1);

    return false;
  }
}

/**
 * @brief Propagation of the covariance matrix from helix-parameter representation to Cartesian coordinate
 * representation
 *
 * @param C_helix        Covariance Matrix in helix-basis (cm and GeV units)
 * @param Position_cm    Inizial Position in cm
 * @param Momentum_gev   Initial Momentum in gev/c
 * @param RefPoint_cm    Reference position (e.g. IP)
 * @param Charge         Charge hypothesis
 * @param Bz             Bz
 *
 * @return Covariance matrix in cartesian-basis
 */
TMatrixDSym GenfitTrack::CovarianceMatrixHelixToCartesian(const TMatrixDSym& C_helix, // 5x5
                                                          TVector3 Position_cm, TVector3 Momentum_gev,
                                                          TVector3 RefPoint_cm, int Charge, double Bz) {

  double x_PCA = Position_cm.X();
  double y_PCA = Position_cm.Y();

  double px = Momentum_gev.X();
  double py = Momentum_gev.Y();
  double pz = Momentum_gev.Z();

  double pt = Momentum_gev.Perp();
  double phi0 = std::atan2(py, px);

  double d0 = -(RefPoint_cm.X() - x_PCA) * sin(phi0) + (RefPoint_cm.Y() - y_PCA) * cos(phi0);

  double omega = (std::abs(ConversionUnits::a_lcio * Bz / pt)) * dd4hep::mm; // in 1/cm

  if (Bz * Charge < 0)
    omega = -omega;

  // --- Jacobian (6x5) ---
  TMatrixD J(6, 5);
  J.Zero();

  // These definitions are taken from
  // https://flc.desy.de/lcnotes/notes/localfsExplorer_read?currentPath=/afs/desy.de/group/flc/lcnotes/LC-DET-2006-004.pdf

  // The transverse impact parameter d0 is defined as the projection of the
  // displacement vector between the reference point Pr and the PCA P0
  // onto the unit normal vector to the track direction at the PCA.
  //
  // d = P0 - Pr = (x_PCA - Xr, y_PCA - Yr)
  //
  // n_PCA = (-sin(phi0), cos(phi0))
  //
  // Therefore:
  //
  // d0 = n_PCA · d
  //    = (x_PCA - Xr)(-sin(phi0))
  //      + (y_PCA - Yr)(cos(phi0))
  //
  // Since the displacement vector from the reference point to the PCA is
  // parallel to the normal direction:
  //
  // P0 - Pr = d0 * n_PCA
  //
  // the PCA coordinates become:
  //
  // x_PCA = Xr - d0 * sin(phi0)
  // y_PCA = Yr + d0 * cos(phi0)
  J(0, 0) = -sin(phi0);      // dx / dd0
  J(0, 1) = -d0 * cos(phi0); // dx / dphi0
  J(0, 2) = 0.0;             // dx / domega
  J(0, 3) = 0.0;             // dx / dz0
  J(0, 4) = 0.0;             // dx / dtanLambda

  J(1, 0) = cos(phi0);       // dy / dd0
  J(1, 1) = -d0 * sin(phi0); // dy / dphi0
  J(1, 2) = 0.0;             // dy / domega
  J(1, 3) = 0.0;             // dy / dz0
  J(1, 4) = 0.0;             // dy / dtanLambda

  // z = z_PCA = P^0_z = z0 + P^r_z
  J(2, 0) = 0.0; // dz / dd0
  J(2, 1) = 0.0; // dz / dphi0
  J(2, 2) = 0.0; // dz / domega
  J(2, 3) = 1.0; // dz / dz0
  J(2, 4) = 0.0; // dz / dtanLambda

  // px = pt * cos(phi0) = a * Bz / omega * cos(phi0)
  J(3, 0) = 0.0;         // dpx / dd0
  J(3, 1) = -py;         // dpx / dphi0
  J(3, 2) = -px / omega; // dpx / domega
  J(3, 3) = 0.0;         // dpx / dz0
  J(3, 4) = 0.0;         // dpx / dtanLambda

  // py = pt * sin(phi0) = a * Bz / omega * sin(phi0)
  J(4, 0) = 0.0;         // dpy / dd0
  J(4, 1) = px;          // dpy / dphi0
  J(4, 2) = -py / omega; // dpy / domega
  J(4, 3) = 0.0;         // dpy / dz0
  J(4, 4) = 0.0;         // dpy / dtanLambda

  // pz = pt * tanLambda = a * Bz / omega * tanLambda = p * cos(theta) = p * cos(cot^-1(tanLambda))
  J(5, 0) = 0.0;         // dpz / dd0
  J(5, 1) = 0.0;         // dpz / dphi0
  J(5, 2) = -pz / omega; // dpz / domega
  J(5, 3) = 0.0;         // dpz / dz0
  J(5, 4) = pt;          // dpz / dtanLambda

  // --- Compute C_cart = J * C_helix * J^T ---
  TMatrixD Jt(TMatrixD::kTransposed, J);
  TMatrixD tmp = J * C_helix * Jt;

  TMatrixDSym C_cart(6);
  C_cart.Zero();
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j <= i; ++j) {
      C_cart(i, j) = tmp(i, j);
    }
  }
  for (int i = 0; i < 6; ++i) {
    for (int j = i + 1; j < 6; ++j) {
      C_cart(i, j) = C_cart(j, i);
    }
  }

  return C_cart;
}

/**
 * @brief Propagation of the covariance matrix from cartesian-parameter representation to Helix coordinate
 * representation
 *
 * @param C_cartesian    Covariance Matrix in cartesian-basis (cm and GeV units)
 * @param Position_cm    Inizial Position in cm
 * @param Momentum_gev   Initial Momentum in gev/c
 * @param RefPoint_cm    Reference position (e.g. IP)
 * @param Charge         Charge hypothesis
 * @param Bz             Bz
 *
 * @return Covariance matrix in helix-basis (mm and GeV units)
 */
TMatrixDSym GenfitTrack::CovarianceMatrixCartesianToHelix(const TMatrixDSym& C_cartesian, // 6x6
                                                          TVector3 Position_cm, TVector3 Momentum_gev,
                                                          TVector3 RefPoint_cm, int Charge, double Bz) {

  double px = Momentum_gev.X();
  double py = Momentum_gev.Y();
  double pz = Momentum_gev.Z();

  double x_PCA_mm = Position_cm.X() / dd4hep::mm;
  double y_PCA_mm = Position_cm.Y() / dd4hep::mm;
  double RefX_mm = RefPoint_cm.X() / dd4hep::mm;
  double RefY_mm = RefPoint_cm.Y() / dd4hep::mm;

  double pt = Momentum_gev.Perp();

  double phi0 = std::atan2(py, px);
  double tanLambda = pz / pt;
  double omega = (std::abs(ConversionUnits::a_lcio * Bz / pt)); // in 1/mm

  if (Bz * Charge < 0)
    omega = -omega;

  // --- Jacobian (5x6) ---
  TMatrixD J(5, 6);
  J.Zero();

  // These definitions are taken from
  // https://flc.desy.de/lcnotes/notes/localfsExplorer_read?currentPath=/afs/desy.de/group/flc/lcnotes/LC-DET-2006-004.pdf

  // d0 = -(Xr - x)*sin(phi0) + (Yr - y)*cos(phi0)
  J(0, 0) = sin(phi0);  // dd0 / dx
  J(0, 1) = -cos(phi0); // dd0 / dy
  J(0, 2) = 0.0;        // dd0 / dz

  // chain rule via phi0 = atan2(py, px)
  double dd0_dphi = -(RefX_mm - x_PCA_mm) * cos(phi0) - (RefY_mm - y_PCA_mm) * sin(phi0);

  J(0, 3) = dd0_dphi * (-py / (pt * pt)); // dd0 / dpx
  J(0, 4) = dd0_dphi * (px / (pt * pt));  // dd0 / dpy
  J(0, 5) = 0.0;                          // dd0 / dpz

  // phi0 = atan2(py, px);
  J(1, 0) = 0.0;                   // dphi0 / dx
  J(1, 1) = 0.0;                   // dphi0 / dy
  J(1, 2) = 0.0;                   // dphi0 / dz
  J(1, 3) = -py / std::pow(pt, 2); // dphi0 / dpx
  J(1, 4) = px / std::pow(pt, 2);  // dphi0 / dpy
  J(1, 5) = 0.0;                   // dphi0 / dpz

  // omega = q * a * Bz / pt, with pt = sqrt(px^2 + py^2)
  double pt2 = px * px + py * py;
  J(2, 0) = 0.0;               // domega / dx_PCA
  J(2, 1) = 0.0;               // domega / dy_PCA
  J(2, 2) = 0.0;               // domega / dz_PCA
  J(2, 3) = -omega * px / pt2; // domega / dpx
  J(2, 4) = -omega * py / pt2; // domega / dpy
  J(2, 5) = 0.0;               // domega / dpz

  // z0 = P^0_z - P^r_z
  J(3, 0) = 0.0; // dz0 / dx_PCA
  J(3, 1) = 0.0; // dz0 / dy_PCA
  J(3, 2) = 1.0; // dz0 / dz_PCA
  J(3, 3) = 0.0; // dz0 / dpx
  J(3, 4) = 0.0; // dz0 / dpy
  J(3, 5) = 0.0; // dz0 / dpz

  // tanLambda = pz / pt
  J(4, 0) = 0.0;                               // dtanLambda / dx_PCA
  J(4, 1) = 0.0;                               // dtanLambda / dy_PCA
  J(4, 2) = 0.0;                               // dtanLambda / dz_PCA
  J(4, 3) = -tanLambda * px / std::pow(pt, 2); // dtanLambda / dpx
  J(4, 4) = -tanLambda * py / std::pow(pt, 2); // dtanLambda / dpy
  J(4, 5) = 1.0 / pt;                          // dtanLambda / dpz

  // --- Compute C_cart = J * C_helix * J^T ---
  TMatrixD Jt(TMatrixD::kTransposed, J);
  TMatrixD tmp = J * C_cartesian * Jt;

  TMatrixDSym C_helix(5);
  C_helix.Zero();
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      C_helix(i, j) = tmp(i, j);
    }
  }
  for (int i = 0; i < 5; ++i) {
    for (int j = i + 1; j < 5; ++j) {
      C_helix(i, j) = C_helix(j, i);
    }
  }

  return C_helix;
}

/**
 * @brief Update an edm4hep::TrackState using the fitted Genfit state parameters.
 *
 * This method extracts position, momentum, and covariance information from a
 * `genfit::MeasuredStateOnPlane` object and converts them into the corresponding
 * `edm4hep::TrackState` representation, including helix parameters.
 *
 * The procedure includes:
 *   - Extraction of position, momentum, and covariance from the Genfit state
 *   - Optional inversion of the momentum direction for states at the Interaction Point (IP)
 *   - Retrieval of the local magnetic field (Bz) at the track position
 *   - Computation of the Point of Closest Approach (PCA) to a reference point
 *   - Derivation of helix parameters:
 *       - d0          : transverse impact parameter (mm)
 *       - z0          : longitudinal impact parameter (mm)
 *       - phi         : azimuthal angle of momentum (rad)
 *       - omega       : curvature (1/mm), signed according to charge hypothesis
 *       - tanLambda   : dip angle (pz / pt)
 *   - Assignment of the reference point and track state location
 *
 * The PCA and associated parameters are computed assuming a helical trajectory
 * in a magnetic field aligned along the z-axis.
 *
 * @param Edm4hepTrackState Output track state to be updated
 * @param MeasuredState Genfit measured state containing fitted position and momentum
 * @param location Location identifier for the track state (e.g. AtIP, AtFirstHit, AtLastHit)
 *
 * @note The momentum is reversed for states at the Interaction Point (AtIP) to ensure
 *       a consistent track parameter convention.
 * @note The magnetic field is retrieved at the current position and converted to Tesla.
 * @note The curvature sign (omega) depends on the assumed particle charge hypothesis.
 */
edm4hep::TrackState GenfitTrack::UpdateTrackState(genfit::MeasuredStateOnPlane MeasuredState, int location) {

  edm4hep::TrackState Edm4hepTrackState;

  TVector3 gen_position, gen_momentum;
  TMatrixDSym covariancePosMom(6);

  MeasuredState.getPosMomCov(gen_position, gen_momentum, covariancePosMom);

  double x_ref = gen_position.X(); // cm
  double y_ref = gen_position.Y(); // cm
  double z_ref = gen_position.Z(); // cm
  double pz = gen_momentum.Z();    // gev
  double pt = gen_momentum.Perp(); // gev

  double Bz = m_fieldMap->getBz(gen_position) / (dd4hep::tesla / dd4hep::kilogauss); // From kilogauss to Tesla
  auto infoComputeD0Z0_firstHit =
      PCAInfo(gen_position, gen_momentum, m_charge_hypothesis, m_VP_referencePoint, Bz); // in cm

  double d0 = ((-(m_VP_referencePoint.X() - infoComputeD0Z0_firstHit.PCA.X())) * sin(infoComputeD0Z0_firstHit.Phi0) +
               (m_VP_referencePoint.Y() - infoComputeD0Z0_firstHit.PCA.Y()) * cos(infoComputeD0Z0_firstHit.Phi0)) /
              dd4hep::mm;                                                                // mm
  double z0 = (infoComputeD0Z0_firstHit.PCA.Z() - m_VP_referencePoint.Z()) / dd4hep::mm; // mm
  double phi = gen_momentum.Phi();                                                       // rad

  double tanLambda = pz / pt;
  double omega = std::abs(ConversionUnits::a_lcio * Bz / pt);
  if (m_charge_hypothesis < 0)
    omega = -omega;

  Edm4hepTrackState.D0 = d0;
  Edm4hepTrackState.Z0 = z0;
  Edm4hepTrackState.phi = phi;
  Edm4hepTrackState.omega = omega;
  Edm4hepTrackState.tanLambda = tanLambda;
  Edm4hepTrackState.time = 0.;

  Edm4hepTrackState.referencePoint = edm4hep::Vector3f(x_ref / dd4hep::mm, y_ref / dd4hep::mm, z_ref / dd4hep::mm);
  Edm4hepTrackState.location = location;

  TMatrixDSym CovHelix = CovarianceMatrixCartesianToHelix(covariancePosMom, gen_position, gen_momentum,
                                                          m_VP_referencePoint, m_charge_hypothesis, Bz);

  // Conversion from TMatrixDSym(5x5) to lower-triangular packed format used in edm4hep::TrackState
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      Edm4hepTrackState.setCovMatrix(static_cast<float>(CovHelix(i, j)), static_cast<edm4hep::TrackParams>(i),
                                     static_cast<edm4hep::TrackParams>(j));
    }
  }

  return Edm4hepTrackState;
}

/**
 * @brief Computes the point of closest approach (PCA) and related parameters
 *        of a charged particle's trajectory in a uniform magnetic field.
 *
 * This function approximates the helical motion of a charged particle in a
 * uniform magnetic field oriented along the z-axis (Bz) and calculates:
 *  - The point on the helix closest to a reference point (`refPoint`) in 3D space.
 *  - The azimuthal angle (`Phi0`) of the trajectory at the PCA.
 *
 * The procedure is as follows:
 * 1. Compute the transverse momentum `pt` and check it is non-zero.
 * 2. Calculate the radius `R` of the helix in the transverse plane using
 *    R = pt / (0.3 * |q| * Bz) * 100 (pt in GeV/c, Bz in Tesla, R in cm).
 * 3. Determine the center of the circular projection of the helix in the XY plane.
 * 4. Find the unit vector pointing from the circle center to the reference point.
 * 5. Compute the closest point on the circle to the reference point in XY.
 * 6. Compute the tangent vector at the closest point in the XY plane.
 * 7. Calculate the azimuthal angle `Phi0` of the tangent.
 * 8. Compute the z-coordinate of the PCA along the helix using the straight-line approximation.
 *
 * @param position  Initial 3D position of the particle (in cm).
 * @param momentum  3D momentum vector of the particle (in GeV/c).
 * @param charge    Particle charge (in e units).
 * @param refPoint  Reference point to which the PCA is calculated (in cm).
 * @param Bz        Magnetic field along the z-axis (Tesla).
 *
 * @return A `PCAInfoHelper` structure containing:
 *         - `PCA`   : The 3D position of the closest approach point (in cm).
 *         - `Phi0`  : Azimuthal angle of the trajectory at the PCA (in rad).
 *
 * @note Assumes a uniform magnetic field along z and a simple helical trajectory.
 */
GenfitTrack::PCAInfoHelper GenfitTrack::PCAInfo(TVector3 position, TVector3 momentum, int charge, TVector3 refPoint,
                                                double Bz) {

  double px = momentum.X();
  double py = momentum.Y();
  double pt = momentum.Perp();

  double R = pt / (0.3 * std::abs(charge) * Bz) * 100; // in cm (assuming pt in GeV/c, Bz in Tesla)

  double tx = px / pt;
  double ty = py / pt;

  double nx = charge * (ty);
  double ny = charge * (-tx);

  double xc = position[0] + R * nx;
  double yc = position[1] + R * ny;
  TVector3 center(xc, yc, 0.0);

  TVector3 v = refPoint - center;

  double vxy = std::sqrt(v.X() * v.X() + v.Y() * v.Y());

  TVector3 u(v.X() / vxy, v.Y() / vxy, 0.0);
  TVector3 closestPoint = center + R * u;
  TVector3 r = closestPoint - center;

  int sign = (charge > 0) ? 1 : -1;
  TVector3 tangent(-sign * r.Y(), sign * r.X(), 0.0);

  tangent = tangent.Unit();

  double angleToX = std::atan2(tangent.Y(), tangent.X());

  double pR = momentum.Perp();
  double pZ = momentum.Z();
  double R0 = position.Perp();
  double Z0 = position.Z();

  double tPCA = -(R0 * pR + Z0 * pZ) / (pR * pR + pZ * pZ);
  double ZPCA = Z0 + pZ * tPCA;

  PCAInfoHelper helper;

  closestPoint.SetZ(ZPCA);
  helper.PCA = closestPoint;

  helper.Phi0 = angleToX;

  return helper;
}

} // namespace GenfitInterface