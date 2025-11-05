#include "utils.hpp"

dd4hep::rec::LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag) {

  dd4hep::rec::LayeredCalorimeterData* theExtension = 0;
  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  if (detector.detectors().empty()) {
    throw std::runtime_error("Detector is not initialized: no DetElements found.");
  }

  const std::vector<dd4hep::DetElement>& theDetectors =
      dd4hep::DetectorSelector(detector).detectors(includeFlag, excludeFlag);

  if (theDetectors.empty()) {
    // this is a standalone function, so we cannot use any message service
    std::cout << "k4RecTracker/Tracking/src/utils.cpp: No detectors found for the given selection:" << std::endl;
    std::cout << "  includeFlag: " << dd4hep::DetType(includeFlag) << " excludeFlag: " << dd4hep::DetType(excludeFlag)
              << std::endl;

    return nullptr;
  }

  int debug_lvl = 0;
  if (debug_lvl > 0) {
    std::cout << " getExtension :  includeFlag: " << dd4hep::DetType(includeFlag)
              << " excludeFlag: " << dd4hep::DetType(excludeFlag) << "  found : " << theDetectors.size()
              << "  - first det: " << theDetectors.at(0).name() << std::endl;
  }

  if (theDetectors.size() != 1) {

    std::stringstream es;
    es << " getExtension: selection is not unique - includeFlag: " << dd4hep::DetType(includeFlag)
       << " excludeFlag: " << dd4hep::DetType(excludeFlag) << " --- found detectors : ";
    for (unsigned i = 0, N = theDetectors.size(); i < N; ++i) {
      es << theDetectors.at(i).name() << ", ";
    }

    throw std::runtime_error(es.str());
  }

  theExtension = theDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();

  return theExtension;
}

edm4hep::TrackState getExtrapolationAtCalorimeter(const pandora::CartesianVector& ecalProjection,
                                                  const HelixClass_double& helixAtLastHit, double Bz) {

  edm4hep::TrackState trackState_AtCalorimeter = edm4hep::TrackState{};

  double posAtCalorimeter[] = {ecalProjection.GetX(), ecalProjection.GetY(), ecalProjection.GetZ()};
  int debug_lvl = 0;
  if (debug_lvl > 0) {
    std::cout << "Projection at calo: x, y, z, r = " << posAtCalorimeter[0] << " " << posAtCalorimeter[1] << " "
              << posAtCalorimeter[2] << " "
              << sqrt(posAtCalorimeter[0] * posAtCalorimeter[0] + posAtCalorimeter[1] * posAtCalorimeter[1])
              << std::endl;
  }

  // get extrapolated momentum from the helix with ref point at last hit
  double momAtCalorimeter[] = {0., 0., 0.};
  helixAtLastHit.getExtrapolatedMomentum(posAtCalorimeter, momAtCalorimeter);

  // produce new helix at calorimeter position
  auto helixAtCalorimeter = HelixClass_double();
  helixAtCalorimeter.Initialize_VP(posAtCalorimeter, momAtCalorimeter, helixAtLastHit.getCharge(), Bz);

  // fill the TrackState parameters
  trackState_AtCalorimeter.location = edm4hep::TrackState::AtCalorimeter;
  trackState_AtCalorimeter.D0 = helixAtCalorimeter.getD0();
  trackState_AtCalorimeter.phi = std::atan2(momAtCalorimeter[1], momAtCalorimeter[0]);
  trackState_AtCalorimeter.omega = helixAtCalorimeter.getOmega();
  trackState_AtCalorimeter.Z0 = helixAtCalorimeter.getZ0();
  trackState_AtCalorimeter.tanLambda = helixAtCalorimeter.getTanLambda();
  trackState_AtCalorimeter.referencePoint =
      edm4hep::Vector3f(posAtCalorimeter[0], posAtCalorimeter[1], posAtCalorimeter[2]);
  return trackState_AtCalorimeter;
}

void FillTrackWithCalorimeterExtrapolation(extension::MutableTrack& edm4hep_track, double m_Bz, int charge, double a,
                                           double m_eCalBarrelInnerR, double m_eCalBarrelMaxZ,
                                           double m_eCalEndCapInnerR, double m_eCalEndCapOuterR,
                                           double m_eCalEndCapInnerZ) {

  auto trackStateLastHit = edm4hep_track.getTrackStates()[2];
  double omega_lastHit = trackStateLastHit.omega;
  double pt_lasthit = a * m_Bz / abs(omega_lastHit);
  double phi_lasthit = trackStateLastHit.phi;
  double pz_lasthit = trackStateLastHit.tanLambda * pt_lasthit;
  double px_lasthit = pt_lasthit * std::cos(phi_lasthit);
  double py_lasthit = pt_lasthit * std::sin(phi_lasthit);
  auto ref_lastHit = trackStateLastHit.referencePoint;

  // produce new helix at last hit position
  double posAtLastHit[] = {ref_lastHit[0], ref_lastHit[1], ref_lastHit[2]};
  double momAtLastHit[] = {px_lasthit, py_lasthit, pz_lasthit};
  auto helixAtLastHit = HelixClass_double();
  helixAtLastHit.Initialize_VP(posAtLastHit, momAtLastHit, charge, m_Bz);

  // Propagation to Endcap
  if (m_eCalBarrelInnerR > 0. || m_eCalEndCapInnerR > 0.) {

    pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
    pandora::CartesianVector secondBestECalProjection(0.f, 0.f, 0.f);
    float minGenericTime(std::numeric_limits<float>::max());

    // create helix to project
    // rather than using parameters at production, better to use those from
    // last hit
    pandora::CartesianVector pos_lasthit(posAtLastHit[0], posAtLastHit[1], posAtLastHit[2]);
    pandora::CartesianVector mom_lasthit(momAtLastHit[0], momAtLastHit[1], momAtLastHit[2]);

    const pandora::Helix helix(pos_lasthit, mom_lasthit, charge, m_Bz);
    const pandora::CartesianVector& referencePoint(helix.GetReferencePoint());
    const int signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);

    // First project to endcap
    pandora::CartesianVector endCapProjection(0.f, 0.f, 0.f);
    if (m_eCalEndCapInnerR > 0) {
      float genericTime(std::numeric_limits<float>::max());
      const pandora::StatusCode statusCode(helix.GetPointInZ(static_cast<float>(signPz) * m_eCalEndCapInnerZ,
                                                             referencePoint, endCapProjection, genericTime));
      float x = endCapProjection.GetX();
      float y = endCapProjection.GetY();
      float r = std::sqrt(x * x + y * y);
      if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime) && (r >= m_eCalEndCapInnerR) &&
          (r <= m_eCalEndCapOuterR)) {
        minGenericTime = genericTime;
        bestECalProjection = endCapProjection;
      }
    }

    // Then project to barrel surface(s), and keep projection
    // if extrapolation is within the z acceptance of the detector
    pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
    if (m_eCalBarrelInnerR > 0) {

      float genericTime(std::numeric_limits<float>::max());
      const pandora::StatusCode statusCode(
          helix.GetPointOnCircle(m_eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));

      if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (std::fabs(barrelProjection.GetZ()) <= m_eCalBarrelMaxZ)) {
        if (genericTime < minGenericTime) {
          minGenericTime = genericTime;
          secondBestECalProjection = bestECalProjection;
          bestECalProjection = barrelProjection;
        } else {
          secondBestECalProjection = barrelProjection;
        }
      }
    }

    // store extrapolation to calo
    // by default, store extrapolation with lower arrival time
    // get extrapolated position
    edm4hep::TrackState trackState_AtCalorimeter =
        getExtrapolationAtCalorimeter(bestECalProjection, helixAtLastHit, m_Bz);

    // attach the TrackState to the track
    edm4hep_track.addToTrackStates(trackState_AtCalorimeter);
  }
}

torch::Tensor find_condpoints(const torch::Tensor& betas, const torch::Tensor& unassigned, float tbeta) {
  torch::Tensor mask_unassigned = torch::zeros_like(betas, torch::dtype(torch::kBool));
  mask_unassigned.index_put_({unassigned}, true);

  torch::Tensor select_condpoints = (betas > tbeta) & mask_unassigned;
  torch::Tensor indices_condpoints = torch::nonzero(select_condpoints).view({-1});

  if (indices_condpoints.numel() == 0) {
    return torch::empty({0}, torch::dtype(torch::kLong));
  }

  torch::Tensor betas_condpoints = betas.index_select(0, indices_condpoints);
  torch::Tensor sorted_indices = std::get<1>(betas_condpoints.sort(/*dim=*/0, /*descending=*/true));

  return indices_condpoints.index_select(0, sorted_indices);
}

torch::Tensor get_clustering(const std::vector<float>& output_vector, int num_rows, float tbeta, float td) {
  // Create tensor from input vector and reshape to (num_rows, 4)
  torch::Tensor output_model_tensor = torch::from_blob(const_cast<float*>(output_vector.data()), {num_rows, 4},
                                                       torch::dtype(torch::kFloat32))
                                          .clone(); // clone to make it contiguous and owned

  torch::Tensor X = output_model_tensor.slice(1, 0, 3);   // columns 0,1,2
  torch::Tensor betas = output_model_tensor.select(1, 3); // column 3

  int64_t n_points = betas.size(0);
  torch::Tensor clustering = torch::zeros({n_points}, torch::dtype(torch::kLong));
  torch::Tensor unassigned = torch::arange(n_points, torch::dtype(torch::kLong));

  int64_t index_assignation = 1;
  torch::Tensor indices_condpoints = find_condpoints(betas, unassigned, tbeta);

  while (indices_condpoints.numel() > 0 && unassigned.numel() > 0) {
    int64_t index_condpoint = indices_condpoints[0].item<int64_t>();
    torch::Tensor coord_cond = X[index_condpoint];

    torch::Tensor dists = torch::norm(X.index_select(0, unassigned) - coord_cond, 2, /*dim=*/1);
    torch::Tensor mask_distance = dists < td;

    if (mask_distance.sum().item<int64_t>() == 0) {
      indices_condpoints = indices_condpoints.slice(0, 1); // remove first
      continue;
    }

    torch::Tensor assigned_points = torch::masked_select(unassigned, mask_distance);
    clustering.index_put_({assigned_points}, index_assignation);

    torch::Tensor mask_keep = ~mask_distance;
    unassigned = torch::masked_select(unassigned, mask_keep);

    indices_condpoints = find_condpoints(betas, unassigned, tbeta);
    index_assignation += 1;
  }

  return clustering;
}

int getHypotesisCharge(int pdg) {
  switch (pdg) {
  case 11:
    return -1; // electron
  case -11:
    return 1; // positron
  case 13:
    return -1; // muon
  case -13:
    return 1; // anti-muon
  case 211:
    return 1; // pion+
  case -211:
    return -1; // pion-
  case 321:
    return 1; // kaon+
  case -321:
    return -1; // kaon-
  case 2212:
    return 1; // proton
  case -2212:
    return -1; // anti-proton
  default:
    return 0; // unknown particle
  }
}

TMatrixDSym computeTrackStateCovMatrix(TVectorD stateTrack, TVectorD params, TVector3 referencePoint, double timeError,
                                       TMatrixDSym statecovMatrix) {

  // NB: we should give to the matrix the following units:
  //      - positions: mm
  //      - momentum: GeV
  //      - time: ns

  TVector3 pos_0(stateTrack[0], stateTrack[1], stateTrack[2]);
  TVector3 mom(stateTrack[3], stateTrack[4], stateTrack[5]);

  double pt = mom.Pt();
  double omega = params[0];
  double phi0 = params[1];
  double tanLambda = params[5];

  TMatrixD J(5, 6);

  double dOmega_dX = 0;
  double dOmega_dY = 0;
  double dOmega_dZ = 0;
  double dOmega_dPx = -omega * pt / pow(pt, 3) * mom.X();
  double dOmega_dPy = -omega * pt / pow(pt, 3) * mom.Y();
  double dOmega_dPz = 0;

  J(0, 0) = dOmega_dX;
  J(0, 1) = dOmega_dY;
  J(0, 2) = dOmega_dZ;
  J(0, 3) = dOmega_dPx;
  J(0, 4) = dOmega_dPy;
  J(0, 5) = dOmega_dPz;

  double dPhi0_dX = 0;
  double dPhi0_dY = 0;
  double dPhi0_dZ = 0;
  double dPhi0_dPx = -mom.Y() / pow(pt, 2);
  double dPhi0_dPy = mom.X() / pow(pt, 2);
  double dPhi0_dPz = 0;

  J(1, 0) = dPhi0_dX;
  J(1, 1) = dPhi0_dY;
  J(1, 2) = dPhi0_dZ;
  J(1, 3) = dPhi0_dPx;
  J(1, 4) = dPhi0_dPy;
  J(1, 5) = dPhi0_dPz;

  double dD0_dX = sin(phi0);
  double dD0_dY = -cos(phi0);
  double dD0_dZ = 0;
  double dD0_dPx = -(referencePoint.X() - pos_0.X()) * sin(phi0) * dPhi0_dPx -
                   (referencePoint.Y() - pos_0.Y()) * cos(phi0) * dPhi0_dPx;
  double dD0_dPy = -(referencePoint.X() - pos_0.X()) * sin(phi0) * dPhi0_dPy -
                   (referencePoint.Y() - pos_0.Y()) * cos(phi0) * dPhi0_dPy;
  double dD0_dPz = 0;

  J(2, 0) = dD0_dX;
  J(2, 1) = dD0_dY;
  J(2, 2) = dD0_dZ;
  J(2, 3) = dD0_dPx;
  J(2, 4) = dD0_dPy;
  J(2, 5) = dD0_dPz;

  double dZ0_dX = 0;
  double dZ0_dY = 0;
  double dZ0_dZ = 1;
  double dZ0_dPx = 0;
  double dZ0_dPy = 0;
  double dZ0_dPz = 0;

  J(3, 0) = dZ0_dX;
  J(3, 1) = dZ0_dY;
  J(3, 2) = dZ0_dZ;
  J(3, 3) = dZ0_dPx;
  J(3, 4) = dZ0_dPy;
  J(3, 5) = dZ0_dPz;

  double dTanLambda_dX = 0;
  double dTanLambda_dY = 0;
  double dTanLambda_dZ = 0;
  double dTanLambda_dPx = -tanLambda * pt / pow(pt, 3) * mom.X();
  double dTanLambda_dPy = -tanLambda * pt / pow(pt, 3) * mom.Y();
  double dTanLambda_dPz = 1 / pt;

  J(4, 0) = dTanLambda_dX;
  J(4, 1) = dTanLambda_dY;
  J(4, 2) = dTanLambda_dZ;
  J(4, 3) = dTanLambda_dPx;
  J(4, 4) = dTanLambda_dPy;
  J(4, 5) = dTanLambda_dPz;

  TMatrixDSym covarianceTrack_temp(5);
  covarianceTrack_temp = statecovMatrix.Similarity(J);

  TMatrixDSym covarianceTrackState(6);
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j <= i; ++j) {
      double val = covarianceTrack_temp(i, j);
      covarianceTrackState(i, j) = val;
      covarianceTrackState(j, i) = val;
    }
  }

  covarianceTrackState(5, 5) = timeError * timeError;

  return covarianceTrackState;
}
