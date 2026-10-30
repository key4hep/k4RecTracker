/*
 * Copyright (c) 2020-2026 Key4hep-Project.
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

#include "BremRecoveryUtils.h"

namespace {
// pT [GeV] = gevPerTeslaMm * Bz [T] * R [mm] = gevPerTeslaMm * Bz [T] / |omega| [1/mm]
// (same value as ConversionUnits::a_lcio in Tracking/include/utils.h)
constexpr double gevPerTeslaMm = 2.998e-4;
} // namespace

edm4hep::TrackState getExtrapolationAtCalorimeter(const pandora::CartesianVector& ecalProjection,
                                                  const HelixClass_double& helixAtLastHit, double Bz) {

  edm4hep::TrackState trackState_AtCalorimeter = edm4hep::TrackState{};

  double posAtCalorimeter[] = {ecalProjection.GetX(), ecalProjection.GetY(), ecalProjection.GetZ()};

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

edm4hep::TrackState helixAtLastTrack(const edm4hep::Track& track, double eCalBarrelInnerR, double eCalBarrelMaxZ,
                                     double eCalEndCapInnerR, double eCalEndCapOuterR, double eCalEndCapInnerZ,
                                     double Bz) {

  edm4hep::TrackState trackState_AtCalorimeter;
  for (const auto& state : track.getTrackStates()) { // start of state loop

    /* track needs to be extrapolated from its last hit */
    if (state.location == edm4hep::TrackState::AtLastHit) { // start of state.location if

      double phi = state.phi;
      double omega = state.omega;
      double tanLambda = state.tanLambda;
      double pt = gevPerTeslaMm * Bz / std::abs(omega);

      double px = pt * std::cos(phi);
      double py = pt * std::sin(phi);
      double pz = pt * tanLambda;
      auto ref = state.referencePoint;
      double posAtLastHit[] = {ref.x, ref.y, ref.z};
      double momAtLastHit[] = {px, py, pz};
      int charge = (omega > 0.f) ? 1 : -1;
      // produce new helix at last hit position
      auto helixAtLastHit = HelixClass_double();
      helixAtLastHit.Initialize_VP(posAtLastHit, momAtLastHit, charge, Bz);
      // TrackState at Calorimeter
      if (eCalBarrelInnerR > 0. || eCalEndCapInnerR > 0.) { // start of eCalBarrel||eCalEndCap if

        pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
        pandora::CartesianVector secondBestECalProjection(0.f, 0.f, 0.f);
        float minGenericTime(std::numeric_limits<float>::max());

        // create helix to project
        // rather than using parameters at production, better to use those from
        // last hit
        pandora::CartesianVector pos(ref.x, ref.y, ref.z);
        pandora::CartesianVector mom(px, py, pz);
        const pandora::Helix helix(pos, mom, charge, Bz);
        const pandora::CartesianVector& referencePoint(helix.GetReferencePoint());
        const int signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);

        // First project to endcap
        pandora::CartesianVector endCapProjection(0.f, 0.f, 0.f);
        // bool hasEndCapProjection(false);
        if (eCalEndCapInnerR > 0) { // start of 1st eCalEndCapInnerR if

          float genericTime(std::numeric_limits<float>::max());
          const pandora::StatusCode statusCode(helix.GetPointInZ(static_cast<float>(signPz) * eCalEndCapInnerZ,
                                                                 referencePoint, endCapProjection, genericTime));
          float x = endCapProjection.GetX();
          float y = endCapProjection.GetY();
          float r = std::sqrt(x * x + y * y);
          if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime) &&
              (r >= eCalEndCapInnerR) && (r <= eCalEndCapOuterR)) { // start of eCalEndCapInnerR if

            minGenericTime = genericTime;
            bestECalProjection = endCapProjection;
            // hasEndCapProjection = true;
          } // end of eCalEndCapInnerR if
        } // end of 1st eCalEndCapInnerR if
        // Then project to barrel surface(s), and keep projection
        // if extrapolation is within the z acceptance of the detector
        pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
        // bool hasBarrelProjection = false;
        if (eCalBarrelInnerR > 0) { // start of eCalBarrelInnerR if

          float genericTime(std::numeric_limits<float>::max());
          const pandora::StatusCode statusCode(
              helix.GetPointOnCircle(eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));
          if ((pandora::STATUS_CODE_SUCCESS == statusCode) &&
              (std::fabs(barrelProjection.GetZ()) <= eCalBarrelMaxZ)) { // start of barrel projection if

            // hasBarrelProjection = true;
            if (genericTime < minGenericTime) { // start of generic time if

              minGenericTime = genericTime;
              secondBestECalProjection = bestECalProjection;
              bestECalProjection = barrelProjection;
            } // end of generic time if
            else { // start of else

              secondBestECalProjection = barrelProjection;
            } // end of else
          } // end of barrel projection if

        } // end of eCalBarrelInnerR if

        trackState_AtCalorimeter = getExtrapolationAtCalorimeter(bestECalProjection, helixAtLastHit, Bz);

      } // end of eCalBarrel||eCalEndCap if
      break;
    } // end of state.location if
  } // end of state loop

  return trackState_AtCalorimeter;
}

double deltaPhi(const double p1, const double p2) {
  double oo2pi = 1. / (2. * M_PI);

  if (std::abs(p1 - p2) <= M_PI)
    return p1 - p2;

  double n = std::round((p1 - p2) * oo2pi);

  return p1 - p2 - n * 2. * M_PI;
}

double getFieldFromCompact() {

  dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
  const double position[3] = {0, 0, 0};      // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3] = {0, 0, 0}; // initialise object to hold magnetic field
  mainDetector.field().magneticField(position, magneticFieldVector); // get the magnetic field vector from DD4hep
  return magneticFieldVector[2] / dd4hep::tesla;                     // z component at (0,0,0)
}

ThreeVector getMomentum(const edm4hep::TrackState& state, double Bz) {

  double phi = state.phi;
  double omega = state.omega;
  double tanLambda = state.tanLambda;

  double pt = gevPerTeslaMm * Bz / std::abs(omega);

  double px = pt * std::cos(phi);
  double py = pt * std::sin(phi);
  double pz = pt * tanLambda;

  return ThreeVector(px, py, pz);
}

bool GetPointInZ(const double zPlane,       // z plane (endcap) coordinate
                 const ThreeVector& refPt,  // reference point
                 const ThreeVector& momVec, // momentum
                 ThreeVector& intersectPt) {

  // protect from zero division
  if (std::abs(momVec.z) < std::numeric_limits<float>::epsilon()) {
    throw std::invalid_argument("BremRecoveryUtility: invalid momentum.z (zero division)");
  }

  // check the assumption
  if (std::abs(refPt.z) > std::abs(zPlane)) {
    throw std::invalid_argument("BremRecoveryUtility: invalid refPoint.z (ref. point outside the endcap)");
  }

  intersectPt.z = momVec.z > 0. ? std::abs(zPlane) : -std::abs(zPlane);
  const double genericTime = (intersectPt.z - refPt.z) / momVec.z;
  intersectPt.x = refPt.x + momVec.x * genericTime;
  intersectPt.y = refPt.y + momVec.y * genericTime;

  return true;
}

bool GetPointOnCircle(const double radius,       // (barrel) radius
                      const ThreeVector& refPt,  // reference point
                      const ThreeVector& momVec, // momentum
                      ThreeVector& intersectPt) {

  double refPtXY = refPt.magXY();
  double momVecXY = momVec.magXY();
  // protect from zero division
  if (std::abs(momVecXY) < std::numeric_limits<float>::epsilon()) {
    throw std::invalid_argument("BremRecoveryUtility: invalid pt (zero division)");
  }

  // check the assumption
  if (std::abs(refPtXY) > radius) {
    throw std::invalid_argument("BremRecoveryUtility: invalid refPoint (ref. point outside the barrel)");
  }
  // suppose a triangle OAP in the XY plane
  //      A        A is the reference point
  //      /\       O is the detector origin
  //     /  \      P is the intersection point
  //    ------     then the length OP = radius
  //   O      P

  // now do some middle school Euclidean geometry
  // cos(O+P) = refPt*momVec / |refPt|*|momVec|
  double normedInnerProd = (refPt.x * momVec.x + refPt.y * momVec.y) / (refPtXY * momVecXY);
  // sine rule : radius/sin(P+O) = |refPt|/sin(P)
  double sinP = refPtXY * std::sqrt(1. - normedInnerProd * normedInnerProd) / radius;
  double angOplusP = std::acos(normedInnerProd);
  double angP = std::asin(sinP); // always smaller than pi/2 by definition
  double angO = angOplusP - angP;
  // sine rule again : lenAP/sin(O) = |refPt|/sin(P)
  double lenAP = refPtXY * std::sin(angO) / sinP;
  double genericTime = lenAP / momVecXY; // this is what we wanted

  // now done with the geometry, come back to 3D
  intersectPt.x = refPt.x + momVec.x * genericTime;
  intersectPt.y = refPt.y + momVec.y * genericTime;
  intersectPt.z = refPt.z + momVec.z * genericTime;

  return true;
}

/* if track reference at the origin */
bool GetPointOnCircleForRefAtOrigin(const double radius,       // (barrel) radius
                                    const ThreeVector& refPt,  // reference point
                                    const ThreeVector& momVec, // momentum
                                    ThreeVector& intersectPt) {

  double refPtXY = refPt.magXY();
  double momVecXY = momVec.magXY();
  // protect from zero division
  if (std::abs(momVecXY) < std::numeric_limits<float>::epsilon()) {
    throw std::invalid_argument("BremRecoveryUtility: invalid pt (zero division)");
  }

  // check the assumption
  if (std::abs(refPtXY) > radius) {
    throw std::invalid_argument("BremRecoveryUtility: invalid refPoint (ref. point outside the barrel)");
  }
  double genericTime = radius / momVecXY; // this is what we wanted

  // now done with the geometry, come back to 3D
  intersectPt.x = refPt.x + momVec.x * genericTime;
  intersectPt.y = refPt.y + momVec.y * genericTime;
  intersectPt.z = refPt.z + momVec.z * genericTime;

  return true;
}
