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

#include "BremsRecoveryUtils.h"

#include "DD4hep/Detector.h"

#include "edm4hep/utils/vector_utils.h"

#include <cmath>
#include <limits>
#include <stdexcept>

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

edm4hep::Vector3f GetPointInZ(const double zPlane, const edm4hep::Vector3f& refPt, const edm4hep::Vector3f& momVec) {
  // check the assumption
  if (std::abs(refPt.z) > std::abs(zPlane)) {
    throw std::invalid_argument("BremsRecoveryUtils: invalid refPoint.z (ref. point outside the endcap)");
  }

  const double z = momVec.z > 0. ? std::abs(zPlane) : -std::abs(zPlane);
  return refPt + ((z - refPt.z) / momVec.z) * momVec;
}

edm4hep::Vector3f GetPointOnCircle(const double radius, const edm4hep::Vector3f& refPt,
                                   const edm4hep::Vector3f& momVec) {
  const double refPtXY = edm4hep::utils::magnitudeTransverse(refPt);
  const double momVecXY = edm4hep::utils::magnitudeTransverse(momVec);
  // protect from zero division
  if (std::abs(momVecXY) < std::numeric_limits<float>::epsilon()) {
    throw std::invalid_argument("BremsRecoveryUtils: invalid pt (zero division)");
  }
  // check the assumption
  if (std::abs(refPtXY) > radius) {
    throw std::invalid_argument("BremsRecoveryUtils: invalid refPoint (ref. point outside the barrel)");
  }
  // suppose a triangle OAP in the XY plane
  //      A        A is the reference point
  //      /\       O is the detector origin
  //     /  \      P is the intersection point
  //    ------     then the length OP = radius

  // cos(O+P) = refPt*momVec / |refPt|*|momVec|, in the transverse plane
  const double normedInnerProd =
      edm4hep::utils::vectorTransverse(refPt) * edm4hep::utils::vectorTransverse(momVec) / (refPtXY * momVecXY);
  // sine rule : radius/sin(P+O) = |refPt|/sin(P)
  const double sinP = refPtXY * std::sqrt(1. - normedInnerProd * normedInnerProd) / radius;
  const double angOplusP = std::acos(normedInnerProd);
  const double angP = std::asin(sinP); // always smaller than pi/2 by definition
  const double angO = angOplusP - angP;
  // sine rule again : lenAP/sin(O) = |refPt|/sin(P)
  const double lenAP = refPtXY * std::sin(angO) / sinP;

  return refPt + (lenAP / momVecXY) * momVec;
}

edm4hep::Vector3f GetPointOnCircleForRefAtOrigin(const double radius, const edm4hep::Vector3f& refPt,
                                                 const edm4hep::Vector3f& momVec) {
  const double refPtXY = edm4hep::utils::magnitudeTransverse(refPt);
  const double momVecXY = edm4hep::utils::magnitudeTransverse(momVec);
  // protect from zero division
  if (std::abs(momVecXY) < std::numeric_limits<float>::epsilon()) {
    throw std::invalid_argument("BremsRecoveryUtils: invalid pt (zero division)");
  }
  // check the assumption
  if (std::abs(refPtXY) > radius) {
    throw std::invalid_argument("BremsRecoveryUtils: invalid refPoint (ref. point outside the barrel)");
  }

  return refPt + (radius / momVecXY) * momVec;
}
