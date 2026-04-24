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

#include "GenfitPlanarMeasurement.h"

namespace GenfitInterface {

PlanarMeasurement::PlanarMeasurement(const edm4hep::TrackerHitPlane& hit, const int det_idx, const int hit_idx,
                                     const int debug_lvl = 0) {

  // NB: dd4hep::mm = 0.1.
  // Therefore:
  //   - val * dd4hep::mm  converts millimeters (mm) to centimeters (cm)
  //   - val / dd4hep::mm  converts centimeters (cm) to millimeters (mm)

  // Global position (convert mm to cm)
  const auto hitPos = hit.getPosition(); // input in mm
  const TVector3 globalPos(hitPos.x * dd4hep::mm, hitPos.y * dd4hep::mm,
                           hitPos.z * dd4hep::mm); // cm

  const TVector3 Origin = globalPos;

  // Local 2D coordinates (initialized to zero, in cm)
  dd4hep::rec::Vector2D localPos{0., 0.};

  TVectorD rawHitCoords(2);
  rawHitCoords[0] = localPos[0]; // cm
  rawHitCoords[1] = localPos[1]; // cm

  // Construct local U and V unit vectors from spherical coordinates
  const auto U_theta = hit.getU()[0];
  const auto U_phi = hit.getU()[1];
  TVector3 U(std::sin(U_theta) * std::cos(U_phi), std::sin(U_theta) * std::sin(U_phi), std::cos(U_theta));
  U = U.Unit();

  const auto V_theta = hit.getV()[0];
  const auto V_phi = hit.getV()[1];
  TVector3 V(std::sin(V_theta) * std::cos(V_phi), std::sin(V_theta) * std::sin(V_phi), std::cos(V_theta));
  V = V.Unit();

  // Measurement uncertainties
  const double sigmaU = hit.getDu() * dd4hep::mm; // cm
  const double sigmaV = hit.getDv() * dd4hep::mm; // cm

  TMatrixDSym rawHitCov(2);
  rawHitCov(0, 0) = sigmaU * sigmaU; // cm^2
  rawHitCov(0, 1) = 0.0;
  rawHitCov(1, 0) = 0.0;
  rawHitCov(1, 1) = sigmaV * sigmaV; // cm^2

  auto cellID0 = hit.getCellID();

  if (debug_lvl > 0) {

    std::cout << "\n========== Planar Measurement Debug ==========\n";

    std::cout << "Global position [cm] : (" << globalPos.X() << ", " << globalPos.Y() << ", " << globalPos.Z() << ")\n";

    std::cout << "Origin [cm]          : (" << Origin.X() << ", " << Origin.Y() << ", " << Origin.Z() << ")\n";

    std::cout << "U direction          : (" << U.X() << ", " << U.Y() << ", " << U.Z() << ")\n";

    std::cout << "V direction          : (" << V.X() << ", " << V.Y() << ", " << V.Z() << ")\n\n";

    std::cout << "Local coordinates [cm]:\n";
    rawHitCoords.Print();

    std::cout << "Covariance matrix [cm^2]:\n";
    rawHitCov.Print();

    std::cout << "\nIdentifiers:\n";
    std::cout << "  detId   : " << det_idx << "\n";
    std::cout << "  hitID   : " << hit_idx << "\n";
    std::cout << "  planeID : " << cellID0 << "\n";

    std::cout << "==============================================\n\n";
  }

  // Create genfit::PlanarMeasurement

  // This pointer is used to create a genfit::TrackPoint, which is then stored
  // inside a genfit::Track. The Track aggregates a collection of TrackPoints
  // representing the hits used for the track fit.
  //
  // GENFIT manages ownership of the full chain:
  // - A genfit::Track owns its TrackPoints
  // - When a Track is deleted, all associated TrackPoints are automatically deleted
  // - When a TrackPoint is deleted, the associated measurement object used to
  //   construct it is also released
  m_genfitHit = new genfit::PlanarMeasurement(rawHitCoords, rawHitCov, det_idx, hit_idx, nullptr);
  genfit::SharedPlanePtr plane(new genfit::DetPlane(Origin, U, V));

  // Add plane
  m_genfitHit->setPlane(plane, cellID0);
}

} // namespace GenfitInterface