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

#include "GenfitWireMeasurement.h"

namespace GenfitInterface {

    WireMeasurement::WireMeasurement(const edm4hep::SenseWireHit& hit, const dd4hep::rec::DCH_info* dch_info, const dd4hep::DDSegmentation::BitFieldCoder* decoder,const int det_idx, const int hit_idx, const int debug_lvl) {


        // NB: dd4hep::mm = 0.1.
        // Therefore:
        //   - val * dd4hep::mm  converts millimeters (mm) to centimeters (cm)
        //   - val / dd4hep::mm  converts centimeters (cm) to millimeters (mm)

        // Cell and geometry parameters
        const int cellid = hit.getCellID();

        const double wireStereoAngle     = hit.getWireStereoAngle();
        const double wireAzimuthalAngle  = hit.getWireAzimuthalAngle();

        const auto hitPos = hit.getPosition();  // mm

        // Convert hit position to cm
        const TVector3 position(
            hitPos.x * dd4hep::mm,
            hitPos.y * dd4hep::mm,
            hitPos.z * dd4hep::mm
        );

        // Measurement uncertainties
        const double positionAlongWireError = hit.getPositionAlongWireError() * dd4hep::mm; // cm
        const double distanceToWire         = hit.getDistanceToWire()         * dd4hep::mm; // cm
        const double distanceToWireError    = hit.getDistanceToWireError()    * dd4hep::mm; // cm

        // Wire direction
        TVector3 direction(0., 0., 1.);
        direction.RotateX(wireStereoAngle);
        direction.RotateZ(wireAzimuthalAngle);
        direction = direction.Unit();

        // Wire extremities
        const int layer       = decoder->get(cellid, "layer");
        const int superlayer  = decoder->get(cellid, "superlayer");
        const int nphi        = decoder->get(cellid, "nphi");

        const int ilayer = dch_info->CalculateILayerFromCellIDFields(layer, superlayer);
        const auto& l    = dch_info->database.at(ilayer);

        const int stereosign = l.StereoSign();
        const double rz0     = l.radius_sw_z0;

        const double dphi  = dch_info->twist_angle;
        const double kappa = (1. / dch_info->Lhalf) * std::tan(dphi / 2.);

        TVector3 p1(
            rz0,
            -stereosign * rz0 * kappa * dch_info->Lhalf,
            -dch_info->Lhalf
        );

        TVector3 p2(
            rz0,
            stereosign * rz0 * kappa * dch_info->Lhalf,
            dch_info->Lhalf
        );

        const double phi_z0 = dch_info->Calculate_wire_phi_z0(ilayer, nphi);

        p1.RotateZ(phi_z0);
        p2.RotateZ(phi_z0);

        // Observables
        const double Rdrift = distanceToWire;           // cm
        const double zreco  = (position - p1).Mag();    // cm

        // Measurement
        TVectorD rawHitCoords(8);
        rawHitCoords[0] = p1.X();
        rawHitCoords[1] = p1.Y();
        rawHitCoords[2] = p1.Z();
        rawHitCoords[3] = p2.X();
        rawHitCoords[4] = p2.Y();
        rawHitCoords[5] = p2.Z();
        rawHitCoords[6] = Rdrift;
        rawHitCoords[7] = zreco;

        // Covariance matrix
        constexpr double sigmaWire = 1e-4;  // cm (almost fixed)

        TMatrixDSym rawHitCov(8);
        rawHitCov.Zero();

        // Diagonal terms
        rawHitCov(0,0) = sigmaWire * sigmaWire;
        rawHitCov(1,1) = sigmaWire * sigmaWire;
        rawHitCov(2,2) = sigmaWire * sigmaWire;
        rawHitCov(3,3) = sigmaWire * sigmaWire;
        rawHitCov(4,4) = sigmaWire * sigmaWire;
        rawHitCov(5,5) = sigmaWire * sigmaWire;

        rawHitCov(6,6) = distanceToWireError    * distanceToWireError;
        rawHitCov(7,7) = positionAlongWireError * positionAlongWireError;


        // Create the genfit::WirePointMeasurement
        m_genfitHit = new genfit::WirePointMeasurement(rawHitCoords, rawHitCov, det_idx, hit_idx, nullptr);

        if (debug_lvl > 0) {

            std::cout << "\n========== Wire Measurement Debug ==========\n";

            std::cout << "Hit global position [cm]      : ("
                    << position.X() << ", "
                    << position.Y() << ", "
                    << position.Z() << ")\n";

            std::cout << "Distance to wire [cm]         : "
                    << distanceToWire
                    << "  ±  " << distanceToWireError << "\n";

            std::cout << "Position along wire error [cm]: "
                    << positionAlongWireError << "\n";

            std::cout << "\nrawHitCoords (8D) [cm]: [ ";
            for (int i = 0; i < rawHitCoords.GetNrows(); ++i) {
                std::cout << rawHitCoords[i];
                if (i != rawHitCoords.GetNrows() - 1)
                    std::cout << ", ";
            }
            std::cout << " ]\n";

            std::cout << "\nUncertainties [cm]:\n";
            std::cout << "  Wire geometry sigma : " << sigmaWire                << "\n";
            std::cout << "  Drift distance sigma: " << distanceToWireError      << "\n";
            std::cout << "  z_reco sigma        : " << positionAlongWireError   << "\n";

            std::cout << "\nIdentifiers:\n";
            std::cout << "  detID : " << det_idx << "\n";
            std::cout << "  hitID : " << hit_idx << "\n";

            std::cout << "============================================\n\n";
        }

    }

} // namespace GENFIT