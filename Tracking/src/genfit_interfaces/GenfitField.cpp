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

#include "GenfitField.h"

namespace GenfitInterface {
    
    GenfitField::GenfitField(dd4hep::OverlayedField dd4hepField)
    :m_dd4hepField(dd4hepField)
    {
        genfit::FieldManager::getInstance()->init(this);
    }

    TVector3 GenfitField::get(const TVector3& pos) const
    {
        double B[3]={1e9,1e9,1e9};
        get(pos.X(),pos.Y(),pos.Z(),B[0],B[1],B[2]);
        return TVector3(B[0],B[1],B[2]);
    }

    
    void GenfitField::get(const double& posX, const double& posY, const double& posZ,
            double& Bx, double& By, double& Bz) const
    {
        /// get field from dd4hepField
        const dd4hep::Direction& field=m_dd4hepField.magneticField(
                dd4hep::Position(
                posX* dd4hep::mm,
                posY* dd4hep::mm,
                posZ* dd4hep::mm));

        Bx=field.X() / dd4hep::kilogauss;
        By=field.Y() / dd4hep::kilogauss;
        Bz=field.Z() / dd4hep::kilogauss;

    }

    double GenfitField::getBz(const TVector3& pos) const
    {
        return get(pos).Z();
    }

}