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


#ifndef GENFIT_FIELD_H
#define GENFIT_FIELD_H

// DD4hep
#include "DD4hep/Fields.h"
#include "DD4hep/DD4hepUnits.h"

//genfit
#include "FieldManager.h"
#include "ConstField.h"
#include "AbsBField.h"

class TVector3;


/** @class GenfitField
 *
 *  Adapter class that provides a magnetic field interface compatible with GENFIT,
 *  using a `dd4hep::OverlayedField` as the underlying field source.
 *
 *  The `GenfitField` implements the `genfit::AbsBField` interface, allowing GENFIT to query
 *  the magnetic field vector at any point in space. It serves as a bridge between the DD4hep
 *  field service and the GENFIT tracking library.
 *
 *  Field values can be retrieved either as a full 3D vector (`TVector3`) or as individual
 *  components (`Bx`, `By`, `Bz`) at a given position. A convenience method is also provided
 *  to directly access the `Bz` component in kilogauss.
 *
 *
 *  Author: Andrea De Vita  
 *  Date  : 2025-11
 *
 */


namespace GenfitInterface {

    class GenfitField : public genfit::AbsBField{
        public:
            GenfitField(dd4hep::OverlayedField dd4hepField);
            virtual ~GenfitField(){;}

            // Get field value from TVector3
            TVector3 get(const TVector3& pos) const override;

            // Get field value from doubles
            void get(const double& posX, const double& posY, const double& posZ,
                    double& Bx, double& By, double& Bz) const override;

            // Get Bz, kilogauss
            double getBz(const TVector3& pos) const;

            // Get dd4hep field
            const dd4hep::OverlayedField field() const {return m_dd4hepField;}

        private:
            dd4hep::OverlayedField m_dd4hepField;
    };
}

#endif