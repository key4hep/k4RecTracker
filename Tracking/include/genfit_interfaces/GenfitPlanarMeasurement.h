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

#ifndef PLANAR_MEASUREMENT_H
#define PLANAR_MEASUREMENT_H

// standard 
#include <iostream>

// ROOT
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TMatrixD.h>
#include "TDecompSVD.h"

// GENFIT
#include <PlanarMeasurement.h>
#include <TrackPoint.h>

// EDM4hep and DDRec
#include "edm4hep/TrackerHitPlane.h"
#include "DDRec/Vector3D.h" 
#include "DDRec/Vector2D.h"
#include "DD4hep/DD4hepUnits.h"


/** @class PlanarMeasurement
 *
 *  Interface class that wraps a TrackerHitPlane object into a GENFIT-compatible PlanarMeasurement.
 *  
 *  This class converts a hit from the EDM4hep TrackerHitPlane collection into a `genfit::PlanarMeasurement`,
 *  allowing it to be used as input for GENFIT's track fitting algorithms.
 *  The constructor initializes the measurement from the given hit and stores the detector and hit indices
 *  for debugging or traceability purposes.
 *
 *
 *  Author: Andrea De Vita  
 *  Date  : 2025-06
 *
*/

namespace GenfitInterface {

    class PlanarMeasurement {
    public:

        // Constructor
        PlanarMeasurement(const edm4hep::TrackerHitPlane& hit, const int det_idx, const int hit_idx,const int debug_lvl);

        // Getter for genfit::WirePointMeasurement
        genfit::PlanarMeasurement* getGenFit() const {return m_genfitHit;};

    private:
    
        genfit::PlanarMeasurement* m_genfitHit;
    };

} 

#endif // PLANAR_MEASUREMENT_H