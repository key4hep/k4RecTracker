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

#ifndef EGAMMARECO_BREMSRECOVERYUTILS_H
#define EGAMMARECO_BREMSRECOVERYUTILS_H

#include "edm4hep/Vector3f.h"

/* p1 - p2, wrapped into [-pi, pi] */
double deltaPhi(const double p1, const double p2);

/* Bz [T] at the origin, from the DD4hep geometry */
double getFieldFromCompact();

/* where the line refPt + t * momVec (t > 0) crosses the endcap plane z = +-zPlane */
edm4hep::Vector3f GetPointInZ(const double zPlane, const edm4hep::Vector3f& refPt, const edm4hep::Vector3f& momVec);

/* where the line refPt + t * momVec (t > 0) crosses the barrel cylinder of the given radius */
edm4hep::Vector3f GetPointOnCircle(const double radius, const edm4hep::Vector3f& refPt,
                                   const edm4hep::Vector3f& momVec);

/* same as GetPointOnCircle, for a reference point at the origin */
edm4hep::Vector3f GetPointOnCircleForRefAtOrigin(const double radius, const edm4hep::Vector3f& refPt,
                                                 const edm4hep::Vector3f& momVec);

#endif // EGAMMARECO_BREMSRECOVERYUTILS_H
