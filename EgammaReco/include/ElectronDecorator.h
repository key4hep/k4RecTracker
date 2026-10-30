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

#ifndef ELECTRON_DECORATOR_H
#define ELECTRON_DECORATOR_H

#include "GaudiKernel/ServiceHandle.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"

#include "edm4hep/ReconstructedParticleCollection.h"

class ElectronDecorator : public k4FWCore::Transformer<edm4hep::ReconstructedParticleCollection(
                              const edm4hep::ReconstructedParticleCollection&)> {

public:
  ElectronDecorator(const std::string& name, ISvcLocator* svc);

  StatusCode initialize() override;

  edm4hep::ReconstructedParticleCollection operator()(const edm4hep::ReconstructedParticleCollection&) const override;

private:
  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Geometry service"};

  Gaudi::Property<double> m_dEtaMargin{this, "dEtaMargin", 0.005, "Eta matching margin"};
  Gaudi::Property<double> m_dPhiMargin{this, "dPhiMargin", 0.05, "Phi matching margin (brem spreads along phi)"};

  // names of the compact-file constants describing the ECAL front face (defaults: CLD)
  Gaudi::Property<std::string> m_eCalBarrelInnerRName{this, "ECalBarrelInnerRadiusName", "ECalBarrel_inner_radius",
                                                      "Constant with the ECAL barrel inner radius"};
  Gaudi::Property<std::string> m_eCalBarrelMaxZName{this, "ECalBarrelHalfLengthName", "ECalBarrel_half_length",
                                                    "Constant with the ECAL barrel half length"};
  Gaudi::Property<std::string> m_eCalEndCapInnerRName{this, "ECalEndcapInnerRadiusName", "ECalEndcap_inner_radius",
                                                      "Constant with the ECAL endcap inner radius"};
  Gaudi::Property<std::string> m_eCalEndCapOuterRName{this, "ECalEndcapOuterRadiusName", "ECalEndcap_outer_radius",
                                                      "Constant with the ECAL endcap outer radius"};
  Gaudi::Property<std::string> m_eCalEndCapInnerZName{this, "ECalEndcapMinZName", "ECalEndcap_min_z",
                                                      "Constant with the ECAL endcap inner z"};

  // ECAL front face [mm] and Bz [T], read from the geometry in initialize()
  double m_eCalBarrelInnerR{0.};
  double m_eCalBarrelMaxZ{0.};
  double m_eCalEndCapInnerR{0.};
  double m_eCalEndCapOuterR{0.};
  double m_eCalEndCapInnerZ{0.};
  double m_Bz{0.};
};

#endif // ELECTRON_DECORATOR_H
