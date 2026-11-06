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

#ifndef ELECTRON_BREMSRECOVERY_H
#define ELECTRON_BREMSRECOVERY_H

#include "GaudiKernel/ServiceHandle.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"

#include "BremsRecoveryUtils.h"

#include "edm4hep/Cluster.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/Vector3f.h"

#include <cstddef>
#include <map>
#include <unordered_set>
#include <vector>

class ElectronBremsRecovery : public k4FWCore::Transformer<edm4hep::ReconstructedParticleCollection(
                                  const edm4hep::ReconstructedParticleCollection&)> {

public:
  ElectronBremsRecovery(const std::string& name, ISvcLocator* svc);

  StatusCode initialize() override;

  edm4hep::ReconstructedParticleCollection operator()(const edm4hep::ReconstructedParticleCollection&) const override;

private:
  // a cluster of an input PFO
  struct ClusterEntry {
    edm4hep::Cluster cluster;
    std::size_t pfoIndex; // of its PFO in the input collection
  };

  // result of the brems matching
  struct BremsMatch {
    std::map<std::size_t, std::vector<std::size_t>> absorbedBy; // electron PFO index -> absorbed neutral PFO indices
    std::unordered_set<std::size_t> absorbed;                   // all absorbed neutral PFO indices
  };

  // one entry per cluster of the input PFOs, sorted by decreasing cluster energy
  std::vector<ClusterEntry> collectClusters(const edm4hep::ReconstructedParticleCollection& pfos) const;
  // the track of the PFO whose state at the calorimeter is closest to the cluster
  const edm4hep::Track& closestTrack(const edm4hep::ReconstructedParticle& pfo, const edm4hep::Cluster& cluster) const;
  // where the tangent to the track at the IP reaches the ECAL front face (endcap plane or barrel cylinder)
  edm4hep::Vector3f bremsTangentAtCalorimeter(const edm4hep::TrackState& stateAtIP) const;
  // match the clusters of neutral PFOs to electrons, inside the eta-phi brems window of each electron
  BremsMatch findBremsPfos(const edm4hep::ReconstructedParticleCollection& pfos,
                           const std::vector<ClusterEntry>& clusters) const;
  // copy the input PFOs, adding the absorbed clusters and energy to their electron and dropping the absorbed PFOs
  static edm4hep::ReconstructedParticleCollection buildOutput(const edm4hep::ReconstructedParticleCollection& pfos,
                                                              const BremsMatch& brems);

  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Geometry service"};

  Gaudi::Property<double> m_dEtaMargin{this, "dEtaMargin", 0.005, "Eta matching margin"};
  Gaudi::Property<double> m_dPhiMargin{this, "dPhiMargin", 0.05, "Phi matching margin (brems spreads along phi)"};

  // names of the compact-file constants describing the ECAL front face (defaults: CLD)
  Gaudi::Property<std::string> m_eCalBarrelInnerRName{this, "ECalBarrelInnerRadiusName", "ECalBarrel_inner_radius",
                                                      "Constant with the ECAL barrel inner radius"};
  Gaudi::Property<std::string> m_eCalEndCapInnerZName{this, "ECalEndcapMinZName", "ECalEndcap_min_z",
                                                      "Constant with the ECAL endcap inner z"};

  // ECAL front face [mm] and Bz [T], read from the geometry in initialize()
  double m_eCalBarrelInnerR{0.};
  double m_eCalEndCapInnerZ{0.};
  double m_Bz{0.};
};

#endif // ELECTRON_BREMSRECOVERY_H
