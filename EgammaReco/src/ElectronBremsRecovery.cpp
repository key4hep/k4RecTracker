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

#include "ElectronBremsRecovery.h"

#include "BremsRecoveryUtils.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"

#include "edm4hep/utils/vector_utils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

DECLARE_COMPONENT(ElectronBremsRecovery)

ElectronBremsRecovery::ElectronBremsRecovery(const std::string& name, ISvcLocator* svc)
    : Transformer(name, svc, {KeyValue{"InputRecoParticles", "PandoraPFOs"}},
                  {KeyValue{"OutputRecoParticles", "RecoveredPFOs"}}) {}

StatusCode ElectronBremsRecovery::initialize() {
  StatusCode sc = Transformer::initialize();
  if (sc.isFailure())
    return sc;

  const auto* detector = m_geoSvc->getDetector();
  try {
    m_eCalBarrelInnerR = detector->constantAsDouble(m_eCalBarrelInnerRName) / dd4hep::mm;
    m_eCalEndCapInnerZ = detector->constantAsDouble(m_eCalEndCapInnerZName) / dd4hep::mm;
  } catch (const std::exception& e) {
    error() << "ECAL constant not found in the geometry (" << e.what() << "): set " << m_eCalBarrelInnerRName.name()
            << " and " << m_eCalEndCapInnerZName.name() << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_eCalBarrelInnerR <= 0.) {
    error() << "ECAL barrel inner radius must be positive, got " << m_eCalBarrelInnerR << " mm" << endmsg;
    return StatusCode::FAILURE;
  }
  debug() << "ECAL barrel inner radius = " << m_eCalBarrelInnerR << " mm, endcap inner z = " << m_eCalEndCapInnerZ
          << " mm" << endmsg;

  m_Bz = getFieldFromCompact();
  return StatusCode::SUCCESS;
}

edm4hep::ReconstructedParticleCollection
ElectronBremsRecovery::operator()(const edm4hep::ReconstructedParticleCollection& input) const {
  std::vector<ClusterEntry> clusters;
  collectClusters(input, clusters);

  BremsMatch brems;
  findBremsPfos(input, clusters, brems);

  edm4hep::ReconstructedParticleCollection output;
  buildOutput(input, brems, output);
  return output;
}

void ElectronBremsRecovery::collectClusters(const edm4hep::ReconstructedParticleCollection& pfos,
                                            std::vector<ClusterEntry>& clusters) const {
  for (std::size_t pfoIndex = 0; pfoIndex < pfos.size(); ++pfoIndex) {
    const auto& pfo = pfos[pfoIndex];
    debug() << "PFO Energy = " << pfo.getEnergy() << " Tracks = " << pfo.getTracks().size()
            << " Clusters = " << pfo.getClusters().size() << endmsg;
    for (const auto& cluster : pfo.getClusters()) {
      if (cluster.isAvailable())
        clusters.push_back({cluster, pfoIndex});
    }
  }
  // the most energetic electrons are matched first
  std::ranges::sort(clusters, std::ranges::greater{},
                    [](const ClusterEntry& entry) { return entry.cluster.getEnergy(); });
}

const edm4hep::Track& ElectronBremsRecovery::closestTrack(const edm4hep::ReconstructedParticle& pfo,
                                                          const edm4hep::Cluster& cluster) const {
  const auto clusterPos = cluster.getPosition();
  const auto distance = [&](const edm4hep::Track& track) {
    // the track state at the calorimeter is filled by the track fit; a track without it is taken last
    const auto state = track.isAvailable() ? track.getTrackState(edm4hep::TrackState::AtCalorimeter) : std::nullopt;
    if (!state.has_value()) {
      debug() << "Track without a track state at the calorimeter" << endmsg;
      return std::numeric_limits<double>::max();
    }
    return edm4hep::utils::magnitude(clusterPos - state->referencePoint);
  };
  return *std::ranges::min_element(pfo.getTracks(), std::ranges::less{}, distance);
}

void ElectronBremsRecovery::bremsTangentAtCalorimeter(const edm4hep::TrackState& stateAtIP,
                                                      edm4hep::Vector3f& tangentPt) const {
  const auto& refPt = stateAtIP.referencePoint;
  const auto momVec = getMomentum(stateAtIP, m_Bz);

  // first the endcap plane
  GetPointInZ(m_eCalEndCapInnerZ, refPt, momVec, tangentPt);
  // the tangent reaches the barrel before the endcap
  if (edm4hep::utils::magnitudeTransverse(tangentPt) > m_eCalBarrelInnerR) {
    const bool refAtOrigin = std::abs(refPt.x) < std::numeric_limits<float>::epsilon() &&
                             std::abs(refPt.y) < std::numeric_limits<float>::epsilon() &&
                             std::abs(refPt.z) < std::numeric_limits<float>::epsilon();
    if (refAtOrigin) {
      GetPointOnCircleForRefAtOrigin(m_eCalBarrelInnerR, refPt, momVec, tangentPt);
    } else {
      GetPointOnCircle(m_eCalBarrelInnerR, refPt, momVec, tangentPt);
    }
  }
}

void ElectronBremsRecovery::findBremsPfos(const edm4hep::ReconstructedParticleCollection& pfos,
                                          const std::vector<ClusterEntry>& clusters, BremsMatch& brems) const {
  const auto hasTrack = [](const edm4hep::ReconstructedParticle& pfo) {
    return std::ranges::any_of(pfo.getTracks(), [](const edm4hep::Track& track) { return track.isAvailable(); });
  };

  for (const auto& seed : clusters) {
    // only an electron (cluster with a track) seeds the brems recovery
    const auto& electron = pfos[seed.pfoIndex];
    if (std::abs(electron.getPDG()) != 11 || !hasTrack(electron))
      continue;

    edm4hep::TrackState stateAtIP;
    if (const auto state = closestTrack(electron, seed.cluster).getTrackState(edm4hep::TrackState::AtIP);
        state.has_value()) {
      stateAtIP = state.value();
    }
    edm4hep::Vector3f tangentPt;
    bremsTangentAtCalorimeter(stateAtIP, tangentPt);

    // brems window: the seed eta (+- dEtaMargin) and the phi band between the seed and the tangent (+- dPhiMargin)
    const auto& seedPos = seed.cluster.getPosition();
    const double seedEta = edm4hep::utils::eta(seedPos);
    const double seedPhi = edm4hep::utils::angleAzimuthal(seedPos);
    const double bremsBandPhi = deltaPhi(edm4hep::utils::angleAzimuthal(tangentPt), seedPhi);
    const double minDPhi = bremsBandPhi > 0. ? -m_dPhiMargin : bremsBandPhi - m_dPhiMargin;
    const double maxDPhi = bremsBandPhi > 0. ? bremsBandPhi + m_dPhiMargin : m_dPhiMargin.value();

    for (const auto& candidate : clusters) {
      // brems candidates: clusters of neutral (track-less) PFOs, as the whole PFO is absorbed by the electron
      if (hasTrack(pfos[candidate.pfoIndex]))
        continue;

      const auto& candidatePos = candidate.cluster.getPosition();
      const double dEta = edm4hep::utils::eta(candidatePos) - seedEta;
      const double dPhi = deltaPhi(edm4hep::utils::angleAzimuthal(candidatePos), seedPhi);
      if (std::abs(dEta) < m_dEtaMargin && minDPhi < dPhi && dPhi < maxDPhi) {
        debug() << "Found a brems candidate" << endmsg;
        // a neutral PFO is absorbed by at most one electron (the most energetic one, which comes first)
        if (brems.absorbed.insert(candidate.pfoIndex).second) {
          brems.absorbedBy[seed.pfoIndex].push_back(candidate.pfoIndex);
        }
      }
    }
  }
}

void ElectronBremsRecovery::buildOutput(const edm4hep::ReconstructedParticleCollection& pfos, const BremsMatch& brems,
                                        edm4hep::ReconstructedParticleCollection& output) {
  for (std::size_t pfoIndex = 0; pfoIndex < pfos.size(); ++pfoIndex) {
    // a neutral PFO identified as brems is now part of its electron: storing it too would double count
    if (brems.absorbed.contains(pfoIndex))
      continue;

    // every other PFO is copied as is: energy, momentum, charge, PDG, mass and its tracks and clusters
    auto particle = pfos[pfoIndex].clone();

    // an electron that recovered brems
    if (const auto it = brems.absorbedBy.find(pfoIndex); it != brems.absorbedBy.end()) {
      for (const auto absorbed : it->second) {
        // takes over the clusters of the absorbed neutral PFO
        for (const auto& cluster : pfos[absorbed].getClusters()) {
          particle.addToClusters(cluster);
        }
        // and its energy, so that the total energy of the event is unchanged. The momentum is not touched:
        // it comes from the track, which measures the electron before it radiated
        particle.setEnergy(particle.getEnergy() + pfos[absorbed].getEnergy());
      }
    }
    output.push_back(particle);
  }
}
