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

#include "ElectronDecorator.h"

#include "BremRecoveryUtils.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_set>
#include <vector>

DECLARE_COMPONENT(ElectronDecorator)

// struct for storing cluster and assocated track from PFO
struct ClusterTrackMap {
  edm4hep::Cluster cluster;
  std::vector<edm4hep::Track> tracks;
  int pdg;              // of the PFO the cluster belongs to
  std::size_t pfoIndex; // index of that PFO in the input collection
};

ElectronDecorator::ElectronDecorator(const std::string& name, ISvcLocator* svc)
    : Transformer(name, svc, {KeyValue{"InputRecoParticles", "PandoraPFOs"}},
                  {KeyValue{"OutputRecoParticles", "RecoveredPFOs"}}) {}

StatusCode ElectronDecorator::initialize() {
  StatusCode sc = Transformer::initialize();
  if (sc.isFailure())
    return sc;

  auto* detector = m_geoSvc->getDetector();
  const auto readConstant = [&](const Gaudi::Property<std::string>& name, double& value) {
    try {
      value = detector->constantAsDouble(name.value()) / dd4hep::mm;
    } catch (const std::exception&) {
      error() << "Constant '" << name.value() << "' not found in the geometry: set the property " << name.name()
              << endmsg;
      return false;
    }
    debug() << name.value() << " = " << value << " mm" << endmsg;
    return true;
  };
  if (!(readConstant(m_eCalBarrelInnerRName, m_eCalBarrelInnerR) &&
        readConstant(m_eCalBarrelMaxZName, m_eCalBarrelMaxZ) &&
        readConstant(m_eCalEndCapInnerRName, m_eCalEndCapInnerR) &&
        readConstant(m_eCalEndCapOuterRName, m_eCalEndCapOuterR) &&
        readConstant(m_eCalEndCapInnerZName, m_eCalEndCapInnerZ)))
    return StatusCode::FAILURE;

  m_Bz = getFieldFromCompact();
  return StatusCode::SUCCESS;
}

edm4hep::ReconstructedParticleCollection
ElectronDecorator::operator()(const edm4hep::ReconstructedParticleCollection& input) const {

  // Forming the cluster->track map
  std::vector<ClusterTrackMap> clusterTrackMaps;

  // brem candidates: clusters of neutral (track-less) PFOs, as the whole PFO is absorbed by the electron
  std::vector<ClusterTrackMap> tempClusters;

  //////////////////////////////////////////////////////////////////////////////////
  /////////////////Step1:loop through input PFOs ///////////////////////////////////
  /////////////////Step2:Unpack the PFOs to extract clusters and tracks ////////////
  /////////////////Step3:Make a map clusters: tracks ///////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////

  for (std::size_t pfoIndex = 0; pfoIndex < input.size(); ++pfoIndex) { // PFO loop starts
    const auto& pfo = input[pfoIndex];

    debug() << "PFO Energy = " << pfo.getEnergy() << " Tracks = " << pfo.getTracks().size()
            << " Clusters = " << pfo.getClusters().size() << endmsg;
    for (const auto& cluster : pfo.getClusters()) { // PFO cluster loop starts

      if (!cluster.isAvailable())
        continue;

      ClusterTrackMap entry;
      entry.cluster = cluster;
      entry.pdg = pfo.getPDG();
      entry.pfoIndex = pfoIndex;

      for (const auto& track : pfo.getTracks()) { // PFO track loop starts
        if (track.isAvailable()) {
          entry.tracks.push_back(track);
        }
      } //  PFO track loop ends

      if (entry.tracks.empty())
        tempClusters.push_back(entry);
      clusterTrackMaps.push_back(std::move(entry));
    } //  PFO cluster loop starts ends

  } // PFO  loop ends

  //////////////////////////////////////////////////////////////////////////////////
  /////////////////Step 1:complete /////////////////////////////////////////////////
  /////////////////Step 2:complete /////////////////////////////////////////////////
  /////////////////Step 3:complete /////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////
  /////////////////Step 4:Sort Cluster:Tracks map wrt cluster Energy ///////////////
  /////////////////Step 5:Sort Temporary Cluster Collection wrt cluster Energy /////
  //////////////////////////////////////////////////////////////////////////////////

  // sorting the cluster track map
  std::sort(clusterTrackMaps.begin(), clusterTrackMaps.end(), [](const ClusterTrackMap& a, const ClusterTrackMap& b) {
    return a.cluster.getEnergy() > b.cluster.getEnergy();
  });

  std::sort(tempClusters.begin(), tempClusters.end(), [](const ClusterTrackMap& a, const ClusterTrackMap& b) {
    return a.cluster.getEnergy() > b.cluster.getEnergy();
  });

  //////////////////////////////////////////////////////////////////////////////////
  /////////////////Step 4:complete /////////////////////////////////////////////////
  /////////////////Step 5:complete /////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////
  /////////////////Step 6:Extrapolate all tracks to clusters with utils function////
  /////////////////Step 7:Arrange the track wrt closest cluster to track distance //
  //////////////////////////////////////////////////////////////////////////////////

  for (auto& entry : clusterTrackMaps) { //  Looping over clusterTrackMaps starts

    if (entry.tracks.size() <= 1)
      continue;

    const auto clusterPos = entry.cluster.getPosition();

    std::sort(entry.tracks.begin(), entry.tracks.end(), [&](const edm4hep::Track& a, const edm4hep::Track& b) {
      auto getDistance = [&](const edm4hep::Track& track) {
        edm4hep::TrackState state = helixAtLastTrack(track, m_eCalBarrelInnerR, m_eCalBarrelMaxZ, m_eCalEndCapInnerR,
                                                     m_eCalEndCapOuterR, m_eCalEndCapInnerZ, m_Bz);

        auto ref = state.referencePoint;
        float dx = clusterPos.x - ref.x;
        float dy = clusterPos.y - ref.y;
        float dz = clusterPos.z - ref.z;

        return std::sqrt(dx * dx + dy * dy + dz * dz);
      };

      return getDistance(a) < getDistance(b); // closest first
    });
  } // Looping over clusterTrackMaps ends

  //////////////////////////////////////////////////////////////////////////////////
  /////////////////Step 6:complete /////////////////////////////////////////////////
  /////////////////Step 7:complete /////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////Step 8: Application of brem interaction finding method (interaction of tangent from track curve to cluster)////
  ////////////8:a) Look for brem interaction point in EndCap Region  ///////////////////////////////////////////////////
  ////////////8:b) check brem interaction point in EndCap Region has x-y radius>barrel radius///////////////////////////
  ////////////8:c) if x-y radius>barrel radius look for brem interaction in barrel region //////////////////////////////
  ///////Step 9: Identify brem clusters from brem interaction point and its eta-phi margin /////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::map<std::size_t, std::vector<std::size_t>> bremPfos; // electron PFO index -> absorbed neutral PFO indices
  std::unordered_set<std::size_t> consumedPfos;

  for (const auto& entry : clusterTrackMaps) { //  Looping over clusterTrackMaps starts for brem recovery

    const auto& cluster = entry.cluster;

    // only an electron (cluster with a track) seeds the brem recovery
    if (entry.tracks.empty() || std::abs(entry.pdg) != 11) {
      continue;
    }

    bool aStatus = false;
    const auto& track1st = entry.tracks[0];
    edm4hep::TrackState trackState1st;
    auto state1st = track1st.getTrackState(edm4hep::TrackState::AtIP);

    if (state1st.has_value()) {
      trackState1st = state1st.value();
    }
    const auto refPt =
        ThreeVector(trackState1st.referencePoint.x, trackState1st.referencePoint.y, trackState1st.referencePoint.z);
    const auto momVec = getMomentum(trackState1st, m_Bz);
    double zPlane = m_eCalEndCapInnerZ;
    double radius = m_eCalBarrelInnerR;
    ThreeVector intersectZ;

    //!!! find brem iteraction in Endacap  region !!!//
    aStatus = GetPointInZ(zPlane, refPt, momVec, intersectZ);

    if (!aStatus) {
      throw std::runtime_error("BremRecoveryAlgorithm: failed to run BremRecoveryUtility::GetPointInZ!");
    }

    ThreeVector finalIntersectPt = intersectZ;

    if (radius <= 0.) { // this should never happen
      throw std::invalid_argument("BremRecoveryAlgorithm: ECAL barrel inner radius < 0!");
    }

    //!!! if brem iteraction in Endacap  region has an x-y coordinate leading a radius greater  !!!//
    //!!! than barrel radius then look for brem interaction in barrel region !!!//

    if (intersectZ.magXY() > radius) { // the tangent intersects the barrel before reaching the endcap - recalculate it

      //!!! handled situation where beam interaction point at x(0),y(0),z(0) !!!//
      if (std::abs(refPt.x) < std::numeric_limits<float>::epsilon() &&
          std::abs(refPt.y) < std::numeric_limits<float>::epsilon() &&
          std::abs(refPt.z) < std::numeric_limits<float>::epsilon()) {

        aStatus = GetPointOnCircleForRefAtOrigin(radius, refPt, momVec, finalIntersectPt);

      } else {
        aStatus = GetPointOnCircle(radius, refPt, momVec, finalIntersectPt);
      }
    }

    if (!aStatus) {
      throw std::runtime_error("BremRecoveryAlgorithm: failed to run BremRecoveryUtility::GetPointOnCircle!");
    }

    double dEtaMargin = m_dEtaMargin;
    double dPhiMargin = m_dPhiMargin;

    const auto& clusterPos1st = cluster.getPosition();
    const auto vecClus1 = ThreeVector(clusterPos1st.x, clusterPos1st.y, clusterPos1st.z);
    const double eta1st = vecClus1.eta();
    const double phi1st = vecClus1.phi();

    //!!! calculation of eta-phi margin for recoverd cluster with input parameters m_dEtaMargin and m_dPhiMargin !!!//
    const double maxEta = eta1st + dEtaMargin;
    const double minEta = eta1st - dEtaMargin;
    const double bremBandPhi = deltaPhi(finalIntersectPt.phi(), phi1st);
    const double maxDPhi = bremBandPhi > 0. ? bremBandPhi + dPhiMargin : dPhiMargin;
    const double minDPhi = bremBandPhi > 0. ? -dPhiMargin : bremBandPhi - dPhiMargin;

    //!!! looping over temporary clusters to check brem counter part of current cluster !!!//
    for (const auto& candidate : tempClusters) { //  Looping over temporary clusters
      const auto& cluster2nd = candidate.cluster;

      const auto& clusterPos2nd = cluster2nd.getPosition();
      const auto vecClus2 = ThreeVector(clusterPos2nd.x, clusterPos2nd.y, clusterPos2nd.z);
      const double eta2nd = vecClus2.eta();
      const double phi2nd = vecClus2.phi();
      const double finalDeltaPhi = deltaPhi(phi2nd, phi1st);

      if (minEta < eta2nd && eta2nd < maxEta && minDPhi < finalDeltaPhi &&
          finalDeltaPhi < maxDPhi) { // start of eta-phi margin inspection
        debug() << "===========found a brem candidate============" << endmsg;

        // a neutral PFO is absorbed by at most one electron (the most energetic one, which comes first)
        if (consumedPfos.insert(candidate.pfoIndex).second) {
          bremPfos[entry.pfoIndex].push_back(candidate.pfoIndex);
        }
      } // end  of eta-phi margin inspection
    } //  end Looping over temporary clusters

  } //  Looping over clusterTrackMaps ends for brem recovery

  //////////////////////////////////////////////////////////////////////////////////
  /////////////////Step 8:complete /////////////////////////////////////////////////
  /////////////////Step 9:complete /////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////
  /////Step 10: output = input PFOs, with the brem clusters added to their electron ////
  /////         and without the absorbed neutral PFOs (no double counting).       ////
  /////         The electron gains their energy but keeps its momentum, charge    ////
  /////         and PDG: brem recovery adds clusters, not tracks.                 ////
  ////////////////////////////////////////////////////////////////////////////////////

  edm4hep::ReconstructedParticleCollection output;
  for (std::size_t pfoIndex = 0; pfoIndex < input.size(); ++pfoIndex) {
    // a neutral PFO identified as brem is now part of its electron: storing it too would double count
    if (consumedPfos.count(pfoIndex))
      continue;

    // every other PFO is copied as is: energy, momentum, charge, PDG, mass and its tracks and clusters
    auto particle = input[pfoIndex].clone();

    // an electron that recovered brem
    if (const auto it = bremPfos.find(pfoIndex); it != bremPfos.end()) {
      for (const auto absorbed : it->second) {
        // takes over the clusters of the absorbed neutral PFO
        for (const auto& cluster : input[absorbed].getClusters()) {
          particle.addToClusters(cluster);
        }
        // and its energy, so that the total energy of the event is unchanged. The momentum is not touched:
        // it comes from the track, which measures the electron before it radiated
        particle.setEnergy(particle.getEnergy() + input[absorbed].getEnergy());
      }
    }
    output.push_back(particle);
  }
  return output;
}
