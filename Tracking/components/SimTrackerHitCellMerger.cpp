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

// Gaudi
#include "Gaudi/Property.h"

// edm4hep
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

// k4FWCore
#include "k4FWCore/Transformer.h"

// C++
#include <algorithm>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

/** @class SimTrackerHitCellMerger
 *
 *  Gaudi transformer that accumulates the Geant4 step lengths (edm4hep::SimTrackerHit::pathLength) of all
 *  simulated hits sharing the same cellID, and writes the result as a new, "merged"
 *  edm4hep::SimTrackerHitCollection with one entry per cell (or per cell and track, see below).
 *
 *  A single Geant4 track usually leaves several SimTrackerHits in one cell (one per step), and several
 *  tracks (the primary plus its delta rays, conversions, ...) can cross the very same cell. What to do
 *  with those different tracks is steered by the MultipleTrackHandling property:
 *
 *   - "All"         : one output hit per cellID. pathLength and eDep are summed over *every* contributing
 *                     hit, irrespective of which track produced it. The MCParticle relation of the output
 *                     hit points to the contributor with the lowest Geant4 track number (see below), i.e.
 *                     the most primary track crossing the cell.
 *   - "PrimaryOnly" : one output hit per cellID, but only the hits of the contributor with the lowest
 *                     Geant4 track number are summed. Contributions of later (secondary) tracks are
 *                     dropped, which is what one wants when the accumulated path length is meant to
 *                     describe the primary particle traversing the cell.
 *   - "PerTrack"    : one output hit per (cellID, track) pair, each summing only the hits of that track.
 *                     The output data type is unchanged; the only difference is that the same cellID may
 *                     now appear several times in the output collection, once per contributing track,
 *                     each entry carrying its own MCParticle relation. Entries of a given cell are
 *                     written in order of increasing Geant4 track number, so the first entry for a cell
 *                     is the most primary one.
 *
 *  Geant4 track number: edm4hep does not store the Geant4 trackID. However, DD4hep/DDG4 fills the
 *  MCParticle collection ordered by increasing trackID, so the index of the MCParticle within its
 *  collection is a monotonic proxy for it: primaries come first, secondaries later. This algorithm
 *  therefore ranks contributors by edm4hep::MCParticle::getObjectID().index. Note that hits flagged
 *  isProducedBySecondary() do not point to the secondary that actually created them (it was not kept in
 *  the MCParticle collection) but to its surviving ancestor, so their step length is booked on that
 *  ancestor. Set ExcludeSecondaryHits to true to skip such hits altogether.
 *
 *  The remaining fields of a merged hit are taken from its *representative* hit, defined as the earliest
 *  (smallest time) of the hits that were summed into it: time, position, momentum and quality. Only
 *  pathLength and eDep are accumulated.
 *
 *  @author Andreas Loeschcke Centeno
 */

struct SimTrackerHitCellMerger final
    : k4FWCore::Transformer<edm4hep::SimTrackerHitCollection(const edm4hep::SimTrackerHitCollection&)> {

  SimTrackerHitCellMerger(const std::string& name, ISvcLocator* svcLoc)
      : Transformer(name, svcLoc, {KeyValues("InputSimTrackerHits", {"SimTrackerHits"})},
                    {KeyValues("OutputSimTrackerHits", {"MergedSimTrackerHits"})}) {}

  /// How to treat several Geant4 tracks crossing the same cell
  enum class TrackHandling { All, PrimaryOnly, PerTrack };

  StatusCode initialize() override {
    if (m_multipleTrackHandling == "All") {
      m_trackHandling = TrackHandling::All;
    } else if (m_multipleTrackHandling == "PrimaryOnly") {
      m_trackHandling = TrackHandling::PrimaryOnly;
    } else if (m_multipleTrackHandling == "PerTrack") {
      m_trackHandling = TrackHandling::PerTrack;
    } else {
      error() << "Unknown MultipleTrackHandling '" << m_multipleTrackHandling.value()
              << "'. Valid values are 'All', 'PrimaryOnly' and 'PerTrack'." << endmsg;
      return StatusCode::FAILURE;
    }
    info() << "Accumulating SimTrackerHit path lengths per cellID with MultipleTrackHandling = "
           << m_multipleTrackHandling.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  edm4hep::SimTrackerHitCollection operator()(const edm4hep::SimTrackerHitCollection& simTrackerHits) const override {

    auto output = edm4hep::SimTrackerHitCollection();

    debug() << "Received SimTrackerHit collection with " << simTrackerHits.size() << " hits" << endmsg;

    // Accumulate per cell and, within a cell, per contributing track. The inner map is ordered so that
    // iterating it yields the contributors sorted by increasing Geant4 track number, i.e. begin() is the
    // most primary track crossing the cell.
    std::unordered_map<std::uint64_t, std::map<TrackKey, Contribution>> cellMap;

    for (const auto& hit : simTrackerHits) {

      if (m_excludeSecondaryHits && hit.isProducedBySecondary())
        continue;

      const auto particle = hit.getParticle();
      if (!particle.isAvailable()) {
        warning() << "SimTrackerHit in cell " << hit.getCellID()
                  << " has no MCParticle relation, its path length cannot be attributed to a track. Skipping it."
                  << endmsg;
        continue;
      }

      const auto objectID = particle.getObjectID();
      auto& contribution = cellMap[hit.getCellID()][TrackKey{objectID.index, objectID.collectionID}];

      contribution.pathLength += hit.getPathLength();
      contribution.eDep += hit.getEDep();
      contribution.nHits++;
      contribution.particle = particle;
      // Keep the earliest hit of this track in this cell as the representative one
      if (!contribution.representative.isAvailable() || hit.getTime() < contribution.representative.getTime())
        contribution.representative = hit;
    }

    // Write the output in order of increasing cellID so that the result does not depend on the hash
    // ordering of the map above
    std::vector<std::uint64_t> cellIDs;
    cellIDs.reserve(cellMap.size());
    for (const auto& [cellID, contributions] : cellMap)
      cellIDs.push_back(cellID);
    std::sort(cellIDs.begin(), cellIDs.end());

    for (const auto cellID : cellIDs) {
      const auto& contributions = cellMap.at(cellID);

      if (m_trackHandling == TrackHandling::PerTrack) {
        // One output hit per (cell, track), ordered by increasing Geant4 track number
        for (const auto& [trackKey, contribution] : contributions)
          addMergedHit(output, cellID, contribution.pathLength, contribution.eDep, contribution.particle,
                       contribution.representative, contribution.nHits);
        continue;
      }

      // The contributor with the lowest Geant4 track number, i.e. the most primary one
      const auto& primaryContribution = contributions.begin()->second;

      if (m_trackHandling == TrackHandling::PrimaryOnly) {
        // Only the primary track's steps are summed, everything else in the cell is dropped
        addMergedHit(output, cellID, primaryContribution.pathLength, primaryContribution.eDep,
                     primaryContribution.particle, primaryContribution.representative, primaryContribution.nHits);
        continue;
      }

      // TrackHandling::All: sum over all tracks, but attribute the merged hit to the most primary one and
      // use the globally earliest hit of the cell as representative
      double pathLength = 0.;
      double eDep = 0.;
      int nHits = 0;
      edm4hep::SimTrackerHit representative{};
      for (const auto& [trackKey, contribution] : contributions) {
        pathLength += contribution.pathLength;
        eDep += contribution.eDep;
        nHits += contribution.nHits;
        if (!representative.isAvailable() || contribution.representative.getTime() < representative.getTime())
          representative = contribution.representative;
      }
      addMergedHit(output, cellID, pathLength, eDep, primaryContribution.particle, representative, nHits);
    }

    debug() << "Wrote " << output.size() << " merged SimTrackerHits for " << cellMap.size() << " cells" << endmsg;

    return output;
  }

private:
  /// Ranking key of a contributing track: the index of its MCParticle within its collection, used as a
  /// proxy for the Geant4 track number, with the collectionID as tie breaker for the (pathological) case
  /// of hits pointing into more than one MCParticle collection.
  using TrackKey = std::pair<int, std::uint32_t>;

  /// Accumulated step lengths and energy deposits of one track in one cell
  struct Contribution {
    double pathLength = 0.;
    double eDep = 0.;
    int nHits = 0;
    edm4hep::MCParticle particle{};
    /// Earliest of the summed hits, provides all fields that are not accumulated
    edm4hep::SimTrackerHit representative{};
  };

  /// Create one merged hit from an accumulated path length and its representative hit
  void addMergedHit(edm4hep::SimTrackerHitCollection& output, std::uint64_t cellID, double pathLength, double eDep,
                    const edm4hep::MCParticle& particle, const edm4hep::SimTrackerHit& representative,
                    int nHits) const {
    if (pathLength < m_minPathLength) {
      debug() << "Dropping merged hit in cell " << cellID << " with path length " << pathLength
              << " mm below MinPathLength" << endmsg;
      return;
    }
    auto mergedHit = output.create();
    mergedHit.setCellID(cellID);
    mergedHit.setPathLength(pathLength);
    mergedHit.setEDep(eDep);
    mergedHit.setParticle(particle);
    mergedHit.setTime(representative.getTime());
    mergedHit.setPosition(representative.getPosition());
    mergedHit.setMomentum(representative.getMomentum());
    mergedHit.setQuality(representative.getQuality());
    verbose() << "Cell " << cellID << ": merged " << nHits << " hits of MCParticle " << particle.getObjectID().index
              << " (PDG " << particle.getPDG() << ") into path length " << pathLength << " mm, eDep " << eDep << " GeV"
              << endmsg;
  }

  /// Configurable property steering what to do when several Geant4 tracks cross the same cell
  Gaudi::Property<std::string> m_multipleTrackHandling{
      this, "MultipleTrackHandling", "All",
      "How to treat several Geant4 tracks in the same cell: 'All' sums the step lengths of every track, "
      "'PrimaryOnly' sums only those of the track with the lowest Geant4 track number, 'PerTrack' writes one "
      "output hit per cell and track"};

  /// Configurable property to skip hits that were created by a secondary which is not kept in the
  /// MCParticle collection, and whose step length would otherwise be booked on its surviving ancestor
  Gaudi::Property<bool> m_excludeSecondaryHits{
      this, "ExcludeSecondaryHits", false,
      "Skip SimTrackerHits flagged isProducedBySecondary() instead of attributing them to the ancestor "
      "MCParticle they point to"};

  /// Configurable property to suppress merged hits with a negligible accumulated path length
  Gaudi::Property<float> m_minPathLength{
      this, "MinPathLength", 0., "Do not write merged hits whose accumulated path length (in mm) is below this value"};

  /// Parsed version of m_multipleTrackHandling, set in initialize()
  TrackHandling m_trackHandling{TrackHandling::All};
};

DECLARE_COMPONENT(SimTrackerHitCellMerger)
