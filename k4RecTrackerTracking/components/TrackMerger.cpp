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

#include "edm4hep/Track.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/TrackerHit.h"
#include "k4FWCore/Transformer.h"

#include "edm4hep/MCParticleCollection.h"

#include <string>

// Type aliases for improved readability
using TrackColl = edm4hep::TrackCollection;
using Track = edm4hep::Track;
using TP = edm4hep::TrackParams;
using TS = edm4hep::TrackState;

struct TrackMerger final : k4FWCore::Transformer<TrackColl(const TrackColl&, const TrackColl&)> {
  TrackMerger(const std::string& name, ISvcLocator* svcLoc)
      : Transformer(name, svcLoc,
                    {
                        KeyValue("InputInnerTracks", {"InnerTracks"}),
                        KeyValue("InputOuterTracks", {"OuterTracks"}),
                    },
                    {KeyValue("OutTracks", {"MyCandidateMergedTracks"})}) {}

  Gaudi::Property<bool> m_greedy{this, "Greedy", true, "If true, each track is used only once."};

  // Per-parameter matching tolerances. A negative value disables that parameter for matching,
  // i.e. it is not considered when deciding whether two tracks belong together.
  // Defaults reproduce the original matching criterion: only D0 and Z0 are considered.
  Gaudi::Property<float> m_d0Tolerance{
      this, "D0Tolerance", 0.5f,
      "Maximum allowed |D0(inner) - D0(outer)| for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_z0Tolerance{
      this, "Z0Tolerance", 2.5f,
      "Maximum allowed |Z0(inner) - Z0(outer)| for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_phiTolerance{
      this, "PhiTolerance", -1.f,
      "Maximum allowed |phi(inner) - phi(outer)| for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_omegaTolerance{
      this, "OmegaTolerance", -1.f,
      "Maximum allowed |omega(inner) - omega(outer)| for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_tanLambdaTolerance{
      this, "TanLambdaTolerance", -1.f,
      "Maximum allowed |tanLambda(inner) - tanLambda(outer)| for a match. Negative disables this criterion."};

  TrackColl operator()(const TrackColl& inputInnerTracks, const TrackColl& inputOuterTracks) const override {
    auto outTracks = TrackColl();

    debug() << "Received InnerTracks collection with " << inputInnerTracks.size() << " tracks" << endmsg;
    debug() << "Received OuterTracks collection with " << inputOuterTracks.size() << " tracks" << endmsg;

    // 1. Check if both input collections have at least one entry
    if (inputInnerTracks.empty() || inputOuterTracks.empty()) {
      warning() << "One of the input collections is empty. InnerTracks: " << inputInnerTracks.size()
                << ", OuterTracks: " << inputOuterTracks.size() << ". Skipping track merging for this event." << endmsg;
      return outTracks; // Returns empty collection
    }

    // Flag to ensure each outer track is only merged once
    std::vector<bool> usedOuterTracks(inputOuterTracks.size(), false);

    // Loop over inner tracks with index for debug output
    for (size_t iInner = 0; iInner < inputInnerTracks.size(); iInner++) {
      const auto trackInner = inputInnerTracks[iInner];
      bool matched = false;

      // Explicit index loop for outer tracks to manage 'usedOuterTracks' and debug indexing
      // LOGIC NOTE: This algorithm accepts the FIRST match found within tolerances.
      // It does not perform a global chi2 minimization or search for the "best" match.
      for (size_t iOuter = 0; iOuter < inputOuterTracks.size(); iOuter++) {
        if (m_greedy && usedOuterTracks[iOuter])
          continue;

        const auto trackOuter = inputOuterTracks[iOuter];

        // compare trackInner and trackOuter at their respective track states (Inner: last hit, Outer: first hit)
        // to determine if they likely originate from the same particle
        if (isMatch(trackInner, TS::AtLastHit, trackOuter, TS::AtFirstHit)) {
          debug() << fmt::format("  [MATCH] Inner track {} matched with Outer track {}. Creating merged track.", iInner,
                                 iOuter)
                  << endmsg;

          auto newTrack = outTracks.create();

          // Combine hits from both tracks
          for (const auto& hit : trackInner.getTrackerHits())
            newTrack.addToTrackerHits(hit);
          for (const auto& hit : trackOuter.getTrackerHits())
            newTrack.addToTrackerHits(hit);

          // Maintain navigation/provenance by linking parent tracks
          newTrack.addToTracks(trackInner);
          newTrack.addToTracks(trackOuter);

          matched = true;
          if (m_greedy) {
            usedOuterTracks[iOuter] = true;
            break; // Exit inner loop: current inner track is satisfied
          }
        }
      }

      if (!matched) {
        debug() << fmt::format("  [INFO] Inner track {} found no matching outer track within tolerances.", iInner)
                << endmsg;
      }
    }

    debug() << fmt::format(
                   "Event processing complete. Created {} merged tracks from {} InnerTracks and {} OuterTracks.",
                   outTracks.size(), inputInnerTracks.size(), inputOuterTracks.size())
            << endmsg;
    return outTracks;
  }

private:
  bool isMatch(const edm4hep::Track& t1, int loc1, const edm4hep::Track& t2, int loc2) const {
    auto ts1 = getTrackState(t1, loc1);
    auto ts2 = getTrackState(t2, loc2);

    if (!ts1.has_value() || !ts2.has_value()) {
      // It's common for some tracks to lack specific states; verbose instead of debug to avoid spam
      warning() << fmt::format("    [SKIP] Missing requested states (Loc1: {}, Loc2: {})", loc1, loc2) << endmsg;
      return false;
    }

    // Define matching criteria based on differences of the individual track parameters.
    // Parameters whose tolerance is negative are not considered (always pass).
    const float d0_diff = std::abs(ts1->D0 - ts2->D0);
    const float z0_diff = std::abs(ts1->Z0 - ts2->Z0);
    const float phi_diff = std::abs(ts1->phi - ts2->phi);
    const float omega_diff = std::abs(ts1->omega - ts2->omega);
    const float tanLambda_diff = std::abs(ts1->tanLambda - ts2->tanLambda);

    const bool match = withinTolerance(d0_diff, m_d0Tolerance) && withinTolerance(z0_diff, m_z0Tolerance) &&
                       withinTolerance(phi_diff, m_phiTolerance) && withinTolerance(omega_diff, m_omegaTolerance) &&
                       withinTolerance(tanLambda_diff, m_tanLambdaTolerance);

    debug() << fmt::format("    Comparing Loc {} vs {}: d0_diff={:.4f}, z0_diff={:.4f}, phi_diff={:.4f}, "
                           "omega_diff={:.4f}, tanLambda_diff={:.4f} -> Match: {}",
                           loc1, loc2, d0_diff, z0_diff, phi_diff, omega_diff, tanLambda_diff, match)
            << endmsg;

    return match;
  }

  // A negative tolerance means the corresponding parameter is not considered for matching.
  static bool withinTolerance(float diff, float tolerance) { return tolerance < 0.f || diff <= tolerance; }

  std::optional<TS> getTrackState(Track track, const int loc) const {
    auto ts = track.getTrackState(loc);
    if (!ts.has_value()) {
      warning() << std::format("No track state at location {} found!", loc) << endmsg;
    }
    return ts;
  }
};

DECLARE_COMPONENT(TrackMerger)
