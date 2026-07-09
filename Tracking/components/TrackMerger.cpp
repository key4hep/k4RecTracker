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
          debug() << fmt::format("  [MATCH] Inner track {} matched with Outer track {}. Creating merged track.",
                                 iInner, iOuter)
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

    // Define matching criteria based on d0 and z0 differences
    const float d0_diff = std::abs(ts1->D0 - ts2->D0);
    const float z0_diff = std::abs(ts1->Z0 - ts2->Z0);
    const bool match = (d0_diff <= 0.5f && z0_diff <= 2.5f);

    debug() << fmt::format("    Comparing Loc {} vs {}: d0_diff={:.4f}, z0_diff={:.4f} -> Match: {}", loc1, loc2,
                           d0_diff, z0_diff, match)
            << endmsg;

    return match;
  }

  std::optional<const TS> getTrackState(const Track& track, const int loc) const {
    const auto& trackStates = track.getTrackStates();

    if (auto it = std::ranges::find(trackStates, loc, &TS::location); it != trackStates.end()) {
      return *it;
    } else {
      warning() << std::format("No track state at location {} found!", loc) << endmsg;
      return std::nullopt;
    }
  }
};

DECLARE_COMPONENT(TrackMerger)
