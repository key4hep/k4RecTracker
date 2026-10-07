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

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"

#include <cmath>
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

  Gaudi::Property<bool> m_useSignificance{
      this, "UseSignificance", true,
      "If true (default), the *Tolerance parameters below are interpreted as significances, i.e. "
      "|diff| / sqrt(sigma(inner)^2 + sigma(outer)^2), using the uncertainties from the track states' covariance "
      "matrices. If false, they are interpreted as absolute differences in the parameter's own units."};

  // Per-parameter matching tolerances. A negative value disables that parameter for matching,
  // i.e. it is not considered when deciding whether two tracks belong together.
  Gaudi::Property<float> m_d0Tolerance{
      this, "D0Tolerance", 3.f,
      "Maximum allowed |D0(inner) - D0(outer)| for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_z0Tolerance{
      this, "Z0Tolerance", 3.f,
      "Maximum allowed |Z0(inner) - Z0(outer)| for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_phiTolerance{
      this, "PhiTolerance", 3.f,
      "Maximum allowed |phi(inner) - phi(outer)| for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_omegaTolerance{this, "OmegaTolerance", 2.f,
                                          "Maximum allowed separation for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_tanLambdaTolerance{
      this, "TanLambdaTolerance", 3.f,
      "Maximum allowed |tanLambda(inner) - tanLambda(outer)| for a match. Negative disables this criterion."};

  StatusCode initialize() override {
    dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
    const double position[3] = {0, 0, 0};
    double magneticFieldVector[3] = {0, 0, 0};
    mainDetector.field().magneticField(position, magneticFieldVector);
    m_Bz = magneticFieldVector[2] / dd4hep::tesla;
    debug() << "B field (T) is : " << m_Bz << endmsg;
    return StatusCode::SUCCESS;
  }

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
  bool isMatch(const edm4hep::Track& t1, edm4hep::TrackState::Location loc1, const edm4hep::Track& t2,
               edm4hep::TrackState::Location loc2) const {
    auto ts1 = t1.getTrackState(loc1);
    auto ts2 = t2.getTrackState(loc2);

    if (!ts1.has_value() || !ts2.has_value()) {
      // It's common for some tracks to lack specific states; verbose instead of debug to avoid spam
      warning() << fmt::format("    [SKIP] Missing requested states (Loc1: {}, Loc2: {})", static_cast<int>(loc1),
                               static_cast<int>(loc2))
                << endmsg;
      return false;
    }

    // Define matching criteria based on differences of the individual track parameters, either as absolute
    // differences or as significances (i.e. normalised by the combined uncertainty of the two states).
    // Parameters whose tolerance is negative are not considered (always pass).
    float d0_value, z0_value, phi_value, omega_value, tanLambda_value;
    if (m_useSignificance) {
      d0_value = significance(ts1->D0, ts2->D0, *ts1, *ts2, TP::d0);
      z0_value = significance(ts1->Z0, ts2->Z0, *ts1, *ts2, TP::z0);
      phi_value = significance(ts1->phi, ts2->phi, *ts1, *ts2, TP::phi);
      tanLambda_value = significance(ts1->tanLambda, ts2->tanLambda, *ts1, *ts2, TP::tanLambda);
      omega_value = ptSignificance(*ts1, *ts2);
    } else {
      d0_value = std::abs(ts1->D0 - ts2->D0);
      z0_value = std::abs(ts1->Z0 - ts2->Z0);
      phi_value = std::abs(ts1->phi - ts2->phi);
      omega_value = std::abs(ts1->omega - ts2->omega);
      tanLambda_value = std::abs(ts1->tanLambda - ts2->tanLambda);
    }

    const bool match = withinTolerance(d0_value, m_d0Tolerance) && withinTolerance(z0_value, m_z0Tolerance) &&
                       withinTolerance(phi_value, m_phiTolerance) && withinTolerance(omega_value, m_omegaTolerance) &&
                       withinTolerance(tanLambda_value, m_tanLambdaTolerance);

    debug() << fmt::format("Comparing Loc {} vs {} ({}): d0={:.4f}, z0={:.4f}, phi={:.4f}, omega={:.4f}, "
                           "tanLambda={:.4f} -> Match: {}",
                           static_cast<int>(loc1), static_cast<int>(loc2),
                           m_useSignificance.value() ? "significance" : "absolute", d0_value, z0_value, phi_value,
                           omega_value, tanLambda_value, match)
            << endmsg;

    return match;
  }

  // A negative tolerance means the corresponding parameter is not considered for matching.
  static bool withinTolerance(float value, float tolerance) { return tolerance < 0.f || value <= tolerance; }

  static float significance(float par1, float par2, const edm4hep::TrackState& ts1, const edm4hep::TrackState& ts2,
                            TP param) {
    const float diff = std::abs(par1 - par2);
    const float sigma = std::sqrt(ts1.getCovMatrix(param, param) + ts2.getCovMatrix(param, param));
    return sigma > 0.f ? diff / sigma : diff;
  }

  // pt derived from the signed curvature omega [1/mm] and the (constant) z-component of the magnetic field [T].
  float ptFromOmega(float omega) const { return 0.3f * m_Bz / (std::abs(omega) * 1000.f); }

  float ptSignificance(const edm4hep::TrackState& ts1, const edm4hep::TrackState& ts2) const {
    const float pt1 = ptFromOmega(ts1.omega);
    const float pt2 = ptFromOmega(ts2.omega);
    const float sigmaPt1 = pt1 * std::sqrt(ts1.getCovMatrix(TP::omega, TP::omega)) / std::abs(ts1.omega);
    const float sigmaPt2 = pt2 * std::sqrt(ts2.getCovMatrix(TP::omega, TP::omega)) / std::abs(ts2.omega);
    const float diff = std::abs(pt1 - pt2);
    const float sigma = std::sqrt(sigmaPt1 * sigmaPt1 + sigmaPt2 * sigmaPt2);
    return sigma > 0.f ? diff / sigma : diff;
  }

  float m_Bz{0.f};
};

DECLARE_COMPONENT(TrackMerger)
