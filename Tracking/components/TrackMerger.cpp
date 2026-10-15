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

#include "Gaudi/Accumulators.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/TrackerHit.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"

#include "edm4hep/MCParticleCollection.h"

#include "DD4hep/Detector.h"

#include <cmath>
#include <string>
#include <utility>

#include <utils.h>

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
      "If true (default), match on the *Significance cuts, i.e. |diff| / sqrt(var(inner) + var(outer)) using the "
      "track states' covariance matrices. If false, match on the *Tolerance cuts (absolute differences)."};

  // Per-parameter matching cuts. Only the set selected by UseSignificance is used. A negative value disables that
  // parameter for matching, i.e. it is not considered when deciding whether two tracks belong together.
  // Absolute cuts (UseSignificance = false):
  Gaudi::Property<float> m_d0Tolerance{
      this, "D0Tolerance", 0.5f,
      "Absolute mode: maximum allowed |D0(inner) - D0(outer)| [mm] for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_z0Tolerance{
      this, "Z0Tolerance", 2.5f,
      "Absolute mode: maximum allowed |Z0(inner) - Z0(outer)| [mm] for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_phiTolerance{this, "PhiTolerance", -1.f,
                                        "Absolute mode: maximum allowed |phi(inner) - phi(outer)| [rad], wrapped into "
                                        "[-pi, pi], for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_omegaTolerance{this, "OmegaTolerance", -1.f,
                                          "Absolute mode: maximum allowed |omega(inner) - omega(outer)| [1/mm] for a "
                                          "match. Negative disables this criterion."};
  Gaudi::Property<float> m_tanLambdaTolerance{this, "TanLambdaTolerance", -1.f,
                                              "Absolute mode: maximum allowed |tanLambda(inner) - tanLambda(outer)| "
                                              "for a match. Negative disables this criterion."};
  // Significance cuts (UseSignificance = true):
  Gaudi::Property<float> m_d0Significance{
      this, "D0Significance", 3.f,
      "Significance mode: maximum allowed D0 significance for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_z0Significance{
      this, "Z0Significance", 3.f,
      "Significance mode: maximum allowed Z0 significance for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_phiSignificance{this, "PhiSignificance", 3.f,
                                           "Significance mode: maximum allowed phi significance (difference wrapped "
                                           "into [-pi, pi]) for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_ptSignificance{
      this, "PtSignificance", 2.f,
      "Significance mode: maximum allowed significance of the signed pT (charge sign times pT, derived from omega "
      "and the magnetic field) for a match. Requires a geometry with a field. Negative disables this criterion."};
  Gaudi::Property<float> m_tanLambdaSignificance{
      this, "TanLambdaSignificance", 3.f,
      "Significance mode: maximum allowed tanLambda significance for a match. Negative disables this criterion."};

  StatusCode initialize() override {
    // The magnetic field is only needed to convert omega to pT for the pT significance.
    if (!m_useSignificance || m_ptSignificance < 0.f) {
      return StatusCode::SUCCESS;
    }
    if (!m_geoSvc.retrieve()) {
      error() << "Unable to retrieve the GeoSvc, which is needed for the magnetic field" << endmsg;
      return StatusCode::FAILURE;
    }
    m_Bz = getBzAtOrigin(*m_geoSvc->getDetector());
    if (m_Bz == 0.f) {
      error() << "Magnetic field is zero, so pT cannot be derived from omega. Load a geometry with a field, or "
                 "disable the pT criterion (PtSignificance < 0 or UseSignificance = false)"
              << endmsg;
      return StatusCode::FAILURE;
    }
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
    const float d0Cut = m_useSignificance ? m_d0Significance : m_d0Tolerance;
    const float z0Cut = m_useSignificance ? m_z0Significance : m_z0Tolerance;
    const float phiCut = m_useSignificance ? m_phiSignificance : m_phiTolerance;
    const float omegaCut = m_useSignificance ? m_ptSignificance : m_omegaTolerance;
    const float tanLambdaCut = m_useSignificance ? m_tanLambdaSignificance : m_tanLambdaTolerance;
    if (m_useSignificance) {
      // Only check parameters whose criterion is enabled. The warning is rate-limited by the counter, which also
      // reports the total number of occurrences at finalize.
      for (const auto& [param, tolerance] :
           {std::pair{TP::d0, d0Cut}, std::pair{TP::z0, z0Cut}, std::pair{TP::phi, phiCut},
            std::pair{TP::omega, omegaCut}, std::pair{TP::tanLambda, tanLambdaCut}}) {
        if (tolerance >= 0.f && ts1->getCovMatrix(param, param) + ts2->getCovMatrix(param, param) <= 0.f) {
          ++m_zeroUncertaintyWarning;
          debug() << "Combined uncertainty is zero for track state parameter " << static_cast<int>(param) << endmsg;
        }
      }
      d0_value = significance(ts1->D0, ts2->D0, ts1->getCovMatrix(TP::d0, TP::d0), ts2->getCovMatrix(TP::d0, TP::d0));
      z0_value = significance(ts1->Z0, ts2->Z0, ts1->getCovMatrix(TP::z0, TP::z0), ts2->getCovMatrix(TP::z0, TP::z0));
      // Pass the wrapped difference against 0 so that states on either side of phi = +-pi are compared correctly.
      phi_value = significance(deltaPhi(ts1->phi, ts2->phi), 0.f, ts1->getCovMatrix(TP::phi, TP::phi),
                               ts2->getCovMatrix(TP::phi, TP::phi));
      tanLambda_value = significance(ts1->tanLambda, ts2->tanLambda, ts1->getCovMatrix(TP::tanLambda, TP::tanLambda),
                                     ts2->getCovMatrix(TP::tanLambda, TP::tanLambda));
      omega_value = ptSignificance(*ts1, *ts2);
    } else {
      d0_value = std::abs(ts1->D0 - ts2->D0);
      z0_value = std::abs(ts1->Z0 - ts2->Z0);
      phi_value = std::abs(deltaPhi(ts1->phi, ts2->phi));
      omega_value = std::abs(ts1->omega - ts2->omega);
      tanLambda_value = std::abs(ts1->tanLambda - ts2->tanLambda);
    }

    const bool match = withinTolerance(d0_value, d0Cut) && withinTolerance(z0_value, z0Cut) &&
                       withinTolerance(phi_value, phiCut) && withinTolerance(omega_value, omegaCut) &&
                       withinTolerance(tanLambda_value, tanLambdaCut);

    debug() << fmt::format("Comparing Loc {} vs {} ({}): d0={:.4f}, z0={:.4f}, phi={:.4f}, omega={:.4f}, "
                           "tanLambda={:.4f} -> Match: {}",
                           static_cast<int>(loc1), static_cast<int>(loc2),
                           m_useSignificance ? "significance" : "absolute", d0_value, z0_value, phi_value, omega_value,
                           tanLambda_value, match)
            << endmsg;

    return match;
  }

  // Difference phi1 - phi2 wrapped into [-pi, pi].
  static float deltaPhi(float phi1, float phi2) { return std::remainder(phi1 - phi2, 2.f * static_cast<float>(M_PI)); }

  // A negative tolerance means the corresponding parameter is not considered for matching.
  static bool withinTolerance(float value, float tolerance) { return tolerance < 0.f || value <= tolerance; }

  // Difference of two values normalised by the combined uncertainty. Falls back to the plain difference if the
  // combined uncertainty vanishes; the caller warns about that case.
  static float significance(float par1, float par2, float var1, float var2) {
    const float diff = std::abs(par1 - par2);
    const float sigma = std::sqrt(var1 + var2);
    return sigma > 0.f ? diff / sigma : diff;
  }

  // Signed pt (charge sign times pt) derived from the signed curvature omega [1/mm] and the (constant) z-component of
  // the magnetic field [T]. Keeping the sign makes segments with opposite charge incompatible.
  float signedPtFromOmega(float omega) const { return static_cast<float>(ConversionUnits::a_lcio) * m_Bz / omega; }

  float ptSignificance(const edm4hep::TrackState& ts1, const edm4hep::TrackState& ts2) const {
    if (ts1.omega == 0.f || ts2.omega == 0.f) {
      debug() << "Track state with omega == 0; using the omega significance instead of the pT significance" << endmsg;
      return significance(ts1.omega, ts2.omega, ts1.getCovMatrix(TP::omega, TP::omega),
                          ts2.getCovMatrix(TP::omega, TP::omega));
    }
    const float pt1 = signedPtFromOmega(ts1.omega);
    const float pt2 = signedPtFromOmega(ts2.omega);
    // Error propagation: var(pt) = (pt / omega)^2 * var(omega)
    const float varPt1 = pt1 * pt1 * ts1.getCovMatrix(TP::omega, TP::omega) / (ts1.omega * ts1.omega);
    const float varPt2 = pt2 * pt2 * ts2.getCovMatrix(TP::omega, TP::omega) / (ts2.omega * ts2.omega);
    return significance(pt1, pt2, varPt1, varPt2);
  }

  mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING> m_zeroUncertaintyWarning{
      this, "Combined uncertainty of a track parameter is zero; using the absolute difference instead", 10};
  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Detector geometry service"};
  float m_Bz{0.f};
};

DECLARE_COMPONENT(TrackMerger)
