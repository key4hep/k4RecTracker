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

//#include "MarlinTrk/Factory.h"
//#include "MarlinTrk/HelixTrack.h"
//#include "MarlinTrk/IMarlinTrack.h"
//#include "MarlinTrk/IMarlinTrkSystem.h"
//#include "MarlinTrk/MarlinTrkUtils.h"

#include "k4Reco/GaudiDDKalTest.h"
#include "k4Reco/GaudiDDKalTestTrack.h"
//#include "k4Reco/GaudiTrkUtils.h"
#include <DDSegmentation/BitFieldCoder.h>

#include "IMPL/TrackImpl.h"
#include "IMPL/TrackerHitImpl.h"

#include <cmath>
#include <string>
#include <utility>

#include <utils.h>

// Type aliases for improved readability
using TrackColl = edm4hep::TrackCollection;
using Track = edm4hep::Track;
using TP = edm4hep::TrackParams;
using TS = edm4hep::TrackState;
using TH = edm4hep::TrackerHit;

struct TrackMerger final : k4FWCore::Transformer<TrackColl(const TrackColl&, const TrackColl&)> {
  TrackMerger(const std::string& name, ISvcLocator* svcLoc)
      : Transformer(name, svcLoc,
                    {
                        KeyValue("InputInnerTracks", {"InnerTracks"}),
                        KeyValue("InputOuterTracks", {"OuterTracks"}),
                    },
                    {KeyValue("OutTracks", {"MyCandidateMergedTracks"})}) {}

  Gaudi::Property<bool> m_greedy{this, "Greedy", true, "If true, each track is used only once."};

  Gaudi::Property<bool> m_fitOutwards{this, "fitOutwards", false, "If true fit inside out; if false fit outside in."};

  Gaudi::Property<float> m_maxChi2PerHit{this, "maxChi2PerHit", 100.0f,
      "Maximum additional chi2 when a hit is added to the fit when testing a combined track."};

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
  Gaudi::Property<float> m_d0Significance{this, "D0Significance", 3.f,
      "Significance mode: maximum allowed D0 significance for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_z0Significance{this, "Z0Significance", 3.f,
      "Significance mode: maximum allowed Z0 significance for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_phiSignificance{this, "PhiSignificance", 3.f,
      "Significance mode: maximum allowed phi significance (difference wrapped "
      "into [-pi, pi]) for a match. Negative disables this criterion."};
  Gaudi::Property<float> m_ptSignificance{this, "PtSignificance", 2.f,
      "Significance mode: maximum allowed significance of the signed pT (charge sign times pT, derived from omega "
      "and the magnetic field) for a match. Requires a geometry with a field. Negative disables this criterion."};
  Gaudi::Property<float> m_tanLambdaSignificance{this, "TanLambdaSignificance", 3.f,
      "Significance mode: maximum allowed tanLambda significance for a match. Negative disables this criterion."};

  Gaudi::Property<std::string> m_trkSystemName{this, "TrkSystemName", "DDKalTest", "Specify the tracking system to be used; default: DDKalTest."};
  Gaudi::Property<bool> m_MSOn{this, "MultipleScatteringOn", true, "Use MultipleScattering in Fit."};
  Gaudi::Property<bool> m_ElossOn{this, "EnergyLossOn", true, "Use Energy Loss in Fit."};
  Gaudi::Property<bool> m_SmoothOn{this, "SmoothOn", true, "Smooth All Mesurement Sites in Fit."};

  Gaudi::Property<std::string> m_encodingStringVariable{
    this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
    "The name of the DD4hep constant that contains the Encoding string for tracking detectors"};
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};



  StatusCode initialize() override {

    // The magnetic field is only needed to convert omega to pT for the pT significance.
    //if (!m_useSignificance || m_ptSignificance < 0.f) {
    //  return StatusCode::SUCCESS;
    //}
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

    std::string cellIDEncodingString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
    m_encoder = dd4hep::DDSegmentation::BitFieldCoder(cellIDEncodingString);
    debug() << cellIDEncodingString << std::endl;

    m_ddkaltest.init();
    m_ddkaltest.setEncoder(m_encoder);

    return StatusCode::SUCCESS;
  }


  TrackColl operator()(const TrackColl& inputInnerTracks, const TrackColl& inputOuterTracks) const override {
    auto outTracks = TrackColl();

    debug() << "Received InnerTracks collection with " << inputInnerTracks.size() << " tracks" << endmsg;
    debug() << "Received OuterTracks collection with " << inputOuterTracks.size() << " tracks" << endmsg;

    // std::string cellIDEncodingString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
    // debug() << "m_encodingStringVariable: " << m_encodingStringVariable << std::endl;
    // debug() << "cellIDEncodingString: " << cellIDEncodingString << std::endl;

    // Check if both input collections are empty
    if (inputInnerTracks.empty() && inputOuterTracks.empty()) {
      warning() << "Input collections are empty. Skipping track merging for this event and returning empty track collection." << endmsg;
      return outTracks;
    }

    // note which tracks are still unmerged
    std::vector<Track> unmergedInnerTracks, unmergedOuterTracks;
    for (const Track& innerTrack : inputInnerTracks)
    {unmergedInnerTracks.push_back(innerTrack);}
    for (const Track& outerTrack : inputOuterTracks)
    {unmergedOuterTracks.push_back(outerTrack);}

    if (checkForNoUnmergedTracks(unmergedInnerTracks, unmergedOuterTracks, outTracks))
    {return outTracks;}

    // check combinations: if their parameters match save the combination and the combined fit chi2/ndf
    std::map<float, std::pair<std::pair<Track, Track>, std::vector<TS>>> chi2OfCombination;

    for (const Track& innerTrack : inputInnerTracks){
      for (const Track& outerTrack : inputOuterTracks){
        if (isMatch(innerTrack, TS::AtLastHit, outerTrack, TS::AtFirstHit)){

          auto marlin_trk = GaudiDDKalTestTrack(this, const_cast<GaudiDDKalTest*>(&m_ddkaltest));
          info() << "GaudiDDKalTestTrack instance loaded" << std::endl;

          //auto combinedTrack = outTracks.create();
          // Combine hits from both tracks
          for (const TH& hit : innerTrack.getTrackerHits()){
            //combinedTrack.addToTrackerHits(hit);

            marlin_trk.addHit(&hit);
          }
          for (const TH& hit : outerTrack.getTrackerHits()){
            //combinedTrack.addToTrackerHits(hit);
            info() << "adding hit " << &hit << "  cellID: " << hit.getCellID() << "  position: " << hit.getPosition()[0] << " " << hit.getPosition()[1] << " " << hit.getPosition()[2] << std::endl;
            marlin_trk.addHit(&hit);
          }
          debug() << "hits added" << std::endl;

          // get track state to initialise MarlinTrk fitter and fit combined track
          auto optState = outerTrack.getTrackState(TS::AtLastHit);
          if (!optState)
          {
            error() << "TrackState is not available" << std::endl;
            continue;
          }
          const auto outerState = *optState;
          debug() << "track state received" << std::endl;

          int ini_status = marlin_trk.initialise(outerState, m_fitOutwards);
          debug() << "ini status: " << ini_status << std::endl;
          int fit_status = marlin_trk.fit(m_maxChi2PerHit);
          debug() << "fit status: " << fit_status << std::endl;

          // get chi2 from the inner most used hit and store track combination with chi2/ndf
          int ndf = 0;
          double chi2 = -DBL_MAX;
          TS tsAtFirstHit;
          bool combined = false;

          for (const TH& hit : innerTrack.getTrackerHits()){
            info() << "inner hit for track state: " << &hit << "  cellID: " << hit.getCellID() << endmsg;
            int state_status = marlin_trk.getTrackState(&hit, tsAtFirstHit , chi2 , ndf);
            info() << " track state returns status code " << state_status << "; NDF = " << ndf << "; chi2 = " << chi2 << endmsg;

            if (state_status == 0){
              // get track state at the IP
              TS tsAtIP, tsAtLastHit;
              const edm4hep::Vector3d ipv(0., 0., 0.);
              marlin_trk.propagate(ipv, &hit, tsAtIP, chi2, ndf);
              tsAtIP.location = TS::AtIP;

              // get track state at the last hit
              for (const TH& ohit : outerTrack.getTrackerHits() | std::views::reverse){
                int state_status2 = marlin_trk.getTrackState(&ohit, tsAtLastHit , chi2 , ndf);
                if (state_status == 0)
                  break;
              }
              tsAtFirstHit.location = TS::AtFirstHit;
              tsAtLastHit.location = TS::AtLastHit;

              std::vector<TS> trackStates{tsAtIP, tsAtFirstHit, tsAtLastHit};
              chi2OfCombination[chi2/ndf] = std::make_pair(std::make_pair(innerTrack, outerTrack), trackStates);
              combined = true;
              break;
            }
          }
          if (!combined){
            warning() << "track combination rejected for all inner hits!" << endmsg;
            printTrack(innerTrack, TS::AtLastHit);
            printTrack(outerTrack, TS::AtFirstHit);
          }
        } // end isMatch, end for loops over tracks
      }
    }

    info() << "Unmerged tracks after matching and fitting:" << endmsg;
    for (Track& track: unmergedInnerTracks)
      printTrack(track,TS::AtLastHit);
    for (Track& track: unmergedOuterTracks)
      printTrack(track,TS::AtFirstHit);


    //iterate over track combinations starting from lowest chi2/ndf
    //if tracks are unmerged accept track combination and remove tracks from unmerged tracks vectors
    for (const auto& [chi2, combination] : chi2OfCombination){
      const Track& innerTrack = combination.first.first;
      const Track& outerTrack = combination.first.second;
      // int i=0, k=0;
      // bool removeTracks = false;
      // for (i = 0; i < unmergedInnerTracks.size(); ++i){
      //   if (unmergedInnerTracks[i] == innerTrack){
      //     for (i = 0; i < unmergedOuterTracks.size(); ++i){
      //   }
      // }

      auto innerIt = std::find(unmergedInnerTracks.begin(), unmergedInnerTracks.end(), innerTrack);
      auto outerIt = std::find(unmergedOuterTracks.begin(), unmergedOuterTracks.end(), outerTrack);
      if (innerIt != unmergedInnerTracks.end() && outerIt != unmergedOuterTracks.end()){
        std::vector<TS> trackStates = combination.second;
        addCombinedTrackToColl(innerTrack, outerTrack, trackStates, outTracks);
        info() << "track combination added. chi2: " << chi2 << endmsg;
        printTrack(*innerIt, TS::AtLastHit);
        printTrack(*outerIt, TS::AtFirstHit);
        unmergedInnerTracks.erase(innerIt);
        unmergedOuterTracks.erase(outerIt);
      }
    }

    info() << "Unmerged tracks after merging:" << endmsg;
    for (Track& track: unmergedInnerTracks)
      printTrack(track,TS::AtLastHit);
    for (Track& track: unmergedOuterTracks)
      printTrack(track,TS::AtFirstHit);

    if (checkForNoUnmergedTracks(unmergedInnerTracks, unmergedOuterTracks, outTracks))
      return outTracks;

    bool remainingMatches = false;
    std::map<Track, bool> innerTrackHasMatch, outerTrackHasMatch;
    for (const Track& outerTrack : unmergedOuterTracks)
    {outerTrackHasMatch[outerTrack] = false;}
    for (const Track& innerTrack : unmergedInnerTracks){
      innerTrackHasMatch[innerTrack] = false;
      for (const Track& outerTrack : unmergedOuterTracks){
        if (isMatch(innerTrack, TS::AtLastHit, outerTrack, TS::AtFirstHit)){
          remainingMatches = true;
          innerTrackHasMatch[innerTrack] = true;
          outerTrackHasMatch[outerTrack] = true;
        }
      }
    }

    for (Track& trk : unmergedInnerTracks)
      if (!innerTrackHasMatch[trk])
        addTrackToColl(trk, outTracks);
    for (Track& trk : unmergedOuterTracks)
      if (!outerTrackHasMatch[trk])
        addTrackToColl(trk, outTracks);

    if (!remainingMatches)
      return outTracks;

    warning() << "There are " << unmergedInnerTracks.size() << " unmerged inner tracks and " << unmergedOuterTracks.size() << " unmerged outer tracks with matches between them!" << endmsg;
    for (Track& track: unmergedInnerTracks)
      printTrack(track,TS::AtLastHit);
    for (Track& track: unmergedOuterTracks)
      printTrack(track,TS::AtFirstHit);

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
                               static_cast<int>(loc2)) << endmsg;
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
        // Negated on purpose: NaN compares false, so `<= 0.f` would not catch it
        if (tolerance >= 0.f && !(ts1->getCovMatrix(param, param) + ts2->getCovMatrix(param, param) > 0.f)) {
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

    info() << fmt::format("Comparing Loc {} vs {} ({}): d0={:.4f}, z0={:.4f}, phi={:.4f}, omega={:.4f}, "
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
  GaudiDDKalTest m_ddkaltest{this};
  dd4hep::DDSegmentation::BitFieldCoder m_encoder;

  void printTrack(const Track& t1, TS::Location loc1) const {
    auto ts1 = t1.getTrackState(loc1);

    info() << fmt::format("Track 1: d0={:.4f}, z0={:.4f}, phi={:.4f}, omega={:.4f}, tanLambda={:.4f}",
                          ts1->D0, ts1->Z0, ts1->phi, ts1->omega, ts1->tanLambda) << endmsg;
    //info() << fmt::format("Track 2: d0={:.4f}, z0={:.4f}, phi={:.4f}, omega={:.4f}, tanLambda={:.4f}",
    //                      ts2->D0, ts2->Z0, ts2->phi, ts2->omega, ts2->tanLambda) << endmsg;
  }

  // check if either or both of the remaining track vectors are empty
  // if one is empty, the other can be fully added to the output collection
  bool checkForNoUnmergedTracks(std::vector<Track> tv1, std::vector<Track> tv2, TrackColl& outTracks) const {
    if (tv1.size()==0 && tv2.size()==0)
      return 1;
    if (tv1.size()==0){
      for (Track& trk : tv2)
        addTrackToColl(trk, outTracks);
      return 1;
    }
    if (tv2.size()==0){
      for (Track& trk : tv1)
        addTrackToColl(trk, outTracks);
      return 1;
    }
    return 0;
  }

  void addTrackToColl(const Track& aTrack, TrackColl& outTracks) const {
    auto newTrack = outTracks.create();
    for (const TH& hit : aTrack.getTrackerHits())
      newTrack.addToTrackerHits(hit);
    newTrack.addToTracks(aTrack);
    for (TS ts : aTrack.getTrackStates())
      newTrack.addToTrackStates(ts);
  }

  void addCombinedTrackToColl(const Track& track1, const Track& track2, std::vector<TS> trackStates, TrackColl& outTracks) const {
    auto newTrack = outTracks.create();
    for (const TH& hit : track1.getTrackerHits())
      newTrack.addToTrackerHits(hit);
    newTrack.addToTracks(track1);
    for (const TH& hit : track2.getTrackerHits())
      newTrack.addToTrackerHits(hit);
    newTrack.addToTracks(track2);
    for (TS ts : trackStates)
      newTrack.addToTrackStates(ts);
  }

};

DECLARE_COMPONENT(TrackMerger)
