/*
 * Copyright (c) 2014-2024 Key4hep-Project.
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

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#include <TGeoManager.h>
#include <TVector3.h>

#include <Exception.h>
#include <FieldManager.h>
#include <MaterialEffects.h>
#include <TGeoMaterialInterface.h>

#include "Gaudi/Property.h"
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/Fields.h"
#include "DD4hep/Readout.h"
#include <DDRec/DetectorData.h>
#include <DDRec/MaterialManager.h>
#include <DDRec/SurfaceManager.h>
#include <DDRec/Vector3D.h>
#include <DDSegmentation/BitFieldCoder.h>

#include "edm4hep/SenseWireHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "GenfitField.h"
#include "GenfitMaterialInterface.h"
#include "GenfitTrack.h"

#include "utils.h"

/** @struct GenfitTrackFitter
 *
 *  Gaudi MultiTransformer that refines the parameters of the reconstructed tracks using the GENFIT library.
 *  The fitting process minimizes the track χ² while accounting for material effects, magnetic field, and detector
 * geometry.
 *
 *  If the fit fails for a given hypothesis, a fallback track is still generated with chi2 = -1 and ndf = -1 to preserve
 * collection integrity. Each fitted track also includes the extrapolated state at the calorimeter (barrel or endcap).
 *
 *  input:
 *    - initial track collection : edm4hep::TrackCollection
 *
 *  output:
 *    - fitted track collection : edm4hep::TrackCollection
 *    - fitted track collection with filtered hits (without L/R ambiguity if drift chamber is present):
 * edm4hep::TrackCollection
 *    - filtered hit collection : edm4hep::TrackerHitPlaneCollection
 *
 *  @author Andrea De Vita
 *  @date   2026-02
 *
 */

struct GenfitTrackFitter final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackCollection, edm4hep::TrackCollection,
                                            edm4hep::TrackerHitPlaneCollection>(const edm4hep::TrackCollection&)> {
  GenfitTrackFitter(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc,
                         {

                             KeyValues("InputTracks", {"InputTracks"})

                         },
                         {

                             KeyValues("OutputFittedTracks", {"Fitted_tracks"}),
                             KeyValues("OutputFittedTracksWithFilteredHits", {"Fitted_tracks_with_filtered_hits"}),
                             KeyValues("OutputFittedHits", {"Fitted_hits"})

                         }) {}

  StatusCode initialize() {

    // Initialize printout level for debug information
    m_printoutLevel = msgLevel();

    if (!gGeoManager) {
      error() << "Error: TGeoManager is not initialized!" << endmsg;
      return StatusCode::FAILURE;
    }

    // Initialize the Genfit: FieldManager and MaterialEffects
    m_detector = m_geoSvc->getDetector();
    auto fieldMap = m_detector->field();
    m_genfitField = std::make_unique<GenfitInterface::GenfitField>(fieldMap);

    m_fieldManager = genfit::FieldManager::getInstance();
    m_fieldManager->init(m_genfitField.get());

    m_geoMaterial = GenfitInterface::GenfitMaterialInterface::getInstance(m_detector);

    if (m_useBrems) {
      genfit::MaterialEffects::getInstance()->setEnergyLossBrems(true);
      genfit::MaterialEffects::getInstance()->setNoiseBrems(true);
      genfit::MaterialEffects::getInstance()->setMscModel("Highland");
    } else {
      genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
      genfit::MaterialEffects::getInstance()->setNoiseBrems(false);
    }

    int debug_lvl_material = 0;
    if (m_printoutLevel == uint(MSG::VERBOSE)) {
      debug_lvl_material = 1;
    }
    genfit::MaterialEffects::getInstance()->setDebugLvl(debug_lvl_material);

    // Retrieve the SurfaceManager, ddd4hep::rec::DCH_info and dd4hep::DDSegmentation::BitFieldCoder
    // These object are necessary to extract the drift chamber hits information, such as positions of the wire
    // extremities
    m_surfMan = m_geoSvc->getDetector()->extension<dd4hep::rec::SurfaceManager>();
    // If the detector doesn't have a drift chamber, this part will be skipped
    try {

      std::string DCH_name(m_DCH_name.value());
      dd4hep::DetElement DCH_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);
      m_dch_info = DCH_DE.extension<dd4hep::rec::DCH_info>();
      dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector(DCH_name);
      dd4hep::Readout dch_readout = dch_sd.readout();
      m_dc_decoder = dch_readout.idSpec().decoder();

    } catch (const std::out_of_range& e) {

      warning() << "No DCH found in detector geometry (name = " << m_DCH_name.value()
                << "). Skipping DCH-related setup." << endmsg;
    }

    // Retrive calorimeter information
    // These parameters are necessary to propagate the track to the calorimeter surface
    dd4hep::Detector& detector = dd4hep::Detector::getInstance();
    const std::vector<dd4hep::DetElement>& isDuaReadout = dd4hep::DetectorSelector(detector).detectors(
        dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL | dd4hep::DetType::ENDCAP, dd4hep::DetType::AUXILIARY);
    if (isDuaReadout.size() > 0) {

      const dd4hep::rec::LayeredCalorimeterData* DualReadoutExtension =
          getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
                       (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));

      // mm (dd4hep::mm = 0.1)
      m_eCalBarrelInnerR = DualReadoutExtension->extent[0] / dd4hep::mm; // barrel rmin
      m_eCalBarrelMaxZ = DualReadoutExtension->extent[2] / dd4hep::mm;   // barrel zmax == endcap zmin

      m_eCalEndCapInnerR = DualReadoutExtension->extent[4] / dd4hep::mm; // endcap rmin
      m_eCalEndCapOuterR = DualReadoutExtension->extent[5] / dd4hep::mm; // endcap rmax
      m_eCalEndCapInnerZ = DualReadoutExtension->extent[2] / dd4hep::mm; // endcap zmin
      m_eCalEndCapOuterZ = DualReadoutExtension->extent[3] / dd4hep::mm; // endcap zmax
    } else {
      const dd4hep::rec::LayeredCalorimeterData* eCalBarrelExtension =
          getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
                       (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));

      // mm (dd4hep::mm = 0.1)
      m_eCalBarrelInnerR = eCalBarrelExtension->extent[0] / dd4hep::mm;
      m_eCalBarrelMaxZ = eCalBarrelExtension->extent[3] / dd4hep::mm;

      const dd4hep::rec::LayeredCalorimeterData* eCalEndCapExtension =
          getExtension((dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP),
                       (dd4hep::DetType::AUXILIARY | dd4hep::DetType::FORWARD));

      // mm (dd4hep::mm = 0.1)
      m_eCalEndCapInnerR = eCalEndCapExtension->extent[0] / dd4hep::mm;
      m_eCalEndCapOuterR = eCalEndCapExtension->extent[1] / dd4hep::mm;
      m_eCalEndCapInnerZ = eCalEndCapExtension->extent[2] / dd4hep::mm;
      m_eCalEndCapOuterZ = eCalEndCapExtension->extent[3] / dd4hep::mm;
    }

    return StatusCode::SUCCESS;
  }

  std::tuple<edm4hep::TrackCollection, edm4hep::TrackCollection, edm4hep::TrackerHitPlaneCollection>
  operator()(const edm4hep::TrackCollection& tracks_input) const override {

    debug() << "Event number: " << event_counter++ << endmsg;

    // This collection stores the output of the fit
    edm4hep::TrackCollection FittedTracks;
    edm4hep::TrackCollection FittedTracksWithFilteredHits;
    edm4hep::TrackerHitPlaneCollection FittedHits;

    // Loop over the tracks created by the pattern recognition step
    for (const auto& track : tracks_input) {

      num_tracks += 1;

      // Skip unmatched tracks if the option is enabled
      // Consider unmatched tracks those with type = 0
      if (m_skip_unmatchedTracks && track.getType() == 0) {
        num_skip += 1;
        warning() << "Track " << num_tracks - 1 << ": unmatched track (type = 0), skipping fit.\n" << endmsg;
        continue; // skip unmatched tracks
      }

      if (track.getTrackerHits().size() < 3) {
        num_skip += 1;
        warning() << "Track " << num_tracks - 1 << ": less than 3 hits, skipping fit.\n" << endmsg;
        continue; // skip tracks with less then 3 hits (seed initialization needs 3 hits)
      }

      num_processed_tracks += 1;

      int isSuccess = 0;
      if (m_singleEvaluation) {

        isSuccess =
            ProcessTrack(track, false, m_particleHypothesis[0], FittedTracks, FittedTracksWithFilteredHits, FittedHits);

      } else {

        int winning_hypothesis = FindBestHypothesis(track, false);

        if (winning_hypothesis == -1) {
          debug() << "Track " << num_tracks - 1 << ": fit failed for all hypotheses, trying with less hits." << endmsg;
        } else {

          isSuccess = 1;
          int pdgCode = winning_hypothesis;
          ProcessTrack(track, false, pdgCode, FittedTracks, FittedTracksWithFilteredHits, FittedHits);
        }
      }

      if (!isSuccess) {

        if (m_singleEvaluation) {

          isSuccess = ProcessTrack(track, true, m_particleHypothesis[0], FittedTracks, FittedTracksWithFilteredHits,
                                   FittedHits);

          if (!isSuccess) {

            number_failures += 1;
            debug() << "Track " << num_tracks - 1 << ": fit failed for single evaluation hypothesis, skipping track."
                    << endmsg;
            auto failedTrack = FittedTracks.create();
            failedTrack.setChi2(-1);
            failedTrack.setNdf(-1);
            continue;
          }

        } else {

          int winning_hypothesis = FindBestHypothesis(track, true);

          if (winning_hypothesis == -1) {

            debug() << "Track " << num_tracks - 1 << ": fit failed for all hypotheses." << endmsg;
            number_failures += 1;
            auto failedTrack = FittedTracks.create();
            failedTrack.setChi2(-1);
            failedTrack.setNdf(-1);
            continue;

          } else {

            int pdgCode = winning_hypothesis;
            ProcessTrack(track, true, pdgCode, FittedTracks, FittedTracksWithFilteredHits, FittedHits);
          }
        }
      }
    }

    return std::make_tuple(std::move(FittedTracks), std::move(FittedTracksWithFilteredHits), std::move(FittedHits));
  }

  StatusCode finalize() {

    info() << "Run report:" << endmsg;
    info() << "Number of tracks: " << num_tracks << endmsg;
    info() << "Number of successes: " << (num_processed_tracks - number_failures) << "/" << num_processed_tracks
           << endmsg;
    info() << "Number of skipped tracks: " << num_skip << "/" << num_tracks << endmsg;
    info() << "----------------\n" << endmsg;

    return StatusCode::SUCCESS;
  }

public:
  // Debug / monitoring counters (not part of core algorithm logic).
  // These are mainly used to keep track of what happens during processing,
  // e.g. for logging, sanity checks, and performance evaluation.

  mutable int event_counter = 0; // Counts how many events have been processed

  // Num_tracks = num_processed_tracks + num_skip
  mutable int num_tracks = 0;           // Total number of tracks seen (including skipped ones)
  mutable int num_skip = 0;             // Number of tracks skipped (e.g. failing pre-selection)
  mutable int num_processed_tracks = 0; // Number of tracks actually processed (i.e. not skipped)
  mutable int number_failures = 0;      // Number of track fits that failed

private:
  // ====================== Debug & Printout ======================
  // Debug level for track fitting and initialization printouts
  // 0       : no printouts
  // INFO    : prints fit results
  // DEBUG   : prints fit results + initial track parameters + track states
  // VERBOSE : prints detailed internal fit information
  uint m_printoutLevel;

  // ====================== Geometry & Detector ======================
  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Detector geometry service"};
  dd4hep::Detector* m_detector{nullptr}; // Detector instance

  std::unique_ptr<GenfitInterface::GenfitField> m_genfitField;
  genfit::FieldManager* m_fieldManager;

  GenfitInterface::GenfitMaterialInterface* m_geoMaterial;

  dd4hep::rec::SurfaceManager* m_surfMan;
  dd4hep::rec::DCH_info* m_dch_info;
  dd4hep::DDSegmentation::BitFieldCoder* m_dc_decoder;

  Gaudi::Property<std::string> m_DCH_name{
      this, "DCHName", "DCH_v2",
      "DCHName in the detector description (used to retrieve DCH geometry and material information)"};

  // ====================== ECAL Geometry Parameters ======================
  double m_eCalBarrelInnerR;
  double m_eCalBarrelMaxZ;
  double m_eCalEndCapInnerR;
  double m_eCalEndCapOuterR;
  double m_eCalEndCapInnerZ;
  double m_eCalEndCapOuterZ;

  // ====================== Fitter Settings ======================

  Gaudi::Property<bool> m_useBrems{this, "UseBrems", false, "Include Bremsstrahlung energy loss and noise in the fit"};

  Gaudi::Property<std::string> m_Fitter_type{
      this, "FitterType", "DAF",
      "Fitter type to use [https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf]: 'DAF', 'KALMAN', "
      "or 'KALMAN_REF'"};
  Gaudi::Property<double> m_Beta_init{this, "BetaInit", 100.,
                                      "Initial annealing parameter (if m_Fitter_type == 'DAF') "
                                      "[https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf]"};
  Gaudi::Property<double> m_Beta_final{this, "BetaFinal", 0.05,
                                       "Final annealing parameter (if m_Fitter_type == 'DAF') "
                                       "[https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf]"};
  Gaudi::Property<int> m_Beta_steps{this, "BetaSteps", 15,
                                    "Number of steps in the annealing schedule (if m_Fitter_type == 'DAF') "
                                    "[https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf]"};

  // ====================== Track Initialization ======================
  Gaudi::Property<std::vector<int>> m_particleHypothesis{
      this,
      "ParticleHypothesisList",
      {211},
      "List of particle hypotheses to consider: {11,13,211,321,2212} -> e, mu, pi, K, p"};

  Gaudi::Property<bool> m_useFirstHitAsReference{
      this, "UseFirstHitAsReference", false,
      "If true, the first hit is used as the reference point for track initialization"};

  Gaudi::Property<int> m_initializationType{
      this, "InitializationType", 1,
      "Method for initializing track parameters before fitting: "
      "0: first two hits (position from first hit, momentum from direction between first and second hit); "
      "1: refined initialization using ComputeInitialParameters(Bz); "
      "2: initialize from a track state (position from reference point, momentum from helix parameters); "
      "3: use user-provided InitPosition and InitMomentum"};

  Gaudi::Property<int> m_trackStateLocation{this, "TrackStateLocation", edm4hep::TrackState::AtFirstHit,
                                            "TrackState location used for initialization when InitializationType = 2. "
                                            "Defines where the reference point and helix parameters are taken from"};

  Gaudi::Property<std::vector<double>> m_init_position{
      this, "InitPosition", {0., 0., 0.}, "User-defined initial position for InitializationType = 3"};

  Gaudi::Property<std::vector<double>> m_init_momentum{
      this, "InitMomentum", {0., 0., 0.}, "User-defined initial momentum for InitializationType = 3"};

  Gaudi::Property<double> m_epsilon{
      this, "Epsilon", 1e-4,
      "Threshold to detect sign changes in the first derivative during track initialization. "
      "Used to identify hits from the first round of the looper"};

  Gaudi::Property<int> m_smoothWindow{
      this, "SmoothWindow", 5,
      "Number of hits used for smoothing before computing the first derivative in track initialization"};

  // ====================== Track Filtering & Evaluation ======================
  Gaudi::Property<bool> m_skipTrackOrdering{this, "SkipTrackOrdering", false, "Skip hit ordering before fitting"};

  Gaudi::Property<bool> m_skip_unmatchedTracks{this, "SkipUnmatchedTracks", true,
                                               "Skip tracks with type = 0 (unmatched) from fitting"};

  Gaudi::Property<bool> m_singleEvaluation{
      this, "RunSingleEvaluation", false,
      "If true, only the first particle hypothesis is evaluated. "
      "If false, all hypotheses are scanned and the best fit is chosen based on chi2/ndf"};

  Gaudi::Property<bool> m_filterTrackHits{
      this, "FilterTrackHits", false,
      "Filter track hits after fitting based on their quality (weights assigned by the fitter). "
      "This also resolves the left/right ambiguity in drift chamber hits: only the hit with the higher weight is "
      "considered."};

  // ====================== Track Classification ======================
  Gaudi::Property<double> m_identifyDisplaced{
      this, "IdentifyDisplacedRadius", 100.,
      "Radius [mm] defining a spherical region around {0,0,0}: "
      "tracks produced within this radius are considered prompt; tracks outside are considered displaced"};

  // Fit implementation
  bool ProcessTrack(const edm4hep::Track& track, bool LimitHits, int particleHypothesis,
                    edm4hep::TrackCollection& FittedTracks, edm4hep::TrackCollection& FittedTracksWithFilteredHits,
                    edm4hep::TrackerHitPlaneCollection& FittedHits) const {

    GenfitInterface::GenfitTrack track_interface(track, m_skipTrackOrdering, m_dch_info, m_dc_decoder,
                                                 m_genfitField.get());

    TVector3 Init_position(m_init_position.value()[0], m_init_position.value()[1], m_init_position.value()[2]);

    TVector3 Init_momentum(m_init_momentum.value()[0], m_init_momentum.value()[1], m_init_momentum.value()[2]);

    track_interface.InitializeTrack(m_identifyDisplaced, m_useFirstHitAsReference, LimitHits, m_initializationType,
                                    m_trackStateLocation.value(), Init_position, Init_momentum, m_epsilon.value(),
                                    m_smoothWindow.value());

    auto track_init = track_interface.GetInitialization();

    debug() << "Track " << num_tracks - 1 << " with " << track.getTrackerHits().size()
            << " hits: initial seed for track fit:" << endmsg;

    debug() << "  Initial position [mm]: (" << track_init.Position.X() / dd4hep::mm << ", "
            << track_init.Position.Y() / dd4hep::mm << ", " << track_init.Position.Z() / dd4hep::mm << ")" << endmsg;

    debug() << "  Initial momentum [GeV/c]: (" << track_init.Momentum.X() << ", " << track_init.Momentum.Y() << ", "
            << track_init.Momentum.Z() << ")" << endmsg;

    debug() << "  Charge hypothesis: " << track_init.Charge << endmsg;
    debug() << "  Max hit for loopers: " << track_init.NumHits << endmsg;

    if (m_printoutLevel == uint(MSG::DEBUG)) {
      debug() << "  Initial covariance matrix:" << endmsg;
      track_init.CovMatrix.Print();
    }

    debug() << endmsg;

    int debug_track = (m_printoutLevel == uint(MSG::DEBUG)) ? 1 : 0;

    track_interface.CreateGenFitTrack(particleHypothesis, debug_track);

    bool isFit = track_interface.Fit(m_Fitter_type.value(), m_printoutLevel, m_Beta_init, m_Beta_final, m_Beta_steps,
                                     m_filterTrackHits);

    if (!isFit) {
      debug() << "Track fit FAILED for track " << num_tracks - 1 << endmsg;
      return false;
    }

    auto edm4hep_track = track_interface.GetTrack_edm4hep();

    double Bz = 0.;
    for (auto ts : edm4hep_track.getTrackStates()) {
      if (ts.location == edm4hep::TrackState::AtLastHit) {

        auto ref = ts.referencePoint;
        Bz = m_genfitField->getBz(TVector3(ref.x * dd4hep::mm, ref.y * dd4hep::mm, ref.z * dd4hep::mm)) / 10.;
      }
    }

    if (Bz > 0) {

      FillTrackWithCalorimeterExtrapolation(edm4hep_track, Bz, track_interface.GetCharge(), m_eCalBarrelInnerR,
                                            m_eCalBarrelMaxZ, m_eCalEndCapInnerR, m_eCalEndCapOuterR,
                                            m_eCalEndCapInnerZ);
    }

    info() << "Fit chi2/ndf = " << edm4hep_track.getChi2() << " / " << edm4hep_track.getNdf() << " = "
           << edm4hep_track.getChi2() / (edm4hep_track.getNdf() * 1.) << endmsg;

    FittedTracks.push_back(edm4hep_track);

    if (m_filterTrackHits) {

      auto edm4hep_track_with_fit = track_interface.GetTrackWithFit_edm4hep();

      if (Bz > 0) {

        FillTrackWithCalorimeterExtrapolation(edm4hep_track_with_fit, Bz, track_interface.GetCharge(),
                                              m_eCalBarrelInnerR, m_eCalBarrelMaxZ, m_eCalEndCapInnerR,
                                              m_eCalEndCapOuterR, m_eCalEndCapInnerZ);
      }

      FittedTracksWithFilteredHits.push_back(edm4hep_track_with_fit);
      FittedHits = std::move(track_interface.GetFittedHits());
    }

    return true;
  }

  int FindBestHypothesis(const edm4hep::Track& track, bool LimitHits) const {

    TVector3 Init_position(m_init_position.value()[0], m_init_position.value()[1], m_init_position.value()[2]);

    TVector3 Init_momentum(m_init_momentum.value()[0], m_init_momentum.value()[1], m_init_momentum.value()[2]);

    int winning_hypothesis = -1;
    double winning_chi2_ndf = std::numeric_limits<double>::max();

    for (int pdgCode : m_particleHypothesis) {

      GenfitInterface::GenfitTrack track_interface(track, m_skipTrackOrdering, m_dch_info, m_dc_decoder,
                                                   m_genfitField.get());

      track_interface.InitializeTrack(m_identifyDisplaced, m_useFirstHitAsReference, LimitHits, m_initializationType,
                                      m_trackStateLocation.value(), Init_position, Init_momentum, m_epsilon.value(),
                                      m_smoothWindow.value());

      track_interface.CreateGenFitTrack(pdgCode, 0);

      bool isFit = track_interface.Fit(m_Fitter_type.value(), 0, m_Beta_init, m_Beta_final, m_Beta_steps, false);

      if (!isFit)
        continue;

      auto edm4hep_track = track_interface.GetTrack_edm4hep();

      float chi2 = edm4hep_track.getChi2();
      int ndf = edm4hep_track.getNdf();

      if (chi2 <= 0 || ndf <= 0)
        continue;

      double chi2_ndf = chi2 / ndf;

      if (chi2_ndf < winning_chi2_ndf) {

        winning_chi2_ndf = chi2_ndf;
        winning_hypothesis = pdgCode;
      }
    }

    return winning_hypothesis;
  }
};

DECLARE_COMPONENT(GenfitTrackFitter)