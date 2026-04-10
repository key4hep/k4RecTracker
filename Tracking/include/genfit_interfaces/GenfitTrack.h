/*
 * Copyright (c) 2020-2024 Key4hep-Project.
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

#ifndef GENFIT_TRACK_H
#define GENFIT_TRACK_H

#include <memory>
#include <ranges>
#include <stdexcept>
#include <string>
#include <vector>

#include <TVector3.h>

#include <AbsBField.h>
#include <AbsFitterInfo.h>
#include <AbsTrackRep.h>
#include <DAF.h>
#include <EventDisplay.h>
#include <FieldManager.h>
#include <GenfitTrack.h>
#include <KalmanFitter.h>
#include <KalmanFitterInfo.h>
#include <KalmanFitterRefTrack.h>
#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <Track.h>

#include "DD4hep/Detector.h"
#include "DD4hep/Fields.h"
#include "DDRec/Vector3D.h"
#include "DDSegmentation/BitFieldCoder.h"

#include "edm4hep/MutableTrack.h"
#include "edm4hep/TrackCollection.h"

#include "FastCircleSeed.h"
#include "GenfitField.h"
#include "GenfitPlanarMeasurement.h"
#include "GenfitWireMeasurement.h"
#include "utils.h"

/** @class GenfitTrack
 *
 *  Internal helper class that bridges an EDM4hep track with the GENFIT track representation.
 *  This class is responsible for preparing, initializing, and managing the GENFIT track object,
 *  starting from an EDM4hep `edm4hep::Track`, and performing the track fit using the GENFIT library.
 *
 *  The `GenfitTrack` encapsulates the logic required to:
 *    - extract and convert the track parameters (position, momentum, charge)
 *    - create and configure the appropriate `genfit::AbsTrackRep` based on the particle hypothesis
 *    - construct the `genfit::Track` object including measurement hits
 *    - execute the GENFIT fitting procedure
 *
 *  The class maintains access to both the EDM4hep and GENFIT representations of the track, and stores
 *  relevant geometry information such as the drift chamber (DCH) description and segmentation decoder.
 *
 *
 *  Author: Andrea De Vita
 *  Date  : 2025-11
 *
 */

namespace GenfitInterface {

class GenfitTrack {
public:
  GenfitTrack(const edm4hep::Track& track, const bool skipTrackOrdering = false,
              const dd4hep::rec::DCH_info* dch_info = nullptr,
              const dd4hep::DDSegmentation::BitFieldCoder* decoder = nullptr,
              const GenfitInterface::GenfitField* fieldMap = nullptr);

  ~GenfitTrack();

  void InitializeTrack(double RadiusForDisplacedTracking, bool UseFirstHitAsReference, bool LimitHits,
                       int InitializationType, std::optional<int> TrackStateLocation,
                       std::optional<TVector3> Init_position, std::optional<TVector3> Init_momentum,
                       std::optional<double> Epsilon, std::optional<int> Window);

  void CreateGenFitTrack(int particle_hypotesis, int debug_lvl);
  bool Fit(std::string FitterType, int debug_lvl, std::optional<double> Beta_init, std::optional<double> Beta_final,
           std::optional<int> Beta_steps, std::optional<bool> FilterHits);

  genfit::Track* GetTrack_genfit() { return m_genfitTrack; }
  genfit::AbsTrackRep* GetRep_genfit() { return m_genfitTrackRep; }
  edm4hep::MutableTrack& GetTrack_edm4hep() { return m_edm4hepTrack; }
  edm4hep::MutableTrack& GetTrackWithFit_edm4hep() { return m_trackWithFit; }
  edm4hep::TrackerHitPlaneCollection& GetFittedHits() { return m_fittedHits; }

  int GetCharge() { return m_charge_hypothesis; }

  void PrintTrack_init() {

    std::cout << "GENFIT Initial position: (" << m_posInit.X() << ", " << m_posInit.Y() << ", " << m_posInit.Z() << ")"
              << std::endl;
    std::cout << "GENFIT Initial momentum: (" << m_momInit.X() << ", " << m_momInit.Y() << ", " << m_momInit.Z() << ")"
              << std::endl;
  }

  struct HelperInitialization {
    TVector3 Position;
    TVector3 Momentum;
    TMatrixDSym CovMatrix;
    int Charge;
    int NumHits;
  };

  HelperInitialization GetInitialization() {

    return {m_posInit, m_momInit, m_covInit, m_charge_hypothesis,
            static_cast<int>(m_edm4hepTrack.getTrackerHits().size())};
  }

private:
  edm4hep::Track m_originalTrack;

  struct PCAInfoHelper {
    TVector3 PCA;
    int Phi0;
  };

  void CheckInitialization();
  void OrderHits(const edm4hep::Track& track, bool skipTrackOrdering);
  void LimitNumberHits(double epsilon, int smoothWindow);

  TMatrixDSym CovarianceMatrixHelixToCartesian(const TMatrixDSym& C_helix, TVector3 Position_cm, TVector3 Momentum_gev,
                                               TVector3 RefPoint_cm, double Bz);

  TMatrixDSym ComputeInitialCovarianceMatrix(double Bz);

  HelperInitialization ComputeInitialParameters(double Bz);
  PCAInfoHelper PCAInfo(TVector3 position, TVector3 momentum, int charge, TVector3 refPoint, double Bz);

  int m_signed_particle_hypothesis = 211;
  int m_charge_hypothesis = 1;

  TVector3 m_posInit = TVector3(0., 0., 0.);
  TVector3 m_momInit = TVector3(0., 0., 0.);
  TMatrixDSym m_covInit;

  genfit::AbsTrackRep* m_genfitTrackRep;
  genfit::Track* m_genfitTrack;
  edm4hep::MutableTrack m_edm4hepTrack;
  edm4hep::MutableTrack m_trackWithFit;
  edm4hep::TrackerHitPlaneCollection m_fittedHits;

  TVector3 m_VP_referencePoint{0., 0., 0.};

  const dd4hep::rec::DCH_info* m_dch_info;
  const dd4hep::DDSegmentation::BitFieldCoder* m_dc_decoder;
  const GenfitInterface::GenfitField* m_fieldMap;
};

} // namespace GenfitInterface

#endif // GenfitTrack