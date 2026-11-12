#ifndef COMBINATORIAL_TRACK_FINDER_H
#define COMBINATORIAL_TRACK_FINDER_H

#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Gaudi/Accumulators.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/ServiceHandle.h"

#include "DD4hep/Detector.h"
#include "DD4hep/Fields.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "DDSegmentation/BitFieldCoder.h"

#include "LayerGeometry.h"
#include "SeedFinder.h"
#include "TrackingMath.h"
#include <Eigen/Dense>

#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"

#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

/** @class CombinatorialTrackFinder
 *
 *  Track finding on the planar hits of one tracking system (inner tracker, silicon wrapper,
 *  muon system, ...), with no input from other detectors. Per event:
 *   1. hits are grouped by layer (barrel / endcap layers from the cellID);
 *   2. candidates are found by
 *        - Exhaustive: all k-hit combinations over the layers, each with a circle fit in xy and a
 *          line fit of z versus arc length, or
 *        - SeedAndExtend: 4-hit then 3-hit seeds, best first, extended layer by layer along a helix
 *          (kink chi2 with multiple scattering; hits claimed as tracks are accepted);
 *      road, outlier and quality cuts reject wrong combinations;
 *   3. the best candidate becomes a track, its hits are removed, and the search repeats;
 *   4. each track stores its helix as AtFirstHit and AtLastHit, optionally AtIP (IPStateMode).
 *  For a Kalman fit, run k4RecTracker's GenfitTrackFitter on the output.
 *
 *  @author Mahmoud Althakeel
 *  @email  mahmoud.althakeel@cern.ch
 */
class CombinatorialTrackFinder final
    : public k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackCollection, edm4hep::TrackerHitPlaneCollection,
                                                   edm4hep::TrackerHitSimTrackerHitLinkCollection>(
          const edm4hep::TrackerHitPlaneCollection&, const edm4hep::EventHeaderCollection&,
          const edm4hep::TrackerHitSimTrackerHitLinkCollection&)> {
public:
  CombinatorialTrackFinder(const std::string& name, ISvcLocator* pSvcLocator);

  StatusCode initialize() override;

  std::tuple<edm4hep::TrackCollection, edm4hep::TrackerHitPlaneCollection,
             edm4hep::TrackerHitSimTrackerHitLinkCollection>
  operator()(const edm4hep::TrackerHitPlaneCollection& hits, const edm4hep::EventHeaderCollection& headers,
             const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks) const override;

private:
  class EventSearch; ///< the track search in one event (CombinatorialTrackFinder.cpp)

  enum class Strategy { Exhaustive, SeedAndExtend };    ///< option Strategy
  enum class IPState { None, Helix, SolenoidCrossing }; ///< option IPStateMode

  /// The outer helix of a track measured outside the solenoid (SolenoidCrossing AtIP state), cm.
  struct OuterHelix {
    double centreXCm = 0.0, centreYCm = 0.0, radiusCm = 0.0;  ///< fitted circle
    bool clockwise = false;                                   ///< seen from +z
    double ptGeV = 0.0;                                       ///< radius and field at the hits' centroid
    double tanLambda = 0.0;                                   ///< kept across the solenoid boundary
    Eigen::Vector3d innermostHitCm = Eigen::Vector3d::Zero(); ///< selects the crossing, anchors z
    bool innermostHitInEndcap = false;                        ///< then: the solenoid's end face
  };

  // ---- the steps of initialize() ----
  StatusCode setupCellIDDecoding();
  StatusCode setupSurfaces();
  StatusCode readSearchOptions();
  StatusCode setupSeedFinder();
  void findSolenoid();

  // ---- detector surfaces and cellID decoding ----
  const dd4hep::rec::Surface* surfaceOf(uint64_t cellID) const; ///< nullptr if none
  int regionOf(uint64_t cellID) const;                          ///< 0 barrel, +1 / -1 endcaps (RegionFieldName)
  int layerNumberOf(uint64_t cellID) const;                     ///< LayerFieldName
  int layerKeyOf(uint64_t cellID) const;                        ///< layout in LayerGeometry.h
  bool readoutHasField(const std::string& name) const;
  std::map<int, std::vector<const dd4hep::rec::Surface*>> groupSurfacesByLayer() const;

  // ---- the solenoid ----
  /// where the central field ends (mm), followed outward from near the IP: radius (alongZ = false)
  /// or half-length (alongZ = true)
  double findSolenoidBoundary(bool alongZ) const;
  double fieldAlongScan(double distanceMm, bool alongZ) const; ///< Bz (T) on the scan line
  /// the field has moved more than half-way from bzCentral towards zero (or beyond)
  bool outsideCentralField(double distanceMm, bool alongZ, double bzCentral) const;
  bool isInsideSolenoid(const Eigen::Vector3d& positionCm) const;
  bool allHitsOutsideSolenoid(const std::vector<edm4hep::TrackerHitPlane>& trackHits) const;
  /// SolenoidCrossing: the outer helix followed back from its innermost hit to the solenoid
  /// boundary, matched there to an inner helix (inner field, same pT and tanLambda) whose perigee is
  /// the AtIP state; no value without a crossing
  std::optional<edm4hep::TrackState> solenoidCrossingState(const OuterHelix& outer) const;
  /// crossing (cm) of the outer circle with the solenoid barrel, first reached going back from the
  /// hit and on its IP side in z (else the other one, with tryOtherCrossing); false if none
  bool solenoidBarrelCrossing(const OuterHelix& outer, bool tryOtherCrossing, Eigen::Vector3d& crossingCm) const;

  void printTrackSummary(const edm4hep::TrackCollection& trackCollection) const; ///< DEBUG, one line per track

  // ---- settings filled in initialize() ----
  k4rt::CommonSearchConfig m_commonConfig; ///< shared by both strategies
  Strategy m_strategyChoice = Strategy::Exhaustive;
  IPState m_ipState = IPState::None;
  std::vector<std::pair<int, int>> m_doubletLayerPairs; ///< from DoubletLayers
  std::map<int, int> m_partnerLayer;                    ///< layer number -> other layer of its doublet
  double m_solenoidRadiusMm = 0.0;
  double m_solenoidHalfLengthMm = 0.0;

  // ---- options used by both strategies ----
  Gaudi::Property<std::vector<std::string>> m_detectorNames{
      this,
      "DetectorNames",
      {},
      "DetElements whose hits are tracked: one (e.g. [\"Muon-System\"]) or several sharing one readout "
      "(e.g. [\"Vertex\", \"SiWrB\", \"SiWrD\"]); their surface maps are merged"};
  Gaudi::Property<std::string> m_encodingStringParameter{this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
                                                         "DD4hep parameter with the cellID encoding string"};
  Gaudi::Property<std::string> m_layerFieldName{this, "LayerFieldName", "layer", "CellID field of the layer index"};
  Gaudi::Property<std::string> m_regionFieldName{
      this, "RegionFieldName", "type",
      "CellID field of barrel (0) / +endcap (1) / -endcap (-1): \"type\" (MuonSystemReadoutID), "
      "\"side\" (GlobalTrackerReadoutID)"};
  Gaudi::Property<std::string> m_strategy{
      this, "Strategy", "Exhaustive",
      "Exhaustive (all k-hit combinations; sparse events, e.g. a muon system) | SeedAndExtend (seeds "
      "extended along their helix; dense events)"};
  Gaudi::Property<std::string> m_ipStateMode{
      this, "IPStateMode", "None",
      "AtIP state: None, Helix (the helix at its closest approach to the origin; uniform field, e.g. "
      "an inner tracker) or SolenoidCrossing (matched at the solenoid boundary; tracks outside it)"};
  Gaudi::Property<int> m_minTrackHits{this, "MinTrackHits", 3, "Fewest hits of a track"};
  Gaudi::Property<double> m_minPt{
      this, "MinComboPT", -1.0,
      "Minimum pT (GeV); Exhaustive: all combinations, SeedAndExtend: tracks of <= ShortTrackMaxHits "
      "hits (random low-pT triplets pass any kink cut); <= 0 = off"};
  Gaudi::Property<double> m_maxPairHitDistance{
      this, "MaxPairHitDistance", 300.0,
      "Largest 3D distance (cm) between two hits of a combination (300 cm allows barrel -> endcap); -1 = off"};
  Gaudi::Property<double> m_maxHitEdepKeV{
      this, "MaxHitEdepKeV", 10.0,
      "Hits above this energy deposit (keV) are not used (delta rays, secondaries); -1 = off"};
  Gaudi::Property<double> m_sigmaHitDefault{this, "SigmaHitDefault", 0.04,
                                            "Hit resolution (cm) when a hit's du/dv is not set"};

  // ---- options of the Exhaustive strategy ----
  Gaudi::Property<int> m_maxCombinatorialHits{
      this, "MaxCombinatorialHits", 7,
      "Exhaustive: most hits per combination (C(10,7) = 120 is fine, C(20,7) = 77520 may be slow)"};
  Gaudi::Property<double> m_maxCircleChi2NDF{this, "NHitMaxChi2NDF", 10.0,
                                             "Exhaustive: largest chi2/ndf of the circle fit of a combination"};
  Gaudi::Property<double> m_maxPt{this, "MaxComboPT", 250.0,
                                  "Exhaustive: largest pT (GeV) of a combination, from its circle radius"};
  Gaudi::Property<double> m_maxConsecutiveDeltaPhi{
      this, "MaxConsecDeltaPhi", 0.5, "Exhaustive: largest |dphi| (rad) between hits of consecutive layers; -1 = off"};
  Gaudi::Property<double> m_maxCombinationPhiSpread{
      this, "MaxComboPhiSpread", 1.0,
      "Exhaustive: largest phi spread (rad) of a combination, relative to its innermost hit; -1 = off"};
  Gaudi::Property<double> m_maxCombinationTanLambdaDeviation{
      this, "MaxComboTanLambdaDev", 0.5,
      "Exhaustive: largest deviation of a consecutive pair's dip (dz / transverse distance) from the "
      "median (0.5 = about 27 degrees); pure barrel or pure endcap combinations; -1 = off"};
  Gaudi::Property<double> m_maxConsecutiveHitDistance{
      this, "MaxConsecutiveHitDistance", 100.0,
      "Exhaustive: largest 3D distance (cm) between consecutive hits (inner -> outer); -1 = off"};
  Gaudi::Property<double> m_hitIsolationCut{
      this, "HitIsolationCut", 2.0,
      "Exhaustive: smallest distance (cm) of a hit from any other unused hit (other layer) not in the "
      "combination (rejects shower hits); 0 = off"};
  Gaudi::Property<double> m_outlierSigma{
      this, "OutlierSigma", 5.0, "Exhaustive: hits beyond this many SigmaHitDefault from the circle are removed"};
  Gaudi::Property<int> m_maxOutlierIterations{this, "MaxOutlierIterations", 10,
                                              "Exhaustive: rounds of outlier removal per track"};
  Gaudi::Property<int> m_maxHitsPerLayer{
      this, "MaxHitsPerCompositeID", 0,
      "Exhaustive: if > 0, only this many hits per layer, best energy-deposit score first (crowded "
      "layers); 0 = all"};
  Gaudi::Property<double> m_edepScoreWeight{
      this, "EdepScoreWeight", 1.0,
      "Exhaustive: weight of the energy-deposit score in the ranking of combinations with the same "
      "number of hits (k = 3: with the radius; k >= 4: chi2/ndf + weight * score / EdepNormKeV); 0 = off"};
  Gaudi::Property<double> m_edepNormKeV{this, "EdepNormKeV", 2.0,
                                        "Exhaustive: energy-deposit deviation (keV) worth one chi2/ndf unit"};
  Gaudi::Property<bool> m_doHitClustering{
      this, "DoHitClustering", false,
      "Exhaustive: hits joined into road clusters (union-find over neighbouring layers with "
      "MaxConsecDeltaPhi and MaxConsecutiveHitDistance), each searched on its own (high occupancy)"};
  Gaudi::Property<bool> m_doCrowdedLayerHitSelection{this, "DoCrowdedLayerHitSelection", false,
                                                     "Exhaustive: guided hit choice in crowded layers"};
  Gaudi::Property<bool> m_requirePairCoincidence{
      this, "RequirePairCoincidence", false,
      "Exhaustive: use the layer doublets (DoubletLayers): pass 1 only hits with a coincident hit in "
      "the paired layer, pass 2 all hits if pass 1 leaves a doublet incomplete; >= 2 complete "
      "doublets required"};
  Gaudi::Property<std::vector<int>> m_doubletLayers{
      this,
      "DoubletLayers",
      {0, 1, 2, 3, 4, 5},
      "Exhaustive (RequirePairCoincidence): layer numbers of the doublets, two by two"};
  Gaudi::Property<double> m_pairCoincidenceMaxDistance{
      this, "PairCoincidenceMaxDist", 25.0,
      "Exhaustive: largest 3D distance (cm) of the two hits of a doublet (25 cm suits the muon system)"};
  Gaudi::Property<double> m_neighbourTrackMaxDist{
      this, "NeighbourTrackMaxDist", 5.0,
      "Exhaustive: a combination with a hit this close (cm) to an accepted track needs "
      "NeighbourTrackMinHits and NeighbourTrackMaxChi2NDF (delta-ray shadows); 0 = off"};
  Gaudi::Property<int> m_neighbourTrackMinHits{this, "NeighbourTrackMinHits", 4,
                                               "Exhaustive: fewest hits of a combination next to an accepted track"};
  Gaudi::Property<double> m_neighbourTrackMaxChi2NDF{
      this, "NeighbourTrackMaxChi2NDF", 3.0, "Exhaustive: largest chi2/ndf of a combination next to an accepted track"};

  // ---- options of the SeedAndExtend strategy ----
  Gaudi::Property<std::string> m_seedLayerPolicy{
      this, "SeedLayerPolicy", "innermost",
      "SeedAndExtend: innermost (first layers from the IP, extended outward), outermost (extended "
      "inward) or consecutive (any run of layers, displaced tracks; extended both ways)"};
  // default tuned on the IDEA and CLD inner trackers (most 3-4 hit tracks in b-bbar are fakes)
  Gaudi::Property<int> m_shortTrackMaxHits{
      this, "ShortTrackMaxHits", 4, "SeedAndExtend: tracks up to this many hits need <= 1 hole and pT >= MinComboPT"};
  // kink chi2 cuts (multiple scattering included), tuned on the IDEA and CLD inner trackers
  Gaudi::Property<double> m_kinkMaxChi2Triplet{this, "KinkMaxChi2Triplet", 40.0,
                                               "SeedAndExtend: a hit is attached below this triplet kink chi2 (2 dof)"};
  Gaudi::Property<double> m_seedKinkMaxChi2NDF{this, "SeedKinkMaxChi2NDF", 9.0,
                                               "SeedAndExtend: seeds are kept below this kink chi2/ndf"};
  Gaudi::Property<double> m_trackKinkMaxChi2NDF{
      this, "TrackKinkMaxChi2NDF", 5.0, "SeedAndExtend: tracks are kept (and claim hits) below this kink chi2/ndf"};
  Gaudi::Property<int> m_maxCandidatesPerLayer{this, "MaxCandidatesPerLayer", 3,
                                               "SeedAndExtend: hits tried per crossed layer, nearest first"};
  Gaudi::Property<double> m_shadowHitWindow{
      this, "ShadowHitWindow", 0.1,
      "SeedAndExtend: unclaimed hits this close (cm) to an accepted track are claimed (overlapping-sensor "
      "double hits); 0 = off"};

  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Detector geometry service"};

  // ---- geometry and field ----
  dd4hep::Detector* m_detector{nullptr};
  dd4hep::OverlayedField m_field;
  const dd4hep::rec::SurfaceMap* m_surfaceMap{nullptr};
  dd4hep::rec::SurfaceMap m_mergedSurfaceMap; ///< when several detectors are combined
  std::vector<const dd4hep::rec::Surface*> m_surfaces;
  std::map<int, std::vector<const dd4hep::rec::Surface*>> m_surfacesByLayer;
  k4rt::LayerGeometry m_layerGeometry;
  std::unordered_map<int, int> m_regionOfLayerKey; ///< layer key -> readout region (0, +1, -1)
  std::unique_ptr<k4rt::SeedFinder> m_seedFinder;  ///< SeedAndExtend only
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_bitFieldCoder;
  bool m_hasRegionField{true};       ///< false: the readout has no barrel/endcap field
  bool m_useSystemInLayerKey{false}; ///< "system" in the layer key (several subdetectors, one readout)

  // ---- run statistics (Gaudi counters; mutable: operator() is const and may run in parallel) ----
  mutable Gaudi::Accumulators::StatCounter<> m_tracksPerEvent{this, "Tracks per event"};
  mutable Gaudi::Accumulators::StatCounter<> m_hitsPerTrack{this, "Hits per track"};
  mutable Gaudi::Accumulators::Counter<> m_combinationsTried{this, "Exhaustive: combinations fitted"};
  mutable Gaudi::Accumulators::Counter<> m_rejectedRoadCuts{this, "Exhaustive: rejected by the phi road cuts"};
  mutable Gaudi::Accumulators::Counter<> m_rejectedDistanceCuts{this, "Exhaustive: rejected by the distance cuts"};
  mutable Gaudi::Accumulators::Counter<> m_rejectedTanLambdaCut{this, "Exhaustive: rejected by the tanLambda cut"};
  mutable Gaudi::Accumulators::Counter<> m_rejectedChi2Cut{this, "Exhaustive: rejected by the chi2/ndf cut"};
  mutable Gaudi::Accumulators::Counter<> m_rejectedPtCuts{this, "Exhaustive: rejected by the pT cuts"};
  mutable Gaudi::Accumulators::BinomialCounter<> m_solenoidCrossingBuilt{this, "SolenoidCrossing AtIP state built"};
};

#endif // COMBINATORIAL_TRACK_FINDER_H
