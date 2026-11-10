#ifndef COMBINATORIAL_TRACK_FINDER_H
#define COMBINATORIAL_TRACK_FINDER_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "GAUDI_VERSION.h"
#include "Gaudi/Accumulators.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ServiceHandle.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Fields.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "DDSegmentation/BitFieldCoder.h"
#include "LayerGeometry.h"
#include "SeedFinder.h"
#include <Eigen/Dense>

#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// EDM4hep includes
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

// Interface includes
#include "k4Interface/IGeoSvc.h"

#include <iomanip>
#include <unordered_map>
#include <unordered_set>

// Forward declarations
namespace dd4hep {
class OverlayedField;
namespace rec {
  class Surface;
}
} // namespace dd4hep

/** @class CombinatorialTrackFinder
 *
 *  Track finding on the planar hits of one tracking system (inner tracker, silicon
 *  wrapper, muon system, ...), with no input from other detectors.
 *
 *  How it works, per event:
 *   1. Hits are grouped by layer (barrel / endcap layers from the cellID).
 *   2. Track candidates are found either
 *        - Exhaustive: all k-hit combinations over the layers, each tested with a
 *          circle fit in xy and a line fit of z versus arc length, or
 *        - SeedAndExtend: 4-hit then 3-hit seeds on neighbouring layers, best seeds
 *          first, extended layer by layer with a helix (kink chi2 with multiple
 *          scattering; hits are claimed as tracks are accepted).
 *      Road, outlier and quality cuts reject wrong combinations.
 *   3. The best candidate is stored as a track; its hits are removed and the search
 *      repeats until no candidate is left.
 *   4. Each track stores the fitted helix as EDM4hep track states: AtFirstHit and
 *      AtLastHit, optionally AtIP (IPStateMode). For a Kalman fit of the found tracks,
 *      run k4RecTracker's GenfitTrackFitter on the output collection.
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
  virtual ~CombinatorialTrackFinder() = default;

  StatusCode initialize() override;

  std::tuple<edm4hep::TrackCollection, edm4hep::TrackerHitPlaneCollection,
             edm4hep::TrackerHitSimTrackerHitLinkCollection>
  operator()(const edm4hep::TrackerHitPlaneCollection& hits, const edm4hep::EventHeaderCollection& headers,
             const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks) const override;

  // The track search in one event (findTracks creates one and runs it; see CombinatorialTrackFinder.cpp)
  class EventSearch;

  // Pattern recognition and track building for one event
  void findTracks(const edm4hep::TrackerHitPlaneCollection* hits, edm4hep::TrackCollection& outputTracks,
                  edm4hep::TrackerHitPlaneCollection& outputHits,
                  const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks,
                  edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinks) const;

  // DEBUG printout of the tracks of one event, one line per track
  void printTrackSummary(const edm4hep::TrackCollection& trackCollection) const;

  // Detector surface of a hit
  const dd4hep::rec::Surface* findSurface(const edm4hep::TrackerHitPlane& hit) const;

private:
  // Detector surface of a cellID
  const dd4hep::rec::Surface* findSurfaceByID(uint64_t cellID) const;
  // Extract type ID from cell ID
  int getTypeID(uint64_t cellID) const;
  // Extract layer ID from cell ID
  int getLayerID(uint64_t cellID) const;
  // Layer key: barrel -> layer, +endcap -> 1000 + layer, -endcap -> -1000 - layer
  int getCompositeID(uint64_t cellID) const;
  // pT from omega, with the field at the state's reference point
  double getPT(const edm4hep::TrackState& state) const;

  // Create a track state with the given parameters
  edm4hep::TrackState createTrackState(double d0, double phi, double omega, double z0, double tanLambda,
                                       int location = edm4hep::TrackState::AtOther) const;

  // Group surfaces by detector layer
  std::map<int, std::vector<const dd4hep::rec::Surface*>> getSurfacesByLayer() const;

  // Circle fit to N >= 3 hits (weighted least squares + Gauss-Newton).
  // fitCov3x3: covariance of (x0, y0, R) in cm^2; residualsCm: per-hit residual (cm).
  // Returns false if the fit is degenerate.
  bool fitCircleNHits(const std::vector<edm4hep::TrackerHitPlane>& hits, double& x0, double& y0, double& radius,
                      double& chi2ndf, Eigen::Matrix3d& fitCov3x3, std::vector<double>& residualsCm) const;

  // ==================== PROPAGATION THROUGH THE SOLENOID ====================
  bool isInsideSolenoid(const Eigen::Vector3d& pos) const;
  bool allHitsOutsideSolenoid(const std::vector<edm4hep::TrackerHitPlane>& trackHits) const;

  // IPStateMode = SolenoidCrossing: the outer circle (centre, radius in cm) is matched
  // at the solenoid boundary to an inner helix with the inner field, which is stored as
  // AtIP. pT and tanLambda are kept across the boundary; the innermost hit (cm) selects
  // the crossing point and anchors z; clockwise = sense of rotation of the outer helix
  // seen from +z. Returns true on success.
  bool analyticalInnerPropagation(double cx_o_cm, double cy_o_cm, double R_outer_cm, bool clockwise, double pT_GeV,
                                  double tanLambda, double x_ref_cm, double y_ref_cm, double z_ref_cm,
                                  bool isEndcap, // true = endcap hit -> use z-cap crossing
                                  edm4hep::MutableTrack& finalTrack) const;

  // Where the central solenoid field ends (mm), found by following the field outward from
  // near the IP: alongZ = false gives the radius, alongZ = true the half-length
  double findSolenoidBoundary(bool alongZ) const;
  // Bz (T) at a distance (mm) along the scan direction of findSolenoidBoundary
  double fieldAlongScan(double distanceMm, bool alongZ) const;
  // true when the field there has moved more than half-way from bzCentral towards zero or beyond
  bool outsideCentralField(double distanceMm, bool alongZ, double bzCentral) const;
  // does the readout's cellID encoding have a field of this name?
  bool readoutHasField(const std::string& name) const;

  // Settings shared by both search strategies (filled in initialize)
  k4rt::CommonSearchConfig m_commonConfig;

  // Solenoid boundary
  double m_solenoidBoundaryMm{0.0}; // set in initialize() by findSolenoidBoundary(false)
  double m_solenoidHalfZmm{0.0};    // set in initialize() by findSolenoidBoundary(true)

  // Properties (options), grouped by the strategy that uses them
  // ---- options used by both strategies ----
  // readout and detector
  Gaudi::Property<std::vector<std::string>> m_detectorNames{
      this,
      "DetectorNames",
      {},
      "Names of the DetElements whose hits are tracked: one (e.g. [\"Muon-System\"]) or several "
      "that share one readout (e.g. [\"Vertex\", \"SiWrB\", \"SiWrD\"] for GlobalTrackerReadoutID); "
      "their surface maps are merged"};
  Gaudi::Property<std::string> m_encodingStringParameter{this, "EncodingStringParameterName", "GlobalTrackerReadoutID",
                                                         "Name of DD4hep parameter with the encoding string"};
  Gaudi::Property<std::string> m_layerFieldName{this, "LayerFieldName", "layer",
                                                "CellID field holding the layer index"};
  // CellID field names: the region field is 0 / +1 / -1 for barrel / +endcap / -endcap
  Gaudi::Property<std::string> m_regionFieldName{
      this, "RegionFieldName", "type",
      "CellID field separating barrel (0) from +endcap (1) and -endcap (-1). "
      "\"type\" for MuonSystemReadoutID, \"side\" for GlobalTrackerReadoutID."};
  // search strategy:
  //   Exhaustive   : all k-hit combinations over the layers (sparse events, e.g. a muon system)
  //   SeedAndExtend: 4-hit then 3-hit seeds, best first, extended along their helix (dense events)
  Gaudi::Property<std::string> m_strategy{this, "Strategy", "Exhaustive", "Exhaustive | SeedAndExtend"};
  Gaudi::Property<std::string> m_ipStateMode{
      this, "IPStateMode", "None",
      "Extra track state AtIP: None (not stored), Helix (the fitted helix at its closest "
      "approach to the origin; for tracks in a uniform field, e.g. an inner tracker) or "
      "SolenoidCrossing (helix matched at the solenoid boundary; only for tracks measured "
      "outside the solenoid, e.g. the muon system)"};
  Gaudi::Property<int> m_minTrackHits{this, "MinTrackHits", 3, "Minimum number of hits required to form a track."};
  Gaudi::Property<double> m_minComboPT{
      this, "MinComboPT", -1.0,
      "Minimum pT in GeV/c. Exhaustive: combinations below it are rejected (no fallback). "
      "SeedAndExtend: tracks with at most ShortTrackMaxHits hits must reach it (random low-pT "
      "triplets pass any kink cut). A value <= 0 disables the cut (default)."};
  Gaudi::Property<double> m_maxPairHitDist{
      this, "MaxPairHitDistance", 300.0,
      "Maximum 3D distance in cm between any pair of hits in a combination. "
      "Catches ghost combos that span widely separated detector regions. "
      "Set to -1 to disable. 300 cm accommodates barrel->endcap transition tracks."};
  Gaudi::Property<double> m_maxHitEdepKeV{
      this, "MaxHitEdepKeV", 10.0,
      "Maximum energy deposit in keV for a hit to be used in tracking. "
      "Hits above this threshold are excluded before combo enumeration. "
      "Rejects high-edep delta-rays and hadronic secondaries. Set to -1 to disable."};
  Gaudi::Property<double> m_sigmaHitDefault{
      this, "SigmaHitDefault", 0.04,
      "Default hit position resolution in cm (0.04 cm = 0.4 mm) used when hit du/dv is not set."};

  // ---- options of the Exhaustive strategy only ----
  Gaudi::Property<int> m_maxCombinatorialHits{
      this, "MaxCombinatorialHits", 7,
      "Exhaustive: Maximum number of hits considered per combination (caps C(N,k) combinatorics). "
      "Increase cautiously - C(10,7)=120 is fine, C(20,7)=77520 may be slow."};
  Gaudi::Property<double> m_nHitMaxChi2NDF{
      this, "NHitMaxChi2NDF", 10.0, "Exhaustive: maximum chi2/NDF of the N-hit circle fit to accept a combination"};
  Gaudi::Property<double> m_maxComboPT{
      this, "MaxComboPT", 250.0,
      "Exhaustive: Maximum allowed combo pT in GeV/c. Combinations whose circle-fit radius "
      "implies pT above this threshold are excluded during the combinatorial search."};
  Gaudi::Property<double> m_maxConsecDeltaPhi{
      this, "MaxConsecDeltaPhi", 0.5,
      "Exhaustive: Maximum |dphi| in radians between consecutive (adjacent-layer) hits in a combination. "
      "Rejects cross-track combos where the azimuthal step between layers is too large. "
      "Set to -1 to disable."};
  Gaudi::Property<double> m_maxComboPhiSpread{
      this, "MaxComboPhiSpread", 1.0,
      "Exhaustive: Maximum total phi spread in radians across all hits in a combination (max - min phi, "
      "relative to the innermost hit). Rejects combos spanning more than this azimuthal window. "
      "Set to -1 to disable."};
  Gaudi::Property<double> m_maxComboTanLambdaDev{
      this, "MaxComboTanLambdaDev", 0.5,
      "Exhaustive: largest deviation of the dip of a pair of consecutive hits (dz over their "
      "transverse distance) from the median of the combination (0.5 = about 27 degrees of dip "
      "angle); pure barrel or pure endcap combinations only. Set to -1 to disable."};
  Gaudi::Property<double> m_maxConsecHitDist{
      this, "MaxConsecutiveHitDistance", 100.0,
      "Exhaustive: Maximum 3D distance in cm between consecutive (inner->outer ordered) hits in a combination. "
      "Rejects combos where adjacent-layer hits are implausibly far apart (e.g. hits from "
      "different physical tracks combined across a large gap). Set to -1 to disable."};
  Gaudi::Property<double> m_hitIsolationCut{
      this, "HitIsolationCut", 2.0,
      "Exhaustive: Minimum distance in cm between a candidate hit and any other unused hit not in the "
      "combination. Set to 0 to disable. Rejects hits inside secondary shower clusters."};
  Gaudi::Property<double> m_outlierSigma{
      this, "OutlierSigma", 5.0,
      "Exhaustive: Outlier rejection threshold in units of SigmaHitDefault. "
      "Hits with residual > OutlierSigma * SigmaHitDefault are iteratively removed."};
  Gaudi::Property<int> m_maxOutlierIterations{this, "MaxOutlierIterations", 10,
                                              "Exhaustive: rounds of outlier removal per track"};
  Gaudi::Property<int> m_maxHitsPerComposite{
      this, "MaxHitsPerCompositeID", 0,
      "Exhaustive: If > 0, at most this many hits per compositeID are kept before combo enumeration, "
      "ranked by |edep - mean_single_layer_edep| (lowest deviation first). "
      "Reduces combinatorial explosion in high-occupancy layers. "
      "Set to 0 to use all hits (default)."};
  Gaudi::Property<double> m_edepScoreWeight{
      this, "EdepScoreWeight", 1.0,
      "Exhaustive: weight of the energy-deposit score when ranking combinations with the same "
      "number of hits (k = 3: tie-break on the radius; k >= 4: chi2/ndf + weight * deviation / "
      "EdepNormKeV); 0 = off"};
  Gaudi::Property<double> m_edepNormKeV{
      this, "EdepNormKeV", 2.0, "Exhaustive: energy-deposit deviation (keV) worth one chi2/ndf unit in the ranking"};
  Gaudi::Property<bool> m_doHitClustering{
      this, "DoHitClustering", false,
      "Exhaustive: Pre-cluster unused hits into spatially coherent road candidates before the combinatorial "
      "search. Uses Union-Find on adjacent compositeID layers with the same MaxConsecDeltaPhi and "
      "MaxConsecutiveHitDistance thresholds. Each cluster becomes an independent hit-pool, "
      "preventing hits from unrelated track regions from ever entering the same combo. "
      "Clusters spanning fewer than MinTrackHits distinct layers are skipped. "
      "Enable for high-occupancy events where delta-ray or multi-track hits inflate the "
      "number of combinations."};
  Gaudi::Property<bool> m_doCrowdedLayerHitSelection{this, "doCrowdedLayerHitSelection", false,
                                                     "Exhaustive: Guided crowded-layer hit selection"};
  Gaudi::Property<bool> m_requirePairCoinc{
      this, "RequirePairCoincidence", false,
      "Exhaustive: If true, exploit paired-layer geometry (doublets at layers (0,1), (2,3), (4,5)) "
      "as a hit-quality filter. Real muons traverse both layers of a doublet within a "
      "tight (dphi, deta) window; delta-rays/secondaries usually fire only one. "
      "Two-pass: pass 1 enumerates combos using only hits with a coincident partner in "
      "the paired layer; pass 2 (if pass 1 finds nothing) uses all hits. Either way, "
      "candidate combos are required to have >=2 complete doublets - exactly 3 is best, "
      "2 is acceptable, <2 rejected."};
  Gaudi::Property<double> m_pairCoincMaxDist{
      this, "PairCoincidenceMaxDist", 25.0,
      "Exhaustive: Max 3D distance (cm) between two hits treated as a paired-layer doublet. "
      "25 cm suits the muon system (paired layers 12-24 cm apart); scale it down "
      "for inner-tracker readouts where layers are far closer."};
  Gaudi::Property<double> m_neighbourTrackMaxDist{
      this, "NeighbourTrackMaxDist", 5.0,
      "Exhaustive: If the best-combo candidate has any hit within this distance (cm) of a hit from "
      "an already-accepted track, apply stricter quality requirements defined by "
      "NeighbourTrackMinHits and NeighbourTrackMaxChi2NDF. "
      "Rejects delta-ray shadow tracks that form next to a real track. Set to 0 to disable."};
  Gaudi::Property<int> m_neighbourTrackMinHits{
      this, "NeighbourTrackMinHits", 4,
      "Exhaustive: Minimum hit count required for a combo flagged as a neighbour track "
      "(any hit lies within NeighbourTrackMaxDist of an accepted track hit). "
      "Stronger than MinTrackHits, to reject 3-hit shadow fakes."};
  Gaudi::Property<double> m_neighbourTrackMaxChi2NDF{
      this, "NeighbourTrackMaxChi2NDF", 3.0,
      "Exhaustive: largest chi2/ndf of a combination next to an already accepted track "
      "(applied with NeighbourTrackMinHits)"};

  // ---- options of the SeedAndExtend strategy only ----
  Gaudi::Property<std::string> m_seedLayerPolicy{
      this, "SeedLayerPolicy", "innermost",
      "SeedAndExtend: where seeds are made, and so which way they are extended: innermost (the "
      "first layers a ray from the IP crosses; extended outward), outermost (the last layers; "
      "extended inward) or consecutive (any run of layers, for displaced tracks; extended both ways)"};
  // Default tuned on the IDEA and CLD inner trackers: in b-bbar events without these checks most
  // 3- and 4-hit tracks are fakes, longer ones almost never are. A larger value removes more
  // duplicates but also loses real tracks and gives more fakes.
  Gaudi::Property<int> m_shortTrackMaxHits{
      this, "ShortTrackMaxHits", 4,
      "SeedAndExtend: tracks with at most this many hits must have at most 1 hole and pT >= MinComboPT"};
  // Kink chi2 cuts of SeedAndExtend (multiple scattering included, calibrated to chi2/ndf ~ 1).
  // Defaults tuned on the IDEA and CLD inner trackers (signal muons and b-bbar / c-cbar events).
  Gaudi::Property<double> m_kinkMaxChi2Triplet{
      this, "KinkMaxChi2Triplet", 40.0,
      "SeedAndExtend: a hit is attached only if the kink chi2 (2 dof) of its triplet is below this"};
  Gaudi::Property<double> m_seedKinkMaxChi2NDF{this, "SeedKinkMaxChi2NDF", 9.0,
                                               "SeedAndExtend: seeds are kept below this kink chi2/ndf"};
  Gaudi::Property<double> m_trackKinkMaxChi2NDF{
      this, "TrackKinkMaxChi2NDF", 5.0,
      "SeedAndExtend: a finished track is kept, and claims its hits, below this kink chi2/ndf"};
  Gaudi::Property<int> m_maxCandPerLayer{this, "MaxCandidatesPerLayer", 3,
                                         "SeedAndExtend: hits tried per crossed layer, nearest first"};
  Gaudi::Property<double> m_shadowHitWindow{
      this, "ShadowHitWindow", 0.1,
      "SeedAndExtend: unclaimed hits within this distance (cm) of an accepted track's trajectory are "
      "claimed (overlapping-sensor double hits), so they cannot seed a duplicate; 0 = off"};

  // Services
  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Detector geometry service"};

  // Geometry and field
  dd4hep::Detector* m_detector{nullptr};
  dd4hep::OverlayedField m_field;
  const dd4hep::rec::SurfaceMap* m_surfaceMap{nullptr};
  std::vector<const dd4hep::rec::Surface*> m_surfaces;
  std::map<int, std::vector<const dd4hep::rec::Surface*>> m_surfacesByLayer;
  k4rt::LayerGeometry m_layerGeo;                 // per-layer radius/z extents, built from m_surfacesByLayer
  std::unordered_map<int, int> m_cidRegion;       // compositeID -> readout region (0, +1, -1); see k4rt::regionOfCid
  bool m_useSeedExtend{false};                    // Strategy == "SeedAndExtend"
  std::unique_ptr<k4rt::SeedFinder> m_seedFinder; // built in initialize() when m_useSeedExtend
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_bitFieldCoder; // cellID decoding
  bool m_hasRegionField{true};          // false when the readout has no barrel/endcap field at all
  bool m_useSystemInCompositeID{false}; // fold "system" into the layer key (multi-subdetector readouts)
  // merged surface map when several detectors are combined (DetectorNames)
  dd4hep::rec::SurfaceMap m_mergedSurfaceMap;

  // Run statistics: thread-safe Gaudi counters, printed by Gaudi at the end of the job
  // (mutable because operator() is const and may run on several events at once)
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

DECLARE_COMPONENT(CombinatorialTrackFinder)
#endif // COMBINATORIAL_TRACK_FINDER_H