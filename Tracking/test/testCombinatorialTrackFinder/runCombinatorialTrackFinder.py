# CombinatorialTrackFinder on the CLD silicon tracker: vertex detector, inner and outer
# trackers (barrel + disks), which share the GlobalTrackerReadoutID readout.
# k4FWCore's CollectionMerger first combines their hit and reco->sim link collections, so a track
# crossing several subdetectors is found as one track. Settings: SeedAndExtend, with the
# values of the inner-tracker studies (the kink chi2 cuts are the defaults, tuned on them).
#
# usage:
#   k4run runCombinatorialTrackFinder.py --inputFile digi.root --outputFile tracks.root
#   k4run runCombinatorialTrackFinder.py --inputFile sim.root --outputFile tracks.root --digitise
#     --digitise: the input holds simulated hits (ddsim output); digitise them first with the
#     CLD reconstruction settings (this is what test_CombinatorialTrackFinder.sh does)
import os
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser
from Configurables import EventDataSvc, GeoSvc, Gaudi__Monitoring__MessageSvcSink
from Configurables import CollectionMerger, CombinatorialTrackFinder

parser.add_argument(
    "--inputFile", required=True, help="digitised (or, with --digitise, simulated) CLD events"
)
parser.add_argument("--outputFile", required=True)
parser.add_argument(
    "--digitise", action="store_true", help="digitise the simulated tracker hits first"
)
args = parser.parse_args()

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
    os.path.join(os.environ["K4GEO"], "FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml")
]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

algorithms = []
# MessageSvcSink prints the run statistics (counters) of the algorithms at the end of the job
services = [EventDataSvc("EventDataSvc"), Gaudi__Monitoring__MessageSvcSink()]

# -- optional digitisation: Gaussian smearing on the sensor planes (CLDConfig TrackingDigi.py) --
# (subdetector, resolution u per layer, resolution v per layer (mm), simulated hits, links, digitised hits)
digitisers = [
    (
        "Vertex",
        [0.003] * 6,
        [0.003] * 6,
        "VertexBarrelCollection",
        "VXDTrackerHitRelations",
        "VXDTrackerHits",
    ),
    (
        "Vertex",
        [0.003] * 6,
        [0.003] * 6,
        "VertexEndcapCollection",
        "VXDEndcapTrackerHitRelations",
        "VXDEndcapTrackerHits",
    ),
    (
        "InnerTrackers",
        [0.007],
        [0.09],
        "InnerTrackerBarrelCollection",
        "InnerTrackerBarrelHitsRelations",
        "ITrackerHits",
    ),
    (
        "InnerTrackers",
        [0.005] + [0.007] * 6,
        [0.005] + [0.09] * 6,
        "InnerTrackerEndcapCollection",
        "InnerTrackerEndcapHitsRelations",
        "ITrackerEndcapHits",
    ),
    (
        "OuterTrackers",
        [0.007] * 3,
        [0.09] * 3,
        "OuterTrackerBarrelCollection",
        "OuterTrackerBarrelHitsRelations",
        "OTrackerHits",
    ),
    (
        "OuterTrackers",
        [0.007] * 5,
        [0.09] * 5,
        "OuterTrackerEndcapCollection",
        "OuterTrackerEndcapHitsRelations",
        "OTrackerEndcapHits",
    ),
]
if args.digitise:
    from Configurables import DDPlanarDigi, RndmGenSvc, UniqueIDGenSvc

    # the digitisers draw random numbers per event
    services += [UniqueIDGenSvc("uidSvc"), RndmGenSvc()]
    for subdetector, resolution_u, resolution_v, sim_hits, links, digis in digitisers:
        digitiser = DDPlanarDigi(digis + "Digitiser")
        digitiser.SubDetectorName = subdetector
        digitiser.IsStrip = False
        digitiser.ResolutionU = resolution_u
        digitiser.ResolutionV = resolution_v
        digitiser.SimTrackHitCollectionName = [sim_hits]
        digitiser.SimTrkHitRelCollection = [links]
        digitiser.TrackerHitCollectionName = [digis]
        algorithms.append(digitiser)

# -- all tracker hits in one collection --
merge_hits = CollectionMerger("MergeInnerTrackerHits")
merge_hits.InputCollections = [d[5] for d in digitisers]
merge_hits.OutputCollection = "InnerTrackerHits"
algorithms.append(merge_hits)
merge_links = CollectionMerger("MergeInnerTrackerHitSimLinks")
merge_links.InputCollections = [d[4] for d in digitisers]
merge_links.OutputCollection = "InnerTrackerHitSimLinks"
algorithms.append(merge_links)

t = CombinatorialTrackFinder("InnerTrackerTracking")
t.OutputLevel = INFO
# -- readout --
t.EncodingStringParameterName = "GlobalTrackerReadoutID"
t.RegionFieldName = "side"  # barrel / endcap field of this readout
# the subdetectors whose sensor surfaces are used
t.DetectorNames = ["Vertex", "InnerTrackers", "OuterTrackers"]

# -- stored track states: the fitted helix AtFirstHit and AtLastHit, plus AtIP --
t.IPStateMode = "Helix"  # uniform field: the same helix at its closest approach to the origin

# -- SeedAndExtend: seeds of 4 then 3 hits, best first, extended along their helix --
t.Strategy = "SeedAndExtend"
# seeds on any run of layers (for displaced tracks), extended both ways
t.SeedLayerPolicy = "consecutive"
# largest distance between consecutive seed hits, cm (CLD: r < 2.1 m, |z| < 2.3 m)
t.MaxPairHitDistance = 100.0
# cm; the per-hit du/dv of the digitiser are used when present
t.SigmaHitDefault = 0.0003
# no energy-deposit filter: silicon hits deposit more than the muon-system default of 10 keV
t.MaxHitEdepKeV = 0
t.MinTrackHits = 3
t.ShortTrackMaxHits = 4  # tracks with <= 4 hits: hole check and minimum pT
t.MinComboPT = 1.0  # GeV, minimum pT of short tracks (random low-pT triplets)

# -- collections --
t.InputHitCollection = ["InnerTrackerHits"]
t.InputRecoSimLinkCollection = ["InnerTrackerHitSimLinks"]
t.OutputTrackCollection = ["InnerTrackerTracks"]
t.OutputHitCollection = ["InnerTrackerTrackHits"]
t.OutputRecoSimLinkCollection = ["InnerTrackerTrackHitSimLinks"]
algorithms.append(t)

iosvc = IOSvc()
iosvc.Input = args.inputFile
iosvc.Output = args.outputFile

ApplicationMgr(TopAlg=algorithms, EvtSel="NONE", EvtMax=-1, ExtSvc=services, OutputLevel=INFO)
