# Steering (pipeline B): digitized hits -> conformal tracking (vertex + silicon
# wrapper hits only) -> RefitFinal (DDKalTest refit with energy loss and MS).
#
# The ConformalTracking configuration follows
# k4Reco/ConformalTracking/options/runIDEAConformalTracking.py (CLD-tuned steps
# with IDEA collections); the RefitFinal one follows
# k4Reco/Tracking/options/runRefitFinal.py.
#
# Usage: k4run runIDEA_ConformalTracking_Refit.py --inputFile digi.edm4hep.root --outputFile tracks_ct_refit.edm4hep.root

import os

from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser
from Configurables import EventDataSvc, UniqueIDGenSvc, RndmGenSvc, GeoSvc
from Configurables import ConformalTracking, RefitFinal
from Configurables import RootHistSvc
from Configurables import Gaudi__Histograming__Sink__Root as RootHistoSink

from conformal_tracking_utils import configure_conformal_tracking_steps

################## Parser
parser.add_argument(
    "--inputFile", default="digi.edm4hep.root", help="InputFile (digitized, IDEA_o1_v04)"
)
parser.add_argument("--outputFile", default="tracks_ct_refit.edm4hep.root", help="OutputFile")
args = parser.parse_args()

################## InputOutput
svc = IOSvc("IOSvc")
svc.Input = args.inputFile
svc.Output = args.outputFile

################## Detector geometry
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
geoservice.detectors = [
    os.path.join(path_to_detector, "FCCee/IDEA/compact/IDEA_o1_v04/IDEA_o1_v04.xml")
]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

############### Conformal Tracking
tracking = ConformalTracking()
tracking.TrackerHitCollectionNames = [
    "VTXBDigis",
    "VTXDDigis",
    "SiWrBDigis",
    "SiWrDDigis",
]
tracking.RelationsNames = [
    "VTXBSimDigiLinks",
    "VTXDSimDigiLinks",
    "SiWrBSimDigiLinks",
    "SiWrDSimDigiLinks",
]
tracking.MCParticleCollectionName = ["MCParticles"]
tracking.SiTrackCollectionName = "SiTracksCT"

tracking.MainTrackerHitCollectionNames = ["SiWrBDigis", "SiWrDDigis"]
tracking.VertexBarrelHitCollectionNames = ["VTXBDigis"]
tracking.VertexEndcapHitCollectionNames = ["VTXDDigis"]

tracking.DebugPlots = False
tracking.DebugTiming = False
tracking.MaxHitInvertedFit = 0
tracking.MinClustersOnTrackAfterFit = 3
tracking.RetryTooManyTracks = False
tracking.SortTreeResults = True
tracking.ThetaRange = 0.05
tracking.TooManyTracks = 100000
tracking.trackPurity = 0.7

CT_MAX_DIST = 0.05

# The keys (VXDBarrel, VXDEndcap...) are simply names and are not passed to ConformalTracking
parameters = {
    "VXDBarrel": {
        "collections": ["VTXBDigis"],
        "params": {
            "MaxCellAngle": 0.01,
            "MaxCellAngleRZ": 0.01,
            "Chi2Cut": 100,
            "MinClustersOnTrack": 4,
            "MaxDistance": CT_MAX_DIST,
            "SlopeZRange": 10.0,
            "HighPTCut": 10.0,
        },
        "flags": ["HighPTFit", "VertexToTracker"],
        "functions": ["CombineCollections", "BuildNewTracks"],
    },
    "VXDEndcap": {
        "collections": ["VTXDDigis"],
        "params": {
            "MaxCellAngle": 0.01,
            "MaxCellAngleRZ": 0.01,
            "Chi2Cut": 100,
            "MinClustersOnTrack": 4,
            "MaxDistance": CT_MAX_DIST,
            "SlopeZRange": 10.0,
            "HighPTCut": 10.0,
        },
        "flags": ["HighPTFit", "VertexToTracker"],
        "functions": ["CombineCollections", "ExtendTracks"],
    },
    "LowerCellAngle1": {
        "collections": ["VTXBDigis", "VTXDDigis"],
        "params": {
            "MaxCellAngle": 0.05,
            "MaxCellAngleRZ": 0.05,
            "Chi2Cut": 100,
            "MinClustersOnTrack": 4,
            "MaxDistance": CT_MAX_DIST,
            "SlopeZRange": 10.0,
            "HighPTCut": 10.0,
        },
        "flags": ["HighPTFit", "VertexToTracker", "RadialSearch"],
        "functions": ["CombineCollections", "BuildNewTracks"],
    },
    "LowerCellAngle2": {
        "collections": [],
        "params": {
            "MaxCellAngle": 0.1,
            "MaxCellAngleRZ": 0.1,
            "Chi2Cut": 2000,
            "MinClustersOnTrack": 4,
            "MaxDistance": CT_MAX_DIST,
            "SlopeZRange": 10.0,
            "HighPTCut": 10.0,
        },
        "flags": ["HighPTFit", "VertexToTracker", "RadialSearch"],
        "functions": ["BuildNewTracks", "SortTracks"],
    },
    "SiWrapper": {
        "collections": ["SiWrBDigis", "SiWrDDigis"],
        "params": {
            "MaxCellAngle": 0.1,
            "MaxCellAngleRZ": 0.1,
            "Chi2Cut": 2000,
            "MinClustersOnTrack": 4,
            "MaxDistance": CT_MAX_DIST,
            "SlopeZRange": 10.0,
            "HighPTCut": 1.0,
        },
        "flags": ["HighPTFit", "VertexToTracker", "RadialSearch"],
        "functions": ["CombineCollections", "ExtendTracks"],
    },
    "Displaced": {
        "collections": ["VTXBDigis", "VTXDDigis", "SiWrBDigis", "SiWrDDigis"],
        "params": {
            "MaxCellAngle": 0.1,
            "MaxCellAngleRZ": 0.1,
            "Chi2Cut": 1000,
            "MinClustersOnTrack": 5,
            "MaxDistance": 0.015,
            "SlopeZRange": 10.0,
            "HighPTCut": 10.0,
        },
        "flags": ["OnlyZSchi2cut", "RadialSearch"],
        "functions": ["CombineCollections", "BuildNewTracks"],
    },
}

configure_conformal_tracking_steps(tracking, parameters)

############### Refit of the conformal tracks (DDKalTest)
refit = RefitFinal(
    "RefitFinal",
    InputTrackCollectionName="SiTracksCT",
    InputRelationCollectionName=[],
    OutputTrackCollectionName="SiTracksCT_Refitted",
    OutputRelationCollectionName="SiTracksCT_Refitted_Relation",
    EnergyLossOn=True,
    Max_Chi2_Incr=1.79769e30,
    MinClustersOnTrackAfterFit=3,
    MultipleScatteringOn=True,
    ReferencePoint=-1,
    SmoothOn=False,
    extrapolateForward=True,
)

############### Application Manager
hps = RootHistSvc("HistogramPersistencySvc")
root_hist_svc = RootHistoSink("RootHistoSink")
root_hist_svc.FileName = "conformal_tracking_idea_hist.root"

ApplicationMgr(
    TopAlg=[tracking, refit],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[
        geoservice,
        EventDataSvc("EventDataSvc"),
        UniqueIDGenSvc("uidSvc"),
        RndmGenSvc(),
        root_hist_svc,
    ],
    OutputLevel=INFO,
)
