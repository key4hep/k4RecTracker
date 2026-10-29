# Steering (pipeline A): digitized hits -> GGTF track finding -> Genfit refit.
#
# The GGTFTrackFinder configuration is a copy of
# /afs/cern.ch/work/a/adevita/public/bib-studies/IDEAsim/runTestTrackFinder.py
# (model path set to smallModelForTest.onnx), extended with the GenfitTrackFitter
# step following the example in k4RecTracker/Tracking/test/testTrackFitter.
#
# Usage: k4run runIDEA_GGTF_Genfit.py --inputFile digi.edm4hep.root --outputFile tracks_ggtf_genfit.edm4hep.root

import os

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc, UniqueIDGenSvc, GeoSvc
from k4FWCore import IOSvc, ApplicationMgr
from k4FWCore.parseArgs import parser

################## Parser
parser.add_argument("--inputFile", default="digi.edm4hep.root", help="InputFile")
parser.add_argument("--outputFile", default="tracks_ggtf_genfit.edm4hep.root", help="OutputFile")
parser.add_argument(
    "--modelPath",
    default="/afs/cern.ch/work/a/adevita/public/trainingFolder/GGTF_trackFinder_training/conversionModel/smallModelForTest.onnx",
    help="model path for the track finder",
)
parser.add_argument("--tbeta", default=0.7, help="tbeta clustering parameter")
parser.add_argument("--td", default=0.1, help="td clustering parameter")
args = parser.parse_args()

################## InputOutput
svc = IOSvc("IOSvc")
svc.Input = args.inputFile
svc.Output = args.outputFile

################ Detector geometry (needed by GenfitTrackFitter)
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
geoservice.detectors = [
    os.path.join(path_to_detector, "FCCee/IDEA/compact/IDEA_o1_v04/IDEA_o1_v04.xml")
]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

############### Track Finder
from Configurables import GGTFTrackFinder

GGTF = GGTFTrackFinder(
    "GGTFTrackFinder",
    InputPlanarHitCollections=["VTXBDigis", "VTXDDigis", "SiWrBDigis", "SiWrDDigis"],
    InputWireHitCollections=["DCH_DigiCollection"],
    OutputTracksGGTF=["GGTFTracks"],
    ModelPath=args.modelPath,
    EmbeddingCoordinates=4,
    Tbeta=args.tbeta,
    Td=args.td,
    OutputLevel=INFO,
)

############### Track Fitter
from Configurables import GenfitTrackFitter

trackFitter = GenfitTrackFitter(
    "GenfitTrackFitter",
    InputTracks=["GGTFTracks"],
    OutputFittedTracks=["FittedTracks"],
    OutputFittedTracksWithFilteredHits=["FittedTracksWithFilteredHits"],
    OutputFittedHits=["FittedHits"],
    OutputLevel=INFO,
)

trackFitter.RunSingleEvaluation = True
trackFitter.UseBrems = False

trackFitter.BetaInit = 100
trackFitter.BetaFinal = 0.05
trackFitter.BetaSteps = 15

trackFitter.InitializationType = 1
trackFitter.SkipTrackOrdering = False
trackFitter.ListOfTypesToSkip = []

trackFitter.FilterTrackHits = True

############### Application Manager
mgr = ApplicationMgr(
    TopAlg=[GGTF, trackFitter],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice, EventDataSvc("EventDataSvc"), UniqueIDGenSvc("uidSvc")],
    OutputLevel=INFO,
)
