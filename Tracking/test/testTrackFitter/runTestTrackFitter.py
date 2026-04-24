import os

from Gaudi.Configuration import INFO, WARNING, DEBUG

from Configurables import AuditorSvc, ChronoAuditor
from Configurables import EventDataSvc


################ parser
from k4FWCore.parseArgs import parser

parser.add_argument("--inputFile", default="output_tracks.root", help="InputFile")
parser.add_argument("--outputFile", default="testFitter.root", help="OutputFile")
parser.add_argument(
    "--Beta_init",
    default=100,
    help="Initial annealing parameter (if FitterType == 'DAF') [https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf]",
)
parser.add_argument(
    "--Beta_final",
    default=0.05,
    help="Final annealing parameter (if FitterType == 'DAF') [https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf]",
)
parser.add_argument(
    "--Beta_steps",
    default=15,
    help="Number of annealing steps (if FitterType == 'DAF') [https://indico.cern.ch/event/258092/papers/1588579/files/4253-genfit.pdf]",
)
args = parser.parse_args()

################ input & output
from k4FWCore import IOSvc

io_svc = IOSvc("IOSvc")
io_svc.Input = args.inputFile
io_svc.Output = args.outputFile

# pattern recognition over digitized hits
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

trackFitter.BetaInit = args.Beta_init
trackFitter.BetaFinal = args.Beta_final
trackFitter.BetaSteps = args.Beta_steps

trackFitter.InitializationType = 1
trackFitter.SkipTrackOrdering = False
trackFitter.SkipUnmatchedTracks = False

trackFitter.FilterTrackHits = True

# ################ Detector geometry
from Configurables import GeoSvc

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [os.environ["K4GEO"] + "/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

################ Application
from k4FWCore import ApplicationMgr

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
ApplicationMgr(
    TopAlg=[trackFitter],
    EvtSel="NONE",
    ExtSvc=[geoservice, EventDataSvc("EventDataSvc"), audsvc],
    StopOnSignal=True,
    EvtMax=-1,
    OutputLevel=INFO,
)
