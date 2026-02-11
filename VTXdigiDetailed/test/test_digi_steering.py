import os
from Gaudi.Configuration import INFO, WARNING, DEBUG

from Gaudi.Configurables import (
    EventDataSvc,
    MarlinProcessorWrapper,
    GeoSvc,
    TrackingCellIDEncodingSvc,
)
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser
from py_utils import SequenceLoader
from k4MarlinWrapper.io_helpers import IOHandlerHelper

parser_group = parser.add_argument_group("CLDReconstruction.py custom options")
# Need the dummy input such that the IOHandlerHelper.add_reader call below does not crash when called with --help
parser_group.add_argument(
    "--inputFiles",
    action="store",
    nargs="+",
    metavar=("file1", "file2"),
    help="One or multiple input files",
    default=["dummy_input.edm4hep.root"],
)
parser_group.add_argument(
    "--outputBasename", help="Basename of the output file(s)", default="output"
)
parser_group.add_argument(
    "--compactFile",
    help="Compact detector file to use",
    type=str,
    default=os.environ["K4GEO"] + "/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml",
)
tracking_group = parser_group.add_mutually_exclusive_group()
reco_args = parser.parse_known_args()[0]

evtsvc = EventDataSvc("EventDataSvc")
iosvc = IOSvc()

svcList = [evtsvc, iosvc]
algList = []

CONFIG = {
    "CalorimeterIntegrationTimeWindow": "10ns",
    "CalorimeterIntegrationTimeWindowChoices": ["10ns", "400ns"],
    "Overlay": "False",
    "OverlayChoices": ["False", "91GeV", "365GeV"],
    "VertexUnconstrained": "OFF",
    "VertexUnconstrainedChoices": ["ON", "OFF"],
    "OutputMode": "EDM4Hep",
    "OutputModeChoices": ["LCIO", "EDM4hep"],  # , "both"] FIXME: both is not implemented yet
}

SRC_PATH = os.environ.get("SRC_TEST", ".")

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [reco_args.compactFile]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False
svcList.append(geoservice)

cellIDSvc = TrackingCellIDEncodingSvc("CellIDSvc")
cellIDSvc.EncodingStringParameterName = "GlobalTrackerReadoutID"
cellIDSvc.GeoSvcName = geoservice.name()
cellIDSvc.OutputLevel = INFO
svcList.append(cellIDSvc)

if len(geoservice.detectors) > 1:
    # we are making assumptions for reconstruction parameters based on the detector option, so we limit the possibilities
    raise RuntimeError(
        "Too many XML files for the detector path, please only specify the main file!"
    )

# from https://github.com/HEP-FCC/FCCeePhysicsPerformance/blob/d6ecee2c2c3ed5d76db55a3ae18ced349b2b914a/General/README.md?plain=1#L457-L467
# for december 2022

sequenceLoader = SequenceLoader(
    algList,
    # global_vars can be used in sequence-loaded modules without explicit import
    global_vars={
        "CONFIG": CONFIG,
        "geoservice": geoservice,
        "reco_args": reco_args,
    },
)

io_handler = IOHandlerHelper(algList, iosvc)
io_handler.add_reader(reco_args.inputFiles)

MyAIDAProcessor = MarlinProcessorWrapper("MyAIDAProcessor")
MyAIDAProcessor.OutputLevel = WARNING
MyAIDAProcessor.ProcessorType = "AIDAProcessor"
MyAIDAProcessor.Parameters = {
    "Compress": ["1"],
    "FileName": [f"{reco_args.outputBasename}_aida"],
    "FileType": ["root"],
}

EventNumber = MarlinProcessorWrapper("EventNumber")
EventNumber.OutputLevel = WARNING
EventNumber.ProcessorType = "Statusmonitor"
EventNumber.Parameters = {"HowOften": ["1"]}

# setup AIDA histogramming and add eventual background overlay
algList.append(MyAIDAProcessor)
sequenceLoader.load(os.path.join(SRC_PATH, "Overlay"))
# tracker hit digitisation
sequenceLoader.load(os.path.join(SRC_PATH, "TrackingDigi"))

# event number processor, down here to attach the conversion back to edm4hep to it
algList.append(EventNumber)


# Output file writing, in EDM4hep
io_handler.add_edm4hep_writer(f"{reco_args.outputBasename}_REC.edm4hep.root", ["keep *"])

# We need to attach all the necessary converters
io_handler.finalize_converters()

ApplicationMgr(
    TopAlg=algList,
    EvtSel="NONE",
    EvtMax=3,  # Overridden by the --num-events switch to k4run
    ExtSvc=svcList,
    OutputLevel=WARNING,
)
