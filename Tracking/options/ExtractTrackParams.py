import os
from pathlib import Path

from commonArgParsing import add_common_args, detModNames, registry
from Configurables import TrackParamExtractor
from Gaudi.Configuration import INFO, VERBOSE
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser

ARGS = add_common_args(parser).parse_known_args()[0]
assert len(ARGS.detectorModels) == 1, (
    f"Only provide one detector model! You provided {ARGS.detectorModels}"
)
ARGS.detectorModels = ARGS.detectorModels[0]

FILE_SUFFIX = ".edm4hep.root"
PROCESSOR_NAME = "TrackParamExtractor"
BASE_PATH = Path(os.getenv("prmDir", Path.home() / "promotion"))
IN_OUT_BASE_PATH = BASE_PATH / "data" / PROCESSOR_NAME
CORE_PATH = f"{ARGS.version}_{detModNames[ARGS.detectorModels]}"

# assert that the input path exists
INPUT_PATH = (IN_OUT_BASE_PATH / "input_data" / f"{CORE_PATH}_REC").with_suffix(
    FILE_SUFFIX
)
assert INPUT_PATH.exists(), f"ERROR: The input path ({INPUT_PATH}) does not exist!"

iosvc = IOSvc()
iosvc.Input = str(INPUT_PATH)
iosvc.Output = str(
    (IN_OUT_BASE_PATH / "out_track_params" / f"{CORE_PATH}_track_params").with_suffix(
        FILE_SUFFIX
    )
)
# iosvc.outputCommands = ["drop *", "keep SiTrackPhi"]

printer = TrackParamExtractor(PROCESSOR_NAME, nStars=40)
printer.OutputLevel = VERBOSE

# the collection name with the SiTracks differs between ILC and FCC models
if registry.get(ARGS.detectorModels).at_fcc:
    printer.InputSiTracks = ["SiTracksCT"]
    SI_TRACK_COLL_NAME = "SiTracksCT"
else:
    SI_TRACK_COLL_NAME = "SiTracks"
iosvc.CollectionNames = ["ClupatraTracks", SI_TRACK_COLL_NAME]


ApplicationMgr(
    TopAlg=[printer],
    EvtSel="NONE",
    EvtMax=10,
    ExtSvc=[iosvc],
    OutputLevel=INFO,
)
