import os
import sys
from pathlib import Path

from Configurables import TrackParamExtractor
from Gaudi.Configuration import INFO, VERBOSE
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser

sys.path.append(os.getenv("trckOptDir"))
from commonArgParsing import add_common_args, detModNames

args = add_common_args(parser).parse_known_args()[0]
assert len(args.detectorModels) == 1, (
    f"Only provide one detector model! You provided {args.detectorModels}"
)
args.detectorModels = args.detectorModels[0]

fileSuffix = ".edm4hep.root"
procName = "TrackParamExtractor"
basePath = Path(os.getenv("prmDir", Path.home() / "promotion"))
corePath = f"{args.version}_{detModNames[args.detectorModels]}"

# assert that the input path exists
inputPath = (
    basePath / "code/ILDConfig/StandardConfig/production/data" / f"{corePath}_REC"
).with_suffix(fileSuffix)
assert inputPath.exists(), f"ERROR: The input path ({inputPath}) does not exist!"

iosvc = IOSvc()
iosvc.Input = str(inputPath)
iosvc.Output = str((basePath / "data" / procName / corePath).with_suffix(fileSuffix))
# iosvc.outputCommands = ["drop *", "keep SiTrackPhi"]

printer = TrackParamExtractor(procName, nStars=40)
printer.OutputLevel = VERBOSE

# the collection name with the SiTracks differs between ILC and FCC models
if "IF" in args.detectorModels:
    printer.InputSiTracks = ["SiTracksCT"]
    siTrackCollName = "SiTracksCT"
else:
    siTrackCollName = "SiTracks"
iosvc.CollectionNames = ["ClupatraTracks", siTrackCollName]


ApplicationMgr(
    TopAlg=[printer],
    EvtSel="NONE",
    EvtMax=10,
    ExtSvc=[iosvc],
    OutputLevel=INFO,
)
