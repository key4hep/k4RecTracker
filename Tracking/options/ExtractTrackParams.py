from os import getenv
from pathlib import Path

from Configurables import TrackParamExtractor
from Gaudi.Configuration import INFO, VERBOSE
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser

detModOpts = ["v02", "if1", "if2"]


class CaseInsensitiveDict(dict):
    def __getitem__(self, key):
        return super().__getitem__(key.lower())

    def __setitem__(self, key, value):
        super().__setitem__(key.lower(), value)


detModNames = CaseInsensitiveDict({el: el for el in detModOpts})

parser.add_argument(
    "--detectorModel",
    "-m",
    help="Which detector model to run reconstruction for",
    choices=detModOpts + [el.upper() for el in detModOpts],
    type=str,
    default="V02",
)
parser.add_argument(
    "--version", type=str, help="str to identify a run through the pipeline"
)
args = parser.parse_known_args()[0]


fileSuffix = ".edm4hep.root"
procName = "TrackParamExtractor"
basePath = Path(getenv("prmDir", Path.home() / "promotion"))
corePath = f"{args.version}_{detModNames[args.detectorModel]}"

iosvc = IOSvc()
iosvc.Input = str(
    (
        basePath / "code/ILDConfig/StandardConfig/production/data" / f"{corePath}_REC"
    ).with_suffix(fileSuffix)
)
iosvc.Output = str((basePath / "data" / procName / corePath).with_suffix(fileSuffix))
# iosvc.outputCommands = ["drop *", "keep SiTrackPhi"]

printer = TrackParamExtractor(procName, nStars=40)
printer.OutputLevel = VERBOSE

# the collection name with the SiTracks differs between ILC and FCC models
if "IF" in args.detectorModel:
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
