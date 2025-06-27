from os import getenv
from pathlib import Path

from Configurables import TrackParamExtractor
from Gaudi.Configuration import INFO, VERBOSE
from k4FWCore import ApplicationMgr, IOSvc

fileSuffix = ".edm4hep.root"
versionName = "v0"
procName = "TrackParamExtractor"
basePath = Path(getenv("prmDir", Path.home() / "promotion"))

iosvc = IOSvc()
iosvc.Input = str(
    (
        basePath
        / "code/ILDConfig/StandardConfig/production/data/test_tracking_3_detmods_V02_REC"
    ).with_suffix(fileSuffix)
)
iosvc.Output = str((basePath / "data" / procName / versionName).with_suffix(fileSuffix))
iosvc.CollectionNames = ["SiTracks", "ClupatraTracks"]
# iosvc.outputCommands = ["drop *", "keep SiTrackPhi"]

printer = TrackParamExtractor(procName, nStars=40)
printer.OutputLevel = VERBOSE

ApplicationMgr(
    TopAlg=[printer],
    EvtSel="NONE",
    EvtMax=10,
    ExtSvc=[iosvc],
    OutputLevel=INFO,
)
