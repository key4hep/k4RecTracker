from pathlib import Path

from Configurables import TrackD0Printer
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc

iosvc = IOSvc()
iosvc.Input = str(
    Path.home()
    / "promotion/code/ILDConfig/StandardConfig/production/data/test_tracking_3_detmods_V02_REC.edm4hep.root"
)

iosvc.CollectionNames = ["SiTracks", "ClupatraTracks"]

printer = TrackD0Printer("TrackD0Printer", nStars=40)

ApplicationMgr(
    TopAlg=[printer],
    EvtSel="NONE",
    EvtMax=10,
    ExtSvc=[iosvc],
    OutputLevel=INFO,
)
