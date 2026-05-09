import os
from Gaudi.Configuration import *

from Configurables import GeoSvc, EventDataSvc
from k4FWCore import ApplicationMgr, IOSvc

geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("LCGEO")
detectors_to_use = [
    "FCCee/CLD/compact/CLD_o3_v01/CLD_o3_v01.xml",
]
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]

iosvc = IOSvc()
iosvc.Input = vars().get("input", "data/arcsim_kaon+_edm4hep.root")
iosvc.CollectionNames = ["ARC_HITS"]
iosvc.Output = vars().get("output", "digi.root")
iosvc.outputCommands = ["keep *"]

from Configurables import ARCdigitizer

arc_digitizer = ARCdigitizer(
    "ARCdigitizer", inputSimHits="ARC_HITS", outputDigiHits="ARC_DIGI_HITS"
)

# CPU information
from Configurables import AuditorSvc, ChronoAuditor

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
arc_digitizer.AuditExecute = True

ApplicationMgr(
    TopAlg=[arc_digitizer],
    EvtSel="NONE",
    EvtMax=10,
    ExtSvc=[geoservice, EventDataSvc("EventDataSvc")],
)
