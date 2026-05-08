from Gaudi.Configuration import *

from Configurables import GeoSvc, EventDataSvc
from k4FWCore import ApplicationMgr, IOSvc

geoservice = GeoSvc("GeoSvc")
geoservice.detectors.append("/path/to/k4geo/FCCee/compact/CLD_ARC/arc_full_v0.xml")

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

from Configurables import EventCounter

event_counter = EventCounter("event_counter")
event_counter.Frequency = 1

ApplicationMgr(
    TopAlg=[event_counter, arc_digitizer],
    EvtSel="NONE",
    EvtMax=10,
    ExtSvc=[geoservice, EventDataSvc("EventDataSvc")],
)
