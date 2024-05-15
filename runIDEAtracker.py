import os

from Gaudi.Configuration import *

# Loading the input SIM file
from Configurables import k4DataSvc, PodioInput

evtsvc = k4DataSvc("EventDataSvc")
evtsvc.input = "output_IDEA_DIGI.root"
inp = PodioInput("InputReader")


# pattern recognition over digitized hits
from Configurables import GenFitter

dch_tracking = GenFitter(
    "GenFitter",
    modelPath="/afs/cern.ch/work/m/mgarciam/private/k4RecTracker_dev_0/Tracking/model_multivector_1_input.onnx",
    inputHits_CDC="CDCHDigis",
    inputHits_VTXD="VTXDDigis",
    inputHits_VTXIB="VTXIBDigis",
    inputHits_VTXOB="VTXOBDigis",
    outputTracks="CDCHTracks",
    OutputLevel=INFO,
)

################ Output
from Configurables import PodioOutput

out = PodioOutput("out", OutputLevel=INFO)
out.outputCommands = ["keep *"]

out.filename = "output_tracking.root"

# CPU information
from Configurables import AuditorSvc, ChronoAuditor

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
out.AuditExecute = True

from Configurables import EventCounter

event_counter = EventCounter("event_counter")
event_counter.Frequency = 1

from Configurables import ApplicationMgr

ApplicationMgr(
    TopAlg=[
        inp,
        event_counter,
        dch_tracking,
        out,
    ],
    EvtSel="NONE",
    EvtMax=1,
    ExtSvc=[evtsvc, audsvc],
    StopOnSignal=True,
)
