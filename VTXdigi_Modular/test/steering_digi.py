from Gaudi.Configuration import *
from k4FWCore import IOSvc

ioSvc = IOSvc("IOSvc")
ioSvc.Input = "SET_FROM_COMMAND_LINE"
ioSvc.Output = "SET_FROM_COMMAND_LINE"
ioSvc.outputCommands = ["keep *"]

from Configurables import GeoSvc

geoSvc = GeoSvc("GeoSvc", OutputLevel=WARNING, detectors=["SET_FROM_COMMAND_LINE"])

from Configurables import VTXdigi_Modular

vtxb_digitizer = VTXdigi_Modular(
    "VTXBdigi_inner",
    SimTrackHitCollectionName=["VertexBarrelCollection"],
    HeaderName=["EventHeader"],
    SimTrkHitRelationsCollection=["VTXBSimDigiLinks"],
    TrackerHitCollectionName=["VTXBDigis"],
    SubDetectorName="VertexBarrel",
    Clusterize=True,
    Layers=[0, 1, 2],  ## for IDEA StaveBased
    ChargeCollectionMethod="LookupTable",  # "Debug" / "SinglePixel" / "LookupTable"
    LookupTableFile="SET_FROM_COMMAND_LINE",
    LookupTableShiftTruthPosition=True,
    ClusterPositionUncertainty=[],
    Threshold=100,
    ThresholdDispersion=5,
    ChargeSmearing=10,
    TimeSmearing=10,
    OutputLevel=VERBOSE,
)

from Configurables import AuditorSvc, ChronoAuditor, UniqueIDGenSvc

chronoAuditor = ChronoAuditor()
auditorSvc = AuditorSvc(Auditors=[chronoAuditor])

from Configurables import HepRndm__Engine_CLHEP__RanluxEngine_ as RndmEngine

rndmEngine = RndmEngine("RndmGenSvc.Engine", SetSingleton=True, Seeds=[1234567])

from Configurables import RndmGenSvc

rndmGenSvc = RndmGenSvc("RndmGenSvc", Engine=rndmEngine.name())

from Configurables import EventDataSvc
from k4FWCore import ApplicationMgr

applicationMgr = ApplicationMgr(
    TopAlg=[
        vtxb_digitizer,
    ],
    EvtSel="NONE",
    ExtSvc=[
        EventDataSvc("EventDataSvc"),
        geoSvc,
        auditorSvc,
        UniqueIDGenSvc("uidSvc"),
        rndmEngine,
        rndmGenSvc,
    ],
    StopOnSignal=True,
)

for algo in applicationMgr.TopAlg:
    algo.AuditExecute = True
