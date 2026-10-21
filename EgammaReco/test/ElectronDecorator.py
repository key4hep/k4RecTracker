from Configurables import ElectronDecorator
from k4FWCore import ApplicationMgr, IOSvc
from Gaudi.Configuration import DEBUG, INFO
from Configurables import GeoSvc
from Configurables import UniqueIDGenSvc
import os

io = IOSvc(
    Input=["output_pandora_electron_1000.root"],
    Output="output_pandora.root",
    CollectionNames=["PandoraPFOs","PandoraClusters","SiTracks_Refitted"]
)

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
    os.environ["K4GEO"] + "/FCCee/CLD/compact/CLD_o2_v08/CLD_o2_v08.xml"
]
geoservice.OutputLevel = INFO

alg = ElectronDecorator(
    "ElectronDecorator",
    InputRecoParticles="PandoraPFOs",
    OutputRecoParticles="RecoveredPFOs",
    dEtaMargin = 0.005,
    dPhiMargin = 0.05
)

ApplicationMgr(
    TopAlg=[alg],
    ExtSvc=[io],
    EvtSel="NONE",
    EvtMax=1000,
    OutputLevel=DEBUG,
)
