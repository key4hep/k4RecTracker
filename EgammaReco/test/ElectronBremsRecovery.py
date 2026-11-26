from Configurables import ElectronBremsRecovery
from k4FWCore import ApplicationMgr, IOSvc
from Gaudi.Configuration import DEBUG, INFO
from Configurables import GeoSvc
import os

io = IOSvc(
    Input=["./CLDConfig/CLDConfig/CLD_SIM.edm4hep.Output_REC.edm4hep.root"],
    Output="output_pandora_bremsRecovered.root",
    CollectionNames=["PandoraPFOs", "PandoraClusters", "SiTracks_Refitted"],
)

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [
    os.environ["K4GEO"] + "/FCCee/CLD/compact/CLD_o2_v08/CLD_o2_v08.xml"
]  # same as the simulation and reconstruction in test_ElectronBremsRecovery.sh
geoservice.OutputLevel = INFO

alg = ElectronBremsRecovery(
    "ElectronBremsRecovery",
    InputRecoParticles="PandoraPFOs",
    OutputRecoParticles="RecoveredPFOs",
    dEtaMargin=0.005,
    dPhiMargin=0.05,
)

ApplicationMgr(
    TopAlg=[alg],
    ExtSvc=[io],
    EvtSel="NONE",
    EvtMax=10,
    OutputLevel=DEBUG,
)
