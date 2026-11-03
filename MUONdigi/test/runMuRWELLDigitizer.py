import os, subprocess

from Gaudi.Configuration import *

from Configurables import EventDataSvc
from k4FWCore import ApplicationMgr, IOSvc

podioevent = EventDataSvc("EventDataSvc")

################## Particle gun setup
momentum = 5  # in GeV
thetaMin = 10  # degrees
thetaMax = 170  # degrees
particle = "mu-"
pdgCode = 13  # 13 = mu+/- (only used in the output file name)
magneticField = True  # ddsim uses the field defined in the detector compact
numberOfEvents = 10

################## Muon-system digi parameters (cluster-based, uRWELL)
detectorName = "Muon-System"
uResolution = 0.4  # mm   surface u (= chamber-local y)
vResolution = 0.4  # mm   surface v (= chamber-local z)
tResolution = -1.0  # ns   <=0 disables time smearing (digi time = earliest SimHit time)
efficiency = 0.98  # per-cluster detection efficiency
timeWindow = 25.0  # ns   SimHits within this window go in one cluster
maxClusterSize = 3  # max distinct strips per view (Y, Z) inside one cluster
forceHitsOntoSurface = (
    True  # if True, project smeared digi back onto the chamber surface if it lands outside
)
printFrequency = 1  # INFO line every N events (set 0 to silence, 100 for production)

################## Detector geometry
from Configurables import GeoSvc

geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
print("K4GEO =", path_to_detector)
detectors_to_use = [
    "FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml",
]
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False  # Geant4 runs in ddsim below, not inside this job

################## Simulation setup
ideaSteeringFile = os.path.join(
    os.environ.get("FCCCONFIG", ""), "FullSim/IDEA/IDEA_o1_v03/SteeringFile_IDEA_o1_v03.py"
)
simFile = "MuonSystemSim_pdgId_" + str(pdgCode) + ".root"
subprocess.run(
    [
        "ddsim",
        "--steeringFile",
        ideaSteeringFile,
        "--compactFile",
        geoservice.detectors[0],
        "--enableGun",
        "--gun.particle",
        particle,
        "--gun.energy",
        str(momentum) + "*GeV",
        "--gun.distribution",
        "uniform",
        "--gun.thetaMin",
        str(thetaMin) + "*deg",
        "--gun.thetaMax",
        str(thetaMax) + "*deg",
        "--numberOfEvents",
        str(numberOfEvents),
        "--random.seed",
        "4242",
        "--outputFile",
        simFile,
    ],
    check=True,  # stop the test if the simulation fails
)

# Random numbers for the digitizer smearing and efficiency
from Configurables import RndmGenSvc

randomService = RndmGenSvc("RndmGenSvc")

# Digitize muon-system hits — cluster-based
from Configurables import MuRWELLDigitizer

muon_digitizer = MuRWELLDigitizer(
    "MuRWELLDigitizer",
    inputSimHits="MuonSystemCollection",
    outputDigiHits="MSTrackerHits",
    outputSimDigiLink="MSTrackerHitRelations",
    detectorName=detectorName,
    readoutName="MuonSystemCollection",
    uResolution=uResolution,
    vResolution=vResolution,
    tResolution=tResolution,
    efficiency=efficiency,
    timeWindow=timeWindow,
    maxClusterSize=maxClusterSize,
    forceHitsOntoSurface=forceHitsOntoSurface,
    printFrequency=printFrequency,
    OutputLevel=INFO,
)

################ Input / Output
iosvc = IOSvc()
iosvc.Input = simFile
iosvc.Output = (
    "output_MuonSystemDigi"
    + "_B"
    + str(magneticField)
    + "_pMin_"
    + str(momentum * 1000)
    + "_MeV"
    + "_ThetaMinMax_"
    + str(thetaMin)
    + "_"
    + str(thetaMax)
    + "_pdgId_"
    + str(pdgCode)
    + ".root"
)
iosvc.outputCommands = ["keep *"]

# CPU information
from Configurables import AuditorSvc, ChronoAuditor

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
muon_digitizer.AuditExecute = True

ApplicationMgr(
    TopAlg=[
        muon_digitizer,
    ],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice, podioevent, randomService, audsvc],
    StopOnSignal=True,
)
