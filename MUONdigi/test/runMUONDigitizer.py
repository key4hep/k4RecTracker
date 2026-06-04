import os, math

from Gaudi.Configuration import *

from Configurables import EventDataSvc
from k4FWCore import ApplicationMgr, IOSvc

podioevent = EventDataSvc("EventDataSvc")

from GaudiKernel.SystemOfUnits import MeV, GeV, tesla

################## Particle gun setup
momentum = 5            # in GeV
thetaMin = 0            # degrees
thetaMax = 180          # degrees
pdgCode  = 13           # 13 = mu+/-
magneticField = True
_pi = 3.14159

################## Muon-system digi parameters (cluster-based, uRWELL)
detectorName         = "Muon-System"
uResolution          = 0.4    # mm   surface u (= chamber-local y)
vResolution          = 0.4    # mm   surface v (= chamber-local z)
tResolution          = -1.0   # ns   <=0 disables time smearing (digi time = earliest SimHit time)
efficiency           = 0.98   # per-cluster detection efficiency
timeWindow           = 25.0   # ns   SimHits within this window go in one cluster
maxClusterSize       = 3      # max distinct strips per view (Y, Z) inside one cluster
forceHitsOntoSurface = True   # if True, project smeared digi back onto the chamber surface if it lands outside
printFrequency       = 1      # INFO line every N events (set 0 to silence, 100 for production)

from Configurables import GenAlg

genAlg = GenAlg()
from Configurables import MomentumRangeParticleGun

pgun = MomentumRangeParticleGun("ParticleGun_Muon")
pgun.PdgCodes    = [pdgCode]
pgun.MomentumMin = momentum * GeV
pgun.MomentumMax = momentum * GeV
pgun.PhiMin      = 0
pgun.PhiMax      = 2 * _pi
pgun.ThetaMin    = thetaMin * _pi / 180.0
pgun.ThetaMax    = thetaMax * _pi / 180.0
genAlg.SignalProvider = pgun
genAlg.hepmc.Path = "hepmc"

from Configurables import HepMCToEDMConverter

hepmc_converter = HepMCToEDMConverter()
hepmc_converter.hepmc.Path = "hepmc"
genParticlesOutputName = "genParticles"
hepmc_converter.GenParticles.Path = genParticlesOutputName
hepmc_converter.hepmcStatusList = []

################## Simulation setup
# Detector geometry
from Configurables import GeoSvc

geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
print("K4GEO =", path_to_detector)
detectors_to_use = [
    "FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml",
]
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

# Geant4
from Configurables import (
    SimG4FullSimActions,
    SimG4Alg,
    SimG4PrimariesFromEdmTool,
    SimG4SaveParticleHistory,
)

actions = SimG4FullSimActions()

# Magnetic field
from Configurables import SimG4ConstantMagneticFieldTool

field = SimG4ConstantMagneticFieldTool(
    "SimG4ConstantMagneticFieldTool",
    FieldComponentZ=-2 * tesla,
    FieldOn=magneticField,
    IntegratorStepper="ClassicalRK4",
)

from Configurables import SimG4Svc

geantservice = SimG4Svc(
    "SimG4Svc",
    detector="SimG4DD4hepDetector",
    physicslist="SimG4FtfpBert",
    magneticField=field,
)
geantservice.randomNumbersFromGaudi = False
geantservice.seedValue = 4242
geantservice.g4PreInitCommands += ["/run/setCut 0.1 mm"]

# Geant4 algorithm
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.GenParticles.Path = genParticlesOutputName

from Configurables import SimG4SaveTrackerHits

SimG4SaveMuonHits = SimG4SaveTrackerHits(
    "SimG4SaveMuonHits", readoutName="MuonSystemCollection"
)
SimG4SaveMuonHits.SimTrackHits.Path = "MuonSystemCollection"

geantsim = SimG4Alg(
    "SimG4Alg",
    outputs=[SimG4SaveMuonHits],
    eventProvider=particle_converter,
    OutputLevel=INFO,
)

# Digitize muon-system hits — cluster-based
from Configurables import MUONDigitizer

muon_digitizer = MUONDigitizer(
    "MUONDigitizer",
    inputSimHits         = SimG4SaveMuonHits.SimTrackHits.Path,
    outputDigiHits       = "MSTrackerHits",
    outputSimDigiLink    = "MSTrackerHitRelations",
    detectorName         = detectorName,
    readoutName          = "MuonSystemCollection",
    uResolution          = uResolution,
    vResolution          = vResolution,
    tResolution          = tResolution,
    efficiency           = efficiency,
    timeWindow           = timeWindow,
    maxClusterSize       = maxClusterSize,
    forceHitsOntoSurface = forceHitsOntoSurface,
    printFrequency       = printFrequency,
    OutputLevel          = INFO,
)

################ Output
iosvc = IOSvc()
iosvc.Output = (
    "output_MuonSystemDigi"
    + "_B" + str(magneticField)
    + "_pMin_" + str(momentum * 1000) + "_MeV"
    + "_ThetaMinMax_" + str(thetaMin) + "_" + str(thetaMax)
    + "_pdgId_" + str(pdgCode)
    + ".root"
)
iosvc.outputCommands = ["keep *"]

# CPU information
from Configurables import AuditorSvc, ChronoAuditor

chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
genAlg.AuditExecute          = True
hepmc_converter.AuditExecute = True
geantsim.AuditExecute        = True
muon_digitizer.AuditExecute  = True

ApplicationMgr(
    TopAlg=[
        genAlg,
        hepmc_converter,
        geantsim,
        muon_digitizer,
    ],
    EvtSel="NONE",
    EvtMax=10,
    ExtSvc=[geoservice, podioevent, geantservice, audsvc],
    StopOnSignal=True,
)
