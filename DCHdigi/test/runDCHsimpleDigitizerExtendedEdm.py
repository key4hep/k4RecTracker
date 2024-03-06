import os

from Gaudi.Configuration import *

from Configurables import FCCDataSvc
podioevent  = FCCDataSvc("EventDataSvc")

from GaudiKernel.SystemOfUnits import MeV, GeV, tesla
################## Particle gun setup
momentum = 10 # in GeV
#thetaMin = 90.25 # degrees
#thetaMax = 90.25 # degrees
thetaMin = 20 # degrees
thetaMax = 130 # degrees
pdgCode = 13 # 11 electron, 13 muon, 22 photon, 111 pi0, 211 pi+
magneticField = True
_pi = 3.14159

from Configurables import GenAlg
genAlg = GenAlg()
from Configurables import  MomentumRangeParticleGun
pgun = MomentumRangeParticleGun("ParticleGun_Electron")
pgun.PdgCodes = [pdgCode]
pgun.MomentumMin = momentum * GeV
pgun.MomentumMax = momentum * GeV
pgun.PhiMin = 0
pgun.PhiMax = 2 * _pi
pgun.ThetaMin = thetaMin * _pi / 180.
pgun.ThetaMax = thetaMax * _pi / 180.
genAlg.SignalProvider = pgun

genAlg.hepmc.Path = "hepmc"

from Configurables import HepMCToEDMConverter
hepmc_converter = HepMCToEDMConverter()
hepmc_converter.hepmc.Path="hepmc"
genParticlesOutputName = "genParticles"
hepmc_converter.GenParticles.Path = genParticlesOutputName
hepmc_converter.hepmcStatusList = []

################## Simulation setup
# Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
# if FCC_DETECTORS is empty, this should use relative path to working directory
path_to_detector = os.environ.get("K4GEO", "")
print(path_to_detector)
detectors_to_use=[
                    'FCCee/IDEA/compact/IDEA_o1_v02/IDEA_o1_v02.xml'
                  ]
# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

# define regions to act on Geant4 max step length to have the expected 1hit/800um in the drift chamber
#from Configurables import SimG4FullSimDCHRegion
#drift_chamber_region = SimG4FullSimDCHRegion("drift_chamber_region", volumeNames = ["CDCH"], max_step_length = 10)
from Configurables import  SimG4UserLimitRegion
regiontool = SimG4UserLimitRegion("limits")
regiontool.volumeNames = ["CDCH"]
regiontool.OutputLevel = DEBUG
from GaudiKernel.SystemOfUnits import mm
#regiontool.maxStep = 0.4*mm

from Configurables import SimG4UserLimitPhysicsList
# create overlay on top of FTFP_BERT physics list, attaching fast sim/parametrization process
physicslisttool = SimG4UserLimitPhysicsList("Physics")
physicslisttool.fullphysics = "SimG4FtfpBert"

# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions

from Configurables import SimG4FullSimActions, SimG4Alg, SimG4PrimariesFromEdmTool, SimG4SaveParticleHistory
actions = SimG4FullSimActions()
# Uncomment if history from Geant4 decays is needed (e.g. to get the photons from pi0) and set actions=actions in SimG4Svc + uncomment saveHistTool in SimG4Alg
#actions.enableHistory=True
#actions.energyCut = 0.2 * GeV
#saveHistTool = SimG4SaveParticleHistory("saveHistory")


# Magnetic field
from Configurables import SimG4ConstantMagneticFieldTool
field = SimG4ConstantMagneticFieldTool("SimG4ConstantMagneticFieldTool", FieldComponentZ = -2 * tesla, FieldOn = magneticField, IntegratorStepper = "ClassicalRK4")

from Configurables import SimG4Svc
#geantservice = SimG4Svc("SimG4Svc", detector = 'SimG4DD4hepDetector', physicslist = "SimG4FtfpBert", actions = actions, magneticField = field)
geantservice = SimG4Svc("SimG4Svc", detector = 'SimG4DD4hepDetector', physicslist = physicslisttool, actions = actions, magneticField = field, regions = [regiontool])
#geantservice = SimG4Svc("SimG4Svc", detector = 'SimG4DD4hepDetector', physicslist = "SimG4FtfpBert", actions = actions, magneticField = field, regions = ["SimG4FullSimDCHRegion/drift_chamber_region"])

# Fixed seed to have reproducible results, change it for each job if you split one production into several jobs
# Mind that if you leave Gaudi handle random seed and some job start within the same second (very likely) you will have duplicates
geantservice.randomNumbersFromGaudi = False
geantservice.seedValue = 4242

# Range cut
geantservice.g4PreInitCommands += ["/run/setCut 0.1 mm"]
#geantservice.g4PreInitCommands += ["/run/setStepMax 0.1 mm"]

# Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools
# and a tool that saves the calorimeter hits

# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
from Configurables import SimG4PrimariesFromEdmTool
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.GenParticles.Path = genParticlesOutputName

from Configurables import SimG4SaveTrackerHits
savetrackertool = SimG4SaveTrackerHits("SimG4SaveTrackerHits", readoutName="CDCHHits")
savetrackertool.SimTrackHits.Path = "DC_simTrackerHits"


from Configurables import SimG4Alg
geantsim = SimG4Alg("SimG4Alg",
                       outputs= [savetrackertool
                                 #saveHistTool
                       ],
                       eventProvider=particle_converter,
                       OutputLevel=DEBUG)
# Digitize tracker hits
from Configurables import DCHsimpleDigitizerExtendedEdm
dch_digitizer = DCHsimpleDigitizerExtendedEdm("DCHsimpleDigitizerExtendedEdm",
    inputSimHits = savetrackertool.SimTrackHits.Path,
    outputDigiHits = savetrackertool.SimTrackHits.Path.replace("sim", "digi"),
    outputSimDigiAssociation = "DC_simDigiAssociation",
    readoutName = "CDCHHits",
    xyResolution = 0.1, # mm
    zResolution = 1, # mm
    debugMode = False,
    OutputLevel = INFO
)

# Derive performance quantities
#from Configurables import DCHsimpleDigitizerExtendedEdmPerformance
#dch_perf = DCHsimpleDigitizerExtendedEdmPerformance("DCHsimpleDigitizerExtendedEdmPerformance",
#        inputSimHits = savetrackertool.SimTrackHits.Path,
#        inputDigiHits = dch_digitizer.outputDigiHits,
#        inputSimDigiAssociation = dch_digitizer.outputSimDigiAssociation
#)



################ Output
from Configurables import PodioOutput
out = PodioOutput("out",
                  OutputLevel=INFO)
out.outputCommands = ["keep *"]
if not dch_digitizer.debugMode:
    out.outputCommands.append("drop *HitSimHitDelta*")
    out.outputCommands.append("drop outputDigiLocalHits")

import uuid
out.filename = "output_simplifiedDriftChamber_MagneticField_"+str(magneticField)+"_pMin_"+str(momentum*1000)+"_MeV"+"_ThetaMinMax_"+str(thetaMin)+"_"+str(thetaMax)+"_pdgId_"+str(pdgCode)+"_stepLength_default.root"

#CPU information
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
genAlg.AuditExecute = True
hepmc_converter.AuditExecute = True
geantsim.AuditExecute = True
out.AuditExecute = True

from Configurables import EventCounter
event_counter = EventCounter('event_counter')
event_counter.Frequency = 1

from Configurables import ApplicationMgr
ApplicationMgr(
    TopAlg = [
              event_counter,
              genAlg,
              hepmc_converter,
              geantsim,
              dch_digitizer,
              #dch_perf,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = 100,
    ExtSvc = [geoservice, podioevent, geantservice, audsvc],
    StopOnSignal = True,
 )
