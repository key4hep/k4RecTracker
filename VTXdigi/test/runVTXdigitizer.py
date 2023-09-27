import os

from Gaudi.Configuration import *

from Configurables import FCCDataSvc
podioevent  = FCCDataSvc("EventDataSvc")

from GaudiKernel.SystemOfUnits import MeV, GeV, tesla
################## Particle gun setup
momentum = 5 # in GeV
thetaMin = 5
thetaMax = 175
pdgCode = 13 # 11 electron, 13 muon, 22 photon, 111 pi0, 211 pi+
magneticField = True
_pi = 3.14159

################## Vertex sensor resolutions
innerVertexResolution_x = 0.005 # [mm], assume 5 µm resolution for ARCADIA sensor
innerVertexResolution_y = 0.005 # [mm], assume 5 µm resolution for ARCADIA sensor
outerVertexResolution_x = 0.014 # [mm], assume ATLASPix3 sensor with 50 µm pitch -> 50/sqrt(12) = 14.4 µm resolution
outerVertexResolution_y = 0.043 # [mm], assume ATLASPix3 sensor with 150 µm pitch -> 150/sqrt(12) = 43.3 µm resolution

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
# path_to_detector = os.environ.get("FCCDETECTORS", "")
# path_to_detector = os.environ.get("FCCDETECTORS", "")
# print(path_to_detector)
# detectors_to_use=[
#                     'FCCee/IDEA/compact/IDEA_o1_v01/FCCee_IDEA_o1_v01.xml',
#                   ]
# prefix all xmls with path_to_detector
# geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.detectors = ["../lcgeo/FCCee/IDEA/compact/IDEA_o1_v01/IDEA_o1_v01.xml"] # IDEA
geoservice.OutputLevel = INFO

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
geantservice = SimG4Svc("SimG4Svc", detector = 'SimG4DD4hepDetector', physicslist = "SimG4FtfpBert", magneticField = field) # , actions = actions)

# Fixed seed to have reproducible results, change it for each job if you split one production into several jobs
# Mind that if you leave Gaudi handle random seed and some job start within the same second (very likely) you will have duplicates
geantservice.randomNumbersFromGaudi = False
geantservice.seedValue = 4242

# Range cut
geantservice.g4PreInitCommands += ["/run/setCut 0.1 mm"]

# Geant4 algorithm
# Translates EDM to G4Event, passes the event to G4, writes out outputs via tools
# and a tool that saves the calorimeter hits

# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
from Configurables import SimG4PrimariesFromEdmTool
particle_converter = SimG4PrimariesFromEdmTool("EdmConverter")
particle_converter.GenParticles.Path = genParticlesOutputName

from Configurables import SimG4SaveTrackerHits

### CLD
# savetrackertool = SimG4SaveTrackerHits("SimG4SaveTrackerHits", readoutName="VertexBarrelCollection")
# savetrackertool.SimTrackHits.Path = "VTXB_simTrackerHits"

# savetrackertool = SimG4SaveTrackerHits("SimG4SaveTrackerHits", readoutName="VertexEndcapCollection")
# savetrackertool.SimTrackHits.Path = "VTXE_simTrackerHits"


### IDEA
SimG4SaveTrackerHitsIB = SimG4SaveTrackerHits("SimG4SaveTrackerHitsIB", readoutName="VTXIBCollection")
SimG4SaveTrackerHitsIB.SimTrackHits.Path = "VTXIB_simTrackerHits"

SimG4SaveTrackerHitsOB = SimG4SaveTrackerHits("SimG4SaveTrackerHitsOB", readoutName="VTXOBCollection")
SimG4SaveTrackerHitsOB.SimTrackHits.Path = "VTXOB_simTrackerHits"

SimG4SaveTrackerHitsD = SimG4SaveTrackerHits("SimG4SaveTrackerHitsD", readoutName="VTXDCollection")
SimG4SaveTrackerHitsD.SimTrackHits.Path = "VTXD_simTrackerHits"


from Configurables import SimG4Alg
geantsim = SimG4Alg("SimG4Alg",
                       outputs= [SimG4SaveTrackerHitsIB, SimG4SaveTrackerHitsOB, SimG4SaveTrackerHitsD
                                 #saveHistTool
                       ],
                       eventProvider=particle_converter,
                       OutputLevel=INFO)
# Digitize tracker hits
from Configurables import VTXdigitizer

### For CLD
# vtx_digitizer = VTXdigitizer("VTXdigitizer",
#     inputSimHits = savetrackertool.SimTrackHits.Path,
#     outputDigiHits = savetrackertool.SimTrackHits.Path.replace("sim", "digi"),
#     readoutName = "VertexBarrelCollection",
#     xResolution = 0.005, # mm
#     yResolution = 0.005, # mm
#     tResolution = 1000, # ns
#     OutputLevel = INFO
# )

# vtx_digitizer = VTXdigitizer("VTXdigitizer",
#     inputSimHits = savetrackertool.SimTrackHits.Path,
#     outputDigiHits = savetrackertool.SimTrackHits.Path.replace("sim", "digi"),
#     readoutName = "VertexEndcapCollection",
#     xResolution = 0.005, # mm
#     yResolution = 0.005, # mm
#     tResolution = 1000, # ns
#     OutputLevel = INFO
# )


### For IDEA
vtxib_digitizer = VTXdigitizer("VTXIBdigitizer",
    inputSimHits = SimG4SaveTrackerHitsIB.SimTrackHits.Path,
    outputDigiHits = SimG4SaveTrackerHitsIB.SimTrackHits.Path.replace("sim", "digi"),
    readoutName = "VTXIBCollection",
    xResolution = innerVertexResolution_x, # mm, r-phi direction
    yResolution = innerVertexResolution_y, # mm, z direction
    tResolution = 1000, # ns
    OutputLevel = INFO
)

vtxob_digitizer = VTXdigitizer("VTXOBdigitizer",
    inputSimHits = SimG4SaveTrackerHitsOB.SimTrackHits.Path,
    outputDigiHits = SimG4SaveTrackerHitsOB.SimTrackHits.Path.replace("sim", "digi"),
    readoutName = "VTXOBCollection",
    xResolution = outerVertexResolution_x, # mm, r-phi direction
    yResolution = outerVertexResolution_y, # mm, z direction
    tResolution = 1000, # ns
    OutputLevel = INFO
)

vtxd_digitizer = VTXdigitizer("VTXDdigitizer",
    inputSimHits = SimG4SaveTrackerHitsD.SimTrackHits.Path,
    outputDigiHits = SimG4SaveTrackerHitsD.SimTrackHits.Path.replace("sim", "digi"),
    readoutName = "VTXDCollection",
    xResolution = outerVertexResolution_x, # mm, r direction
    yResolution = outerVertexResolution_y, # mm, phi direction
    tResolution = 1000, # ns
    OutputLevel = INFO
)

# run the genfit tracking 
# from Configurables import GenFitter
# genfitter = GenFitter("GenFitter", inputHits = savetrackertool.SimTrackHits.Path.replace("sim", "digi"), outputTracks = "genfit_tracks") 

################ Output
from Configurables import PodioOutput
out = PodioOutput("out",
                  OutputLevel=INFO)
out.outputCommands = ["keep *"]

import uuid
out.filename = "output_vertex_"+str(magneticField)+"_pMin_"+str(momentum*1000)+"_MeV"+"_ThetaMinMax_"+str(thetaMin)+"_"+str(thetaMax)+"_pdgId_"+str(pdgCode)+".root"

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
              #vtx_digitizer,
              vtxib_digitizer,
              vtxob_digitizer,
              vtxd_digitizer,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = 1000,
    ExtSvc = [geoservice, podioevent, geantservice, audsvc],
    StopOnSignal = True
 )
