import os

from Gaudi.Configuration import *

from Configurables import FCCDataSvc
podioevent  = FCCDataSvc("EventDataSvc")

from GaudiKernel.SystemOfUnits import MeV, GeV, tesla
################## Particle gun setup
momentum = 5 # in GeV
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
path_to_detector = os.environ.get("FCCDETECTORS", "")
print(path_to_detector)
detectors_to_use=[
                    'Detector/DetFCCeeIDEA/compact/FCCee_DectMaster.xml',
                  ]
# prefix all xmls with path_to_detector
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
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
geantservice = SimG4Svc("SimG4Svc", detector = 'SimG4DD4hepDetector', physicslist = "SimG4FtfpBert", actions = actions, magneticField = field)

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
saveDCHsimHitTool = SimG4SaveTrackerHits("saveDCHsimHitTool", readoutNames=["SimplifiedDriftChamberCollection"])
saveDCHsimHitTool.SimTrackHits.Path = "DC_simTrackerHits"

saveVTXBsimHitTool = SimG4SaveTrackerHits("saveVTXBsimHitTool", readoutNames=["VertexBarrelCollection"])
saveVTXBsimHitTool.SimTrackHits.Path = "VTXB_simTrackerHits"

saveVTXEsimHitTool = SimG4SaveTrackerHits("saveVTXEsimHitTool", readoutNames=["VertexEndcapCollection"])
saveVTXEsimHitTool.SimTrackHits.Path = "VTXE_simTrackerHits"

from Configurables import SimG4Alg
geantsim = SimG4Alg("SimG4Alg",
                       outputs= [saveVTXBsimHitTool, saveVTXEsimHitTool, saveDCHsimHitTool,
                                 #saveHistTool
                       ],
                       eventProvider = particle_converter,
                       OutputLevel = INFO)

# Digitize tracker hits (for now digitization is reconstruction)
vtxb_reco_hit_name = saveVTXBsimHitTool.SimTrackHits.Path.replace("sim", "reco")
from Configurables import VTXdigitizer
vtxb_digitizer = VTXdigitizer("VTXBdigitizer",
    inputSimHits = saveVTXBsimHitTool.SimTrackHits.Path,
    outputDigiHits = vtxb_reco_hit_name
)

vtxe_reco_hit_name = saveVTXEsimHitTool.SimTrackHits.Path.replace("sim", "reco")
vtxe_digitizer = VTXdigitizer("VTXEdigitizer",
    inputSimHits = saveVTXEsimHitTool.SimTrackHits.Path,
    outputDigiHits = vtxe_reco_hit_name
)

dch_reco_hit_name = saveDCHsimHitTool.SimTrackHits.Path.replace("sim", "reco")
from Configurables import DCHdigitizer
dch_digitizer = DCHdigitizer("DCHdigitizer",
    inputSimHits = saveDCHsimHitTool.SimTrackHits.Path,
    outputDigiHits = dch_reco_hit_name
)

# run the genfit tracking 
from Configurables import GenFitter
genfitter = GenFitter("GenFitter", inputHits = dch_reco_hit_name, outputTracks = "genfit_tracks") 

################ Output
from Configurables import PodioOutput
out = PodioOutput("out",
                  OutputLevel=INFO)
out.outputCommands = ["keep *"]

import uuid
out.filename = "output_simplifiedDriftChamber_MagneticField_"+str(magneticField)+"_pMin_"+str(momentum*1000)+"_MeV"+"_ThetaMinMax_"+str(thetaMin)+"_"+str(thetaMax)+"_pdgId_"+str(pdgCode)+".root"

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
              vtxb_digitizer,
              vtxe_digitizer,
              dch_digitizer,
              genfitter,
              out
              ],
    EvtSel = 'NONE',
    EvtMax   = 4,
    ExtSvc = [geoservice, podioevent, geantservice, audsvc],
    StopOnSignal = True,
 )
