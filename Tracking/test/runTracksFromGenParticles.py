from Configurables import ApplicationMgr
from Configurables import EventCounter
from Gaudi.Configuration import INFO, WARNING, DEBUG
import os

if not os.path.isfile("ddsim_output_edm4hep.root"):
    os.system("ddsim --enableGun --gun.distribution uniform --gun.energy '10*GeV' --gun.particle e- --numberOfEvents 100 --outputFile ddsim_output_edm4hep.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/IDEA/compact/IDEA_o1_v02/IDEA_o1_v02.xml")

# Loading the output of the SIM step
from Configurables import k4DataSvc, PodioInput, PodioOutput
evtsvc = k4DataSvc('EventDataSvc')
evtsvc.input = "ddsim_output_edm4hep.root"
#evtsvc.input = "/afs/cern.ch/user/b/brfranco/work/public/MIT_tutorial/CLDConfig/CLDConfig/wzp6_ee_mumuH_ecm240_CLD_RECO_edm4hep.root"
Nevts = -1
input_reader = PodioInput('InputReader')

# Calling TracksFromGenParticles
from Configurables import TracksFromGenParticles
tracksFromGenParticles = TracksFromGenParticles("TracksFromGenParticles",
                                               InputGenParticles = "MCParticles",
                                               OutputTracks = "TracksFromGenParticles",
                                               OutputMCRecoTrackParticleAssociation = "TracksFromGenParticlesAssociation",
                                               Bz = 2.0,
                                               OutputLevel = INFO)

# produce a TH1 with distances between tracks and simTrackerHits
from Configurables import PlotTrackHitDistances
plotTrackHitDistances = PlotTrackHitDistances("PlotTrackHitDistances",
                                             InputSimTrackerHits = "CDCHHits",
                                             InputTracksFromGenParticlesAssociation = tracksFromGenParticles.OutputMCRecoTrackParticleAssociation, 
                                             Bz = 2.0)

# Output
from Configurables import PodioOutput
out = PodioOutput("out",
                  OutputLevel=INFO)
out.outputCommands = ["keep *"]
out.filename = "tracks_from_genParticle_output.root"

# Set auditor service
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]

ApplicationMgr(
    TopAlg= [input_reader, tracksFromGenParticles, plotTrackHitDistances, out],
    EvtSel='NONE',
    EvtMax=Nevts,
    ExtSvc=[evtsvc, audsvc],
    StopOnSignal=True,
)
