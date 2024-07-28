from k4FWCore import ApplicationMgr
from Gaudi.Configuration import INFO, WARNING, DEBUG
import os

if not os.path.isfile("ddsim_output_edm4hep.root"):
    os.system("ddsim --enableGun --gun.distribution uniform --gun.energy '10*GeV' --gun.particle e- --numberOfEvents 100 --outputFile ddsim_output_edm4hep.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/IDEA/compact/IDEA_o1_v02/IDEA_o1_v02.xml")

# Loading the output of the SIM step
from k4FWCore import IOSvc
io_svc = IOSvc("IOSvc")
io_svc.input = "ddsim_output_edm4hep.root"
io_svc.output = "tracks_from_genParticle_output.root"

# Calling TracksFromGenParticles
from Configurables import TracksFromGenParticles
tracksFromGenParticles = TracksFromGenParticles("TracksFromGenParticles",
                                               InputGenParticles = ["MCParticles"],
                                               OutputTracks = ["TracksFromGenParticles"],
                                               OutputMCRecoTrackParticleAssociation = ["TracksFromGenParticlesAssociation"],
                                               Bz = 2.0,
                                               OutputLevel = INFO)

# produce a TH1 with distances between tracks and simTrackerHits
from Configurables import PlotTrackHitDistances
plotTrackHitDistances = PlotTrackHitDistances("PlotTrackHitDistances",
                                             InputSimTrackerHits = ["CDCHHits"],
                                             InputTracksFromGenParticlesAssociation = tracksFromGenParticles.OutputMCRecoTrackParticleAssociation, 
                                             Bz = 2.0)

# Set auditor service
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
tracksFromGenParticles.AuditExecute = True
plotTrackHitDistances.AuditExecute = True

# event counter
from Configurables import EventCounter
event_counter = EventCounter('event_counter')
event_counter.Frequency = 10

from Configurables import EventDataSvc
ApplicationMgr(
    TopAlg= [event_counter, tracksFromGenParticles, plotTrackHitDistances],
    EvtSel='NONE',
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), audsvc],
    StopOnSignal=True,
)
