from k4FWCore import ApplicationMgr
from Gaudi.Configuration import INFO, WARNING, DEBUG
import os

if not os.path.isfile("ddsim_output_edm4hep.root"):
    os.system("ddsim --enableGun --gun.distribution uniform --gun.energy '10*GeV' --gun.particle e- --numberOfEvents 100 --outputFile ddsim_output_edm4hep.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml")

# Loading the output of the SIM step
from k4FWCore import IOSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = "ddsim_output_edm4hep.root"
io_svc.Output = "tracks_from_genParticle_output.root"

# Calling TracksFromGenParticles
from Configurables import TracksFromGenParticles
tracksFromGenParticles = TracksFromGenParticles("TracksFromGenParticles",
                                               InputGenParticles = ["MCParticles"],
                                               OutputTracks = ["TracksFromGenParticles"],
                                               OutputMCRecoTrackParticleAssociation = ["TracksFromGenParticlesAssociation"],
                                               Bz = 2.0,
                                               OutputLevel = INFO)

# produce a TH1 with distances between tracks and simTrackerHits
from Configurables import PlotTrackHitDistances, RootHistSvc
from Configurables import Gaudi__Histograming__Sink__Root as RootHistoSink
plotTrackHitDistances = PlotTrackHitDistances("PlotTrackHitDistances",
                                             InputSimTrackerHits = ["DCHCollection"],
                                             InputTracksFromGenParticlesAssociation = tracksFromGenParticles.OutputMCRecoTrackParticleAssociation, 
                                             Bz = 2.0)

hps = RootHistSvc("HistogramPersistencySvc")
root_hist_svc = RootHistoSink("RootHistoSink")
root_hist_svc.FileName = "TrackHitDistances.root"

# Set auditor service
from Configurables import AuditorSvc, ChronoAuditor
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]
tracksFromGenParticles.AuditExecute = True
plotTrackHitDistances.AuditExecute = True

from Configurables import EventDataSvc
ApplicationMgr(
    TopAlg= [tracksFromGenParticles, plotTrackHitDistances],
    EvtSel='NONE',
    EvtMax=-1,
    ExtSvc=[root_hist_svc, EventDataSvc("EventDataSvc"), audsvc],
    StopOnSignal=True,
)
