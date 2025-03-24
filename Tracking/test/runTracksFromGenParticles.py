from k4FWCore import ApplicationMgr
from Gaudi.Configuration import INFO, DEBUG
import os

# Loading the output of the SIM step
from k4FWCore import IOSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = "ddsim_output_edm4hep.root"
io_svc.Output = "tracks_from_genParticle_output.root"

# Geometry service
from Configurables import GeoSvc
import os
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
detectors_to_use = [
    'FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml'
]
geoservice.detectors = [
    os.path.join(path_to_detector, _det) for _det in detectors_to_use
]
geoservice.OutputLevel = INFO

# Calling TracksFromGenParticles
from Configurables import TracksFromGenParticles
tracksFromGenParticles = TracksFromGenParticles("TracksFromGenParticles",
                                                InputGenParticles=["MCParticles"],
                                                InputSimTrackerHits=["VertexBarrelCollection",
                                                                     "VertexEndcapCollection",
                                                                     "DCHCollection",
                                                                     "SiWrBCollection",
                                                                     "SiWrDCollection"],
                                                OutputTracks=["TracksFromGenParticles"],
                                                OutputMCRecoTrackParticleAssociation=["TracksFromGenParticlesAssociation"],
                                                TrackerIDs=[1,2,3,23,24],
                                                ExtrapolateToECal=True,
                                                OutputLevel=DEBUG)

# produce a TH1 with distances between tracks and simTrackerHits
from Configurables import PlotTrackHitDistances, RootHistSvc
from Configurables import Gaudi__Histograming__Sink__Root as RootHistoSink
plotTrackHitDistances = PlotTrackHitDistances("PlotTrackHitDistances",
                                              InputSimTrackerHits=["DCHCollection"],
                                              InputTracksFromGenParticlesAssociation=["TracksFromGenParticlesAssociation"],
                                              Bz=2.0)

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
    TopAlg=[tracksFromGenParticles, plotTrackHitDistances],
    EvtSel='NONE',
    EvtMax=-1,
    ExtSvc=[root_hist_svc, EventDataSvc("EventDataSvc"), audsvc, geoservice],
    StopOnSignal=True,
)
