# CombinatorialTrackFinder on the IDEA muon system: tracks are found from the muon-chamber
# hits alone (no inner tracker), so muons from decays far from the IP are found too.
# Settings: Exhaustive, with the values of the HNL muon-system production
# (IDEA_muon_optimization/steering_tracking/track_gun_<layout>_100k.py).
#
# Stored track states (all EDM4hep convention):
#   AtFirstHit  the fitted helix at the innermost muon hit
#   AtLastHit   the same helix at the outermost muon hit
#   AtIP        the track at its closest approach to the origin: the outer helix is followed
#               back to the solenoid and matched there to a helix in the inner field
#               (IPStateMode = SolenoidCrossing)
# Note: the muon chambers sit in the return field (B < 0), so the sign of omega is MINUS
# the charge of the muon (omega > 0 means clockwise seen from +z).
#
# usage:
#   k4run runMuonSystemTracking.py --inputFile digi.root --outputFile tracks.root \
#         --compact <compact the input was digitised with> [--numEvents N] [--genfit]
#   --genfit: the found tracks are then fitted with k4RecTracker's GenfitTrackFitter
#             (output collection MuonSystemGenfitTracks)
#   compact examples: <k4geo>/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml
#                     <k4geo>/FCCee/IDEA/compact/IDEA_o1_v05/IDEA_o1_v05.xml
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser
from Configurables import EventDataSvc, GeoSvc, Gaudi__Monitoring__MessageSvcSink
from Configurables import CombinatorialTrackFinder

parser.add_argument("--inputFile", required=True, help="digitised file with MSTrackerHits")
parser.add_argument("--outputFile", required=True)
parser.add_argument(
    "--compact",
    required=True,
    help="geometry compact; MUST be the one the input was digitised with",
)
parser.add_argument("--numEvents", type=int, default=-1, help="-1 = all events")
parser.add_argument(
    "--genfit", action="store_true", help="also fit the found tracks with GenfitTrackFitter"
)
args = parser.parse_args()

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [args.compact]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

t = CombinatorialTrackFinder("MuonSystemTracking")
t.OutputLevel = INFO

# -- detector and readout --
t.DetectorNames = ["Muon-System"]
t.EncodingStringParameterName = "MuonSystemReadoutID"
t.RegionFieldName = "type"  # 0 barrel, +1/-1 endcaps

# -- stored track states --
t.IPStateMode = "SolenoidCrossing"  # AtIP from the helix matched at the solenoid boundary

# -- Exhaustive: all k-hit combinations over the layers --
t.Strategy = "Exhaustive"
t.MinTrackHits = 3
t.MaxCombinatorialHits = 6
t.NHitMaxChi2NDF = 10.0  # circle-fit chi2/ndf cut
t.MaxComboPT = 250.0  # GeV; higher-pT combinations only as fallback
t.SigmaHitDefault = 0.04  # cm, used when a hit has no resolution
t.OutlierSigma = 2.5  # outlier removal at 2.5 x SigmaHitDefault
t.HitIsolationCut = 0  # off
t.DoHitClustering = True  # independent hit pools per road

# road cuts: hits of one track stay in a narrow azimuth window
t.MaxConsecDeltaPhi = 0.15  # rad, between consecutive layers
t.MaxComboPhiSpread = 0.3  # rad, over the whole track
t.MaxConsecutiveHitDistance = 200.0  # cm
t.MaxPairHitDistance = 450.0  # cm, allows barrel -> endcap tracks

# energy-deposit cuts against delta rays and hadronic secondaries
t.MaxHitEdepKeV = 10  # hits above 10 keV are not used
t.MaxHitsPerCompositeID = 10  # at most 10 hits per layer, best energy-deposit score first

# stricter cuts for a track next to an accepted one (delta-ray shadows)
t.NeighbourTrackMaxDist = 10.0  # cm
t.NeighbourTrackMinHits = 4

# -- collections --
t.InputHitCollection = ["MSTrackerHits"]
t.InputRecoSimLinkCollection = ["MSTrackerHitRelations"]
t.OutputTrackCollection = ["MuonSystemTracks"]
t.OutputHitCollection = ["MuonSystemTrackHits"]
t.OutputRecoSimLinkCollection = ["MuonSystemTrackHitSimLinks"]

algorithms = [t]

# -- optional: Kalman fit of the found tracks with GenfitTrackFitter --
if args.genfit:
    from Configurables import GenfitTrackFitter

    fitter = GenfitTrackFitter("MuonSystemGenfitTrackFitter")
    fitter.InputTracks = ["MuonSystemTracks"]
    fitter.OutputFittedTracks = ["MuonSystemGenfitTracks"]
    fitter.OutputFittedTracksWithFilteredHits = ["MuonSystemGenfitTracksWithFilteredHits"]
    fitter.OutputFittedHits = ["MuonSystemGenfitHits"]
    fitter.ParticleHypothesisList = [13]  # muons
    fitter.RunSingleEvaluation = True
    fitter.InitializationType = 1
    fitter.ListOfTypesToSkip = []  # fit every track type
    fitter.RunCalorimeterExtrapolation = False  # the tracks are already outside the calorimeter
    fitter.OutputLevel = INFO
    algorithms.append(fitter)

iosvc = IOSvc()
iosvc.Input = args.inputFile
iosvc.Output = args.outputFile

ApplicationMgr(
    TopAlg=algorithms,
    EvtSel="NONE",
    EvtMax=args.numEvents,
    # MessageSvcSink prints the run statistics (counters) of the algorithms at the end of the job
    ExtSvc=[EventDataSvc("EventDataSvc"), Gaudi__Monitoring__MessageSvcSink()],
    OutputLevel=INFO,
)
