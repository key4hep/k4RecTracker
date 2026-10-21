import os
import math
from Gaudi.Configuration import INFO, DEBUG
from Configurables import EventDataSvc, UniqueIDGenSvc
from Configurables import RndmGenSvc
from Configurables import SimG4SaveTrackerHits
from Configurables import GeoSvc

from k4FWCore import IOSvc, ApplicationMgr
from k4FWCore.parseArgs import parser


################## Parser
parser.add_argument("--inputFile", help="InputFile")
parser.add_argument("--outputFile", help="OutputFile")
parser.add_argument("--modelPath", default="", help="model path for the track finder")
parser.add_argument("--tbeta", default=0.7, help="tbeta clustering parameter")
parser.add_argument("--td", default=0.1, help="td clustering parameter")
parser.add_argument("--inferenceThreads", default=4, type=int, help="ONNX Runtime intra-operation threads")
parser.add_argument("--sessionRecreateEvents", default=10, type=int, help="Recreate the ONNX session every N events (0 disables)")
args = parser.parse_args()

# ################## InputOutput
svc = IOSvc("IOSvc")
svc.Input = args.inputFile
svc.Output = args.outputFile

################ Detector geometry
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
detectors_to_use = ["FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml"]
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

# different sensors for inner/outer barrel layers
# see https://indico.cern.ch/event/1244371/contributions/5350233
innerVertexResolution_x = 0.003  # [mm], assume 3 µm resolution for ARCADIA sensor
innerVertexResolution_y = 0.003  # [mm], assume 3 µm resolution for ARCADIA sensor
innerVertexResolution_t = 1000  # [ns]
outerVertexResolution_x = 0.050 / math.sqrt(12)  # [mm], assume ATLASPix3 sensor with 50 µm pitch
outerVertexResolution_y = 0.150 / math.sqrt(12)  # [mm], assume ATLASPix3 sensor with 150 µm pitch
outerVertexResolution_t = 1000  # [ns]

# silicon wrapper hits parameters
siWrapperResolution_x = 0.050 / math.sqrt(12)  # [mm]
siWrapperResolution_y = 1.0 / math.sqrt(12)  # [mm]
siWrapperResolution_t = 0.040  # [ns], assume 40 ps timing resolution for a single layer -> Should lead to <30 ps resolution when >1 hit

# Define arguments for digitizers
vxd_barrel_digi_args = {
    "IsStrip": False,
    "ResolutionU": [innerVertexResolution_x]*3 + [outerVertexResolution_x]*2,
    "ResolutionV": [innerVertexResolution_y]*3 + [outerVertexResolution_y]*2,
    "ResolutionT": [innerVertexResolution_t]*3 + [outerVertexResolution_t]*2,
    "SimTrackHitCollectionName": ["VertexBarrelCollection"],
    "SimTrkHitRelCollection": ["VTXBSimDigiLinks"],
    # IDEA_o1_v04 exposes both barrel and disks below the single top-level
    # DD4hep detector element named "Vertex". The collection name identifies
    # which hits this digitizer processes.
    "SubDetectorName": "Vertex",
    "TrackerHitCollectionName": ["VTXBDigis"],
    "ForceHitsOntoSurface": True,
    "CellIDBits": 32,
}

vxd_endcap_digi_args = {
    "IsStrip": False,
    "ResolutionU": [outerVertexResolution_x]*3,
    "ResolutionV": [outerVertexResolution_y]*3,
    "ResolutionT": [outerVertexResolution_t]*3,
    "SimTrackHitCollectionName": ["VertexEndcapCollection"],
    "SimTrkHitRelCollection": ["VTXDSimDigiLinks"],
    "SubDetectorName": "Vertex",
    "TrackerHitCollectionName": ["VTXDDigis"],
    "ForceHitsOntoSurface": True,
    "CellIDBits": 32,
}

siWr_barrel_digi_args = {
    "IsStrip": False,
    "ResolutionU": [siWrapperResolution_x]*4,
    "ResolutionV": [siWrapperResolution_y]*4,
    "ResolutionT": [siWrapperResolution_t]*4,
    "SimTrackHitCollectionName": ["SiWrBCollection"],
    "SimTrkHitRelCollection": ["SiWrBSimDigiLinks"],
    "SubDetectorName": "SiWrB",
    "TrackerHitCollectionName": ["SiWrBDigis"],
    "ForceHitsOntoSurface": True,
    "CellIDBits": 32,
}

siWr_endcap_digi_args = {
    "IsStrip": False,
    "ResolutionU": [siWrapperResolution_x]*4,
    "ResolutionV": [siWrapperResolution_y]*4,
    "ResolutionT": [siWrapperResolution_t]*4,
    "SimTrackHitCollectionName": ["SiWrDCollection"],
    "SimTrkHitRelCollection": ["SiWrDSimDigiLinks"],
    "SubDetectorName": "SiWrD",
    "TrackerHitCollectionName": ["SiWrDDigis"],
    "ForceHitsOntoSurface": True,
    "CellIDBits": 32,
}

# digitize vertex hits through "native" DDPlanarDigi
from Configurables import DDPlanarDigi

vtxb_digitizer = DDPlanarDigi(
    "VXDBarrelDigitizer",
    **vxd_barrel_digi_args,
    OutputLevel=INFO
)

vtxd_digitizer = DDPlanarDigi(
    "VXDEndcapDigitizer",
    **vxd_endcap_digi_args,
    OutputLevel=INFO
)

siwrb_digitizer = DDPlanarDigi(
    "SiWrBarrelDigitizer",
    **siWr_barrel_digi_args,
    OutputLevel=INFO
)

siwrd_digitizer = DDPlanarDigi(
    "SiWrEndcapDigitizer",
    **siWr_endcap_digi_args,
    OutputLevel=INFO
)

from Configurables import DCHdigi_v02
dch_digitizer = DCHdigi_v02(
    "DCHdigi2",
    InputSimHitCollection=["DCHCollection"],
    OutputDigihitCollection = ["DCH_DigiCollection"],
    OutputLinkCollection = ["DCH_DigiSimAssociationCollection"],
    DCH_name="DCH_v2",
    zResolution_mm = 30.,               # in mm
    xyResolution_mm = 0.1,              # in mm
    Deadtime_ns = 400.0,                # in ns
    GasType=0,                          # 0: He(90%)-Isobutane(10%), 1: pure He, 2: Ar(50%)-Ethane(50%), 3: pure Ar
    ReadoutWindowStartTime_ns=1.0,      # in ns (taking into account time of flight, drift, and signal travel)
    ReadoutWindowDuration_ns=450.0,     # in ns
    DriftVelocity_um_per_ns=-1.0,       # in um/ns, if negative, automatically chosen based on GasType
    SignalVelocity_mm_per_ns=200.0,     # in mm/ns (Default: 2/3 of the speed of light)
    OutputLevel=INFO,
)

############### Track Finder
from Configurables import GGTFTrackFinder

GGTF = GGTFTrackFinder(
    "GGTFTrackFinder",
    InputPlanarHitCollections=["VTXBDigis", "VTXDDigis", "SiWrBDigis", "SiWrDDigis"],
    InputWireHitCollections=["DCH_DigiCollection"],
    OutputTracksGGTF=["GGTFTracks"],
    ModelPath=args.modelPath,
    EmbeddingCoordinates=4,
    Tbeta=args.tbeta,
    Td=args.td,
    InferenceThreads=args.inferenceThreads,
    SessionRecreateEvents=args.sessionRecreateEvents,
    OutputLevel=INFO,
)

############### Application Manager
mgr = ApplicationMgr(
    TopAlg=[dch_digitizer, vtxb_digitizer, vtxd_digitizer, siwrb_digitizer, siwrd_digitizer, GGTF],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice, EventDataSvc("EventDataSvc"), UniqueIDGenSvc("uidSvc"), RndmGenSvc()],
    OutputLevel=INFO,
)
