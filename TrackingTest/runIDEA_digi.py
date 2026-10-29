# Steering: IDEA_o1_v04 simulation file (ddsim) -> digitized hits.
# Produces the shared input for the tracking comparison:
#   planar digis: VTXBDigis, VTXDDigis, SiWrBDigis, SiWrDDigis (+ SimDigiLinks)
#   wire digis:   DCH_DigiCollection (+ association)
# Digitization settings as in the IDEA v4 o1 production steering:
# /afs/cern.ch/work/a/adevita/public/trainingFolder/GGTF_trackFinder_training/data_creation/condor_pipeline/IDEA/noBackground_parquet/utils/runIDEA_v4o1_trackerDigitizer.py
#
# Usage: k4run runIDEA_digi.py --inputFile sim.edm4hep.root --outputFile digi.edm4hep.root

import os
import math

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc, UniqueIDGenSvc, RndmGenSvc, GeoSvc
from Configurables import DDPlanarDigi
from Configurables import DCHdigi_v02
from k4FWCore import IOSvc, ApplicationMgr
from k4FWCore.parseArgs import parser

################## Parser
parser.add_argument(
    "--inputFile", default="sim.edm4hep.root", help="InputFile (ddsim output, IDEA_o1_v04)"
)
parser.add_argument("--outputFile", default="digi.edm4hep.root", help="OutputFile")
args = parser.parse_args()

################## InputOutput
svc = IOSvc("IOSvc")
svc.Input = args.inputFile
svc.Output = args.outputFile

################ Detector geometry
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
detectors_to_use = ["FCCee/IDEA/compact/IDEA_o1_v04/IDEA_o1_v04.xml"]
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

# different sensors for inner/outer barrel layers
# see https://indico.cern.ch/event/1244371/contributions/5350233
innerVertexResolution_x = 0.003  # [mm], assume 3 um resolution for ARCADIA sensor
innerVertexResolution_y = 0.003  # [mm], assume 3 um resolution for ARCADIA sensor
innerVertexResolution_t = 25  # [ns]
outerVertexResolution_x = 0.050 / math.sqrt(12)  # [mm], assume ATLASPix3 sensor with 50 um pitch
outerVertexResolution_y = 0.150 / math.sqrt(12)  # [mm], assume ATLASPix3 sensor with 150 um pitch
outerVertexResolution_t = 25  # [ns]

# silicon wrapper hits parameters
siWrapperResolution_x = 0.050 / math.sqrt(12)  # [mm]
siWrapperResolution_y = 1.0 / math.sqrt(12)  # [mm]
siWrapperResolution_t = (
    0.040  # [ns], assume 40 ps timing resolution for a single layer -> <30 ps when >1 hit
)

############### Vertex Digitizer

vtxb_digitizer = DDPlanarDigi("VTXBdigitizer")
vtxb_digitizer.SubDetectorName = "Vertex"
vtxb_digitizer.CellIDBits = 32
vtxb_digitizer.IsStrip = False
vtxb_digitizer.ResolutionU = [
    innerVertexResolution_x,
    innerVertexResolution_x,
    innerVertexResolution_x,
    outerVertexResolution_x,
    outerVertexResolution_x,
]
vtxb_digitizer.ResolutionV = [
    innerVertexResolution_y,
    innerVertexResolution_y,
    innerVertexResolution_y,
    outerVertexResolution_y,
    outerVertexResolution_y,
]
vtxb_digitizer.ResolutionT = [
    innerVertexResolution_t,
    innerVertexResolution_t,
    innerVertexResolution_t,
    outerVertexResolution_t,
    outerVertexResolution_t,
]
vtxb_digitizer.SimTrackHitCollectionName = ["VertexBarrelCollection"]
vtxb_digitizer.SimTrkHitRelCollection = ["VTXBSimDigiLinks"]
vtxb_digitizer.TrackerHitCollectionName = ["VTXBDigis"]
vtxb_digitizer.ForceHitsOntoSurface = True

vtxd_digitizer = DDPlanarDigi("VTXDdigitizer")
vtxd_digitizer.SubDetectorName = "Vertex"
vtxd_digitizer.CellIDBits = 32
vtxd_digitizer.IsStrip = False
vtxd_digitizer.ResolutionU = [
    outerVertexResolution_x,
    outerVertexResolution_x,
    outerVertexResolution_x,
]
vtxd_digitizer.ResolutionV = [
    outerVertexResolution_y,
    outerVertexResolution_y,
    outerVertexResolution_y,
]
vtxd_digitizer.ResolutionT = [
    outerVertexResolution_t,
    outerVertexResolution_t,
    outerVertexResolution_t,
]
vtxd_digitizer.SimTrackHitCollectionName = ["VertexEndcapCollection"]
vtxd_digitizer.SimTrkHitRelCollection = ["VTXDSimDigiLinks"]
vtxd_digitizer.TrackerHitCollectionName = ["VTXDDigis"]
vtxd_digitizer.ForceHitsOntoSurface = True

############### Wrapper Digitizer

siwrb_digitizer = DDPlanarDigi("SiWrBdigitizer")
siwrb_digitizer.SubDetectorName = "SiWrB"
siwrb_digitizer.CellIDBits = 32
siwrb_digitizer.IsStrip = False
siwrb_digitizer.ResolutionU = [siWrapperResolution_x] * 4
siwrb_digitizer.ResolutionV = [siWrapperResolution_y] * 4
siwrb_digitizer.ResolutionT = [siWrapperResolution_t] * 4
siwrb_digitizer.SimTrackHitCollectionName = ["SiWrBCollection"]
siwrb_digitizer.SimTrkHitRelCollection = ["SiWrBSimDigiLinks"]
siwrb_digitizer.TrackerHitCollectionName = ["SiWrBDigis"]
siwrb_digitizer.ForceHitsOntoSurface = True

siwrd_digitizer = DDPlanarDigi("SiWrDdigitizer")
siwrd_digitizer.SubDetectorName = "SiWrD"
siwrd_digitizer.CellIDBits = 32
siwrd_digitizer.IsStrip = False
siwrd_digitizer.ResolutionU = [siWrapperResolution_x] * 4
siwrd_digitizer.ResolutionV = [siWrapperResolution_y] * 4
siwrd_digitizer.ResolutionT = [siWrapperResolution_t] * 4
siwrd_digitizer.SimTrackHitCollectionName = ["SiWrDCollection"]
siwrd_digitizer.SimTrkHitRelCollection = ["SiWrDSimDigiLinks"]
siwrd_digitizer.TrackerHitCollectionName = ["SiWrDDigis"]
siwrd_digitizer.ForceHitsOntoSurface = True

############### DCH Digitizer

dch_digitizer = DCHdigi_v02(
    "DCHdigi2",
    InputSimHitCollection=["DCHCollection"],
    OutputDigihitCollection=["DCH_DigiCollection"],
    OutputLinkCollection=["DCH_DigiSimAssociationCollection"],
    DCH_name="DCH_v2",
    zResolution_mm=30.0,  # in mm
    xyResolution_mm=0.1,  # in mm
    Deadtime_ns=400.0,  # in ns
    GasType=0,  # 0: He(90%)-Isobutane(10%), 1: pure He, 2: Ar(50%)-Ethane(50%), 3: pure Ar
    ReadoutWindowStartTime_ns=1.0,  # in ns (taking into account time of flight, drift, and signal travel)
    ReadoutWindowDuration_ns=450.0,  # in ns
    DriftVelocity_um_per_ns=-1.0,  # in um/ns, if negative, automatically chosen based on GasType
    SignalVelocity_mm_per_ns=200.0,  # in mm/ns (Default: 2/3 of the speed of light)
    OutputLevel=INFO,
)

############### Application Manager

mgr = ApplicationMgr(
    TopAlg=[dch_digitizer, vtxd_digitizer, vtxb_digitizer, siwrb_digitizer, siwrd_digitizer],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice, EventDataSvc("EventDataSvc"), UniqueIDGenSvc("uidSvc"), RndmGenSvc()],
    OutputLevel=INFO,
)
