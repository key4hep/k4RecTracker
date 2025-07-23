import os
import math 
from Gaudi.Configuration import INFO,DEBUG
from Configurables import EventDataSvc,UniqueIDGenSvc
from Configurables import RndmGenSvc
from Configurables import SimG4SaveTrackerHits
from Configurables import GeoSvc

from k4FWCore import IOSvc,ApplicationMgr
from k4FWCore.parseArgs import parser


################## Parser
parser.add_argument("--inputFile", default="ddsim_output_edm4hep.root", help="InputFile")
parser.add_argument("--outputFile", default="output_digi.root", help="OutputFile")
parser.add_argument("--modelPath", default="", help="model path for the track finder")
parser.add_argument("--tbeta", default=0.6, help="tbeta clustering parameter")
parser.add_argument("--td", default=0.3, help="td clustering parameter")
args = parser.parse_args()

# ################## InputOutput
svc = IOSvc("IOSvc")
svc.Input = args.inputFile
svc.Output = args.outputFile

################ Detector geometry
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
detectors_to_use=['FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml']
geoservice.detectors = [os.path.join(path_to_detector, _det) for _det in detectors_to_use]
geoservice.OutputLevel = INFO

from Configurables import DDPlanarDigi
from Configurables import DCHdigi_v01

############### Vertex Digitizer

innerVertexResolution_x = 0.003 # [mm], assume 3 µm resolution for ARCADIA sensor
innerVertexResolution_y = 0.003 # [mm], assume 3 µm resolution for ARCADIA sensor
innerVertexResolution_t = 1000 # [ns]
outerVertexResolution_x = 0.050/math.sqrt(12) # [mm], assume ATLASPix3 sensor with 50 µm pitch
outerVertexResolution_y = 0.150/math.sqrt(12) # [mm], assume ATLASPix3 sensor with 150 µm pitch
outerVertexResolution_t = 1000 # [ns]

vtxb_digitizer = DDPlanarDigi("VTXBdigitizer")
vtxb_digitizer.SubDetectorName = "Vertex"
vtxb_digitizer.IsStrip = False
vtxb_digitizer.ResolutionU = [innerVertexResolution_x, innerVertexResolution_x, innerVertexResolution_x, outerVertexResolution_x, outerVertexResolution_x]
vtxb_digitizer.ResolutionV = [innerVertexResolution_y, innerVertexResolution_y, innerVertexResolution_y, outerVertexResolution_y, outerVertexResolution_y]
vtxb_digitizer.ResolutionT = [innerVertexResolution_t, innerVertexResolution_t, innerVertexResolution_t, outerVertexResolution_t, outerVertexResolution_t]
vtxb_digitizer.SimTrackHitCollectionName = ["VertexBarrelCollection"]
vtxb_digitizer.SimTrkHitRelCollection = ["VTXBSimDigiLinks"]
vtxb_digitizer.TrackerHitCollectionName = ["VTXBDigis"]
vtxb_digitizer.ForceHitsOntoSurface = True

vtxd_digitizer = DDPlanarDigi("VTXDdigitizer")
vtxd_digitizer.SubDetectorName = "Vertex"
vtxd_digitizer.IsStrip = False
vtxd_digitizer.ResolutionU = [outerVertexResolution_x, outerVertexResolution_x, outerVertexResolution_x]
vtxd_digitizer.ResolutionV = [outerVertexResolution_y, outerVertexResolution_y, outerVertexResolution_y]
vtxd_digitizer.ResolutionT = [outerVertexResolution_t, outerVertexResolution_t, outerVertexResolution_t]
vtxd_digitizer.SimTrackHitCollectionName = ["VertexEndcapCollection"]
vtxd_digitizer.SimTrkHitRelCollection = ["VTXDSimDigiLinks"]
vtxd_digitizer.TrackerHitCollectionName = ["VTXDDigis"]
vtxd_digitizer.ForceHitsOntoSurface = True

############### Wrapper Digitizer
siWrapperResolution_x   = 0.050/math.sqrt(12) # [mm]
siWrapperResolution_y   = 1.0/math.sqrt(12) # [mm]
siWrapperResolution_t   = 0.040 # [ns], assume 40 ps timing resolution for a single layer -> Should lead to <30 ps resolution when >1 hit

siwrb_digitizer = DDPlanarDigi("SiWrBdigitizer")
siwrb_digitizer.SubDetectorName = "SiWrB"
siwrb_digitizer.IsStrip = False
siwrb_digitizer.ResolutionU = [siWrapperResolution_x, siWrapperResolution_x]
siwrb_digitizer.ResolutionV = [siWrapperResolution_y, siWrapperResolution_y]
siwrb_digitizer.ResolutionT = [siWrapperResolution_t, siWrapperResolution_t]
siwrb_digitizer.SimTrackHitCollectionName = ["SiWrBCollection"]
siwrb_digitizer.SimTrkHitRelCollection = ["SiWrBSimDigiLinks"]
siwrb_digitizer.TrackerHitCollectionName = ["SiWrBDigis"]
siwrb_digitizer.ForceHitsOntoSurface = True

siwrd_digitizer = DDPlanarDigi("SiWrDdigitizer")
siwrd_digitizer.SubDetectorName = "SiWrD"
siwrd_digitizer.IsStrip = False
siwrd_digitizer.ResolutionU = [siWrapperResolution_x, siWrapperResolution_x]
siwrd_digitizer.ResolutionV = [siWrapperResolution_y, siWrapperResolution_y]
siwrd_digitizer.ResolutionT = [siWrapperResolution_t, siWrapperResolution_t]
siwrd_digitizer.SimTrackHitCollectionName = ["SiWrDCollection"]
siwrd_digitizer.SimTrkHitRelCollection = ["SiWrDSimDigiLinks"]
siwrd_digitizer.TrackerHitCollectionName = ["SiWrDDigis"]
siwrd_digitizer.ForceHitsOntoSurface = True


############### DCH Digitizer

dch_digitizer = DCHdigi_v01("DCHdigi",
    DCH_simhits = ["DCHCollection"],
    DCH_name = "DCH_v2",
    fileDataAlg = "DataAlgFORGEANT.root",
    calculate_dndx = False, # cluster counting disabled (to be validated, see FCC-config#239)
    create_debug_histograms = False,
    zResolution_mm = 30., # in mm - Note: At this point, the z resolution comes without the stereo measurement
    xyResolution_mm = 0.1 # in mm
)

############### Track Finder
from Configurables import GGTF_tracking

GGTF = GGTF_tracking(
    "GGTF_tracking",
    inputWireHits=["DCH_DigiCollection"],
    inputPlanarHits=["VTXBDigis", "VTXDDigis", "SiWrBDigis", "SiWrDDigis"],
    outputTracks=["CDCHTracks"],
    OutputLevel=INFO,
)
GGTF.modelPath = args.modelPath
GGTF.tbeta = args.tbeta
GGTF.td = args.td

############### Application Manager
import subprocess

ifilename = "https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/IDEA/DataAlgFORGEANT.root"
subprocess.run(["wget", "--no-clobber", ifilename])

mgr = ApplicationMgr(TopAlg=[dch_digitizer, vtxb_digitizer, vtxd_digitizer, siwrb_digitizer, siwrd_digitizer,GGTF],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice,EventDataSvc("EventDataSvc"),UniqueIDGenSvc("uidSvc"),RndmGenSvc()],
    OutputLevel=INFO,
    )