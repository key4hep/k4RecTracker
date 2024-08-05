#
# gaudi steering file that runs DCHdigi
#
# to execute:
# k4run runDCHdigi.py

from Gaudi.Configuration import INFO,DEBUG
from Configurables import EventDataSvc, UniqueIDGenSvc
from k4FWCore import ApplicationMgr, IOSvc

svc = IOSvc("IOSvc")
svc.input = [ "dch_proton_10GeV.root"]
svc.output = "dch_proton_10GeV_digi.root"

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = ['/home/alvarotd/work/mydch/DCH_detector/compact/DCH_standalone_o1_v02.xml']

from Configurables import DCHdigi
DCHdigi = DCHdigi("DCHdigi")
DCHdigi.DCH_simhits=["DCHCollection"]
DCHdigi.DCH_name="DCH_v2"
DCHdigi.fileDataAlg="DataAlgFORGEANT.root"
DCHdigi.create_debug_histograms=True
DCHdigi.zResolution_mm=1.23
DCHdigi.xyResolution_mm=4.56


DCHdigi.OutputLevel=INFO

mgr = ApplicationMgr(
    TopAlg=[DCHdigi],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice,EventDataSvc("EventDataSvc"),UniqueIDGenSvc("uidSvc")],
    OutputLevel=INFO,
)
