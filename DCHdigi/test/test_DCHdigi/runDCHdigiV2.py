#
# gaudi steering file that runs DCHdigi
#
# to execute:
# k4run runDCHdigi.py

from Gaudi.Configuration import INFO,DEBUG
from Configurables import EventDataSvc, UniqueIDGenSvc
from k4FWCore import ApplicationMgr, IOSvc

svc = IOSvc("IOSvc")
svc.Input = [ "dch_proton_10GeV.root"]
svc.Output = "dch_proton_10GeV_digi.root"

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = ['./compact/DCH_standalone_o1_v02.xml']

from Configurables import DCHdigi_v02
DCHdigi = DCHdigi_v02("DCHdigi2",
                      InputSimHitCollection = ["DCHCollection"],
                      DCH_name = "DCH_v2",
                      zResolution_mm = 1., # in mm 
                      xyResolution_mm = 0.1, # in mm
                      GasSel = 0, 
                      OutputLevel = INFO)

mgr = ApplicationMgr(
    TopAlg=[DCHdigi],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice,EventDataSvc("EventDataSvc"),UniqueIDGenSvc("uidSvc")],
    OutputLevel=INFO,
)
