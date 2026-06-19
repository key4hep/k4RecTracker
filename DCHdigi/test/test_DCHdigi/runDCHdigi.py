#
# gaudi steering file that runs DCHdigi
#
# to execute:
# k4run runDCHdigi.py

import os

from Gaudi.Configuration import INFO, DEBUG
from Configurables import EventDataSvc, UniqueIDGenSvc
from k4FWCore import ApplicationMgr, IOSvc

CI_TYPE = os.environ.get('CI_TYPE', '')

svc = IOSvc("IOSvc")
svc.Input = [f"dch_proton_10GeV_{CI_TYPE}.root"]
svc.Output = f"dch_proton_10GeV_digi_{CI_TYPE}.root"

from Configurables import GeoSvc

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = ["./compact/DCH_standalone_o1_v02.xml"]

from Configurables import DCHdigi_v01

DCHdigi = DCHdigi_v01("DCHdigi")
DCHdigi.DCH_simhits = ["DCHCollection"]
DCHdigi.DCH_name = "DCH_v2"
DCHdigi.fileDataAlg = "DataAlgFORGEANT.root"
DCHdigi.calculate_dndx = True
DCHdigi.create_debug_histograms = True
DCHdigi.out_debug_filename = f"dch_digi_alg_debug_{CI_TYPE}.root"
DCHdigi.zResolution_mm = 1
DCHdigi.xyResolution_mm = 0.1


DCHdigi.OutputLevel = INFO

mgr = ApplicationMgr(
    TopAlg=[DCHdigi],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice, EventDataSvc("EventDataSvc"), UniqueIDGenSvc("uidSvc")],
    OutputLevel=INFO,
)
