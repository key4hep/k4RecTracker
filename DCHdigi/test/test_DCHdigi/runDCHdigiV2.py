#
# gaudi steering file that runs DCHdigi
#
# to execute:
# k4run runDCHdigiV2.py

from Gaudi.Configuration import INFO, DEBUG
from Configurables import EventDataSvc, UniqueIDGenSvc
from k4FWCore import ApplicationMgr, IOSvc

svc = IOSvc("IOSvc")
svc.Input = ["dch_proton_10GeV.root"]
svc.Output = "dch_proton_10GeV_digi.root"

from Configurables import GeoSvc

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = ["./compact/DCH_standalone_o1_v02.xml"]

from Configurables import DCHdigi_v02

DCHdigi = DCHdigi_v02(
    "DCHdigi2",
    InputSimHitCollection=["DCHCollection"],
    DCH_name="DCH_v2",
    zResolution_mm=1.0,  # in mm
    xyResolution_mm=0.1,  # in mm
    Deadtime_ns=400.0,  # in ns
    GasType=0,  # 0: He(90%)-Isobutane(10%), 1: pure He, 2: Ar(50%)-Ethane(50%), 3: pure Ar
    ReadoutWindowStartTime_ns=1.0,  # in ns (taking into account time of flight, drift, and signal travel)
    ReadoutWindowDuration_ns=450.0,  # in ns
    DriftVelocity_um_per_ns=-1.0,  # in um/ns, if negative, automatically chosen based on GasType
    SignalVelocity_mm_per_ns=200.0,  # in mm/ns (Default: 2/3 of the speed of light)
    OutputLevel=INFO,
)

mgr = ApplicationMgr(
    TopAlg=[DCHdigi],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice, EventDataSvc("EventDataSvc"), UniqueIDGenSvc("uidSvc")],
    OutputLevel=INFO,
)
