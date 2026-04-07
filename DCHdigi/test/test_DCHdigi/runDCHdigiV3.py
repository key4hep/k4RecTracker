#
# gaudi steering file that runs DCHdigi
#
# to execute:
# k4run runDCHdigi.py

from Gaudi.Configuration import INFO,DEBUG
from Configurables import EventDataSvc, UniqueIDGenSvc
from k4FWCore import ApplicationMgr, IOSvc
from Gaudi.Configuration import *

svc = IOSvc("IOSvc")
svc.Input = [ "Muon_data/dch_muon_10GeV_N100.root"]
svc.Output = "Muon_data/dch_muon_10GeV_digi_N100.root"

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = ['./compact/DCH_standalone_o1_v02.xml']

from Configurables import DCHdigi_v03

DCHdigi = DCHdigi_v03("DCHdigi")
DCHdigi.DCH_simhits=["DCHCollection"]
DCHdigi.DCH_name="DCH_v2"
DCHdigi.create_debug_histograms=True

DCHdigi.zResolution_mm=1
DCHdigi.xyResolution_mm=0.1

#Cluster Calculation
DCHdigi.MeanExcitationEnergy_eV = 48.48
DCHdigi.GasDensity_g_cm3 = 3.984e-4
DCHdigi.W_eff_eV = 35.0

#LUT
DCHdigi.XTFileName = "par.root"

#Avalanche
DCHdigi.GasGain = 2.0e5 
DCHdigi.PolyaTheta = 0.5

#Analog Waveform
DCHdigi.SignalPulseWindow_ns = 200.0
DCHdigi.PulseBinSize_ns = 1.0
DCHdigi.PulseAmplitudeScale = 5.e-6
DCHdigi.PulseRiseTime_ns = 1.0
DCHdigi.PulseFallTime_ns = 7.0
DCHdigi.WaveformTimeWindow_ns = 3000.0
DCHdigi.WaveformBinSize_ns    = 0.5
DCHdigi.EnableWireAttenuation = True

#Electronics
DCHdigi.FrontEndGain = 1.0
DCHdigi.EnableImpedanceMismatch = True
DCHdigi.MatchingRes_Ohm = 330.0
DCHdigi.TubeImpedance_Ohm = 50.0
DCHdigi.EnableElectronicsTF = True
DCHdigi.ElectronicsTauRise_ns = 3.0
DCHdigi.ElectronicsTauFall1_ns = 9.0
DCHdigi.ElectronicsTauFall2_ns = 25.0
DCHdigi.ElectronicsMixFraction = 0.5
DCHdigi.ElectronicsKernelLength_ns = 200.0

#Noise
DCHdigi.EnableFFTNoise = True
DCHdigi.NoiseFileName = "par.root"
DCHdigi.NoiseDirName = "fft_noise"
DCHdigi.NoiseScale = 1.e-3
DCHdigi.NoiseRemoveDC = True

#Digitization
DCHdigi.DoWaveformDigitization = True
DCHdigi.ADCSamplePeriod_ns   = 0.5
DCHdigi.ADCPolarity = 1.0
DCHdigi.ADCBits = 12
DCHdigi.ADCLSB_mV = 1.0
DCHdigi.ADCBaseline_mV = 50.0

DCHdigi.OutputLevel=INFO

THistSvc().Output = ["rec DATAFILE='DCHdigi_debug.root' TYP='ROOT' OPT='RECREATE'"]
THistSvc().PrintAll = True
THistSvc().AutoSave = True
THistSvc().AutoFlush = True
THistSvc().OutputLevel = INFO

mgr = ApplicationMgr(
    TopAlg=[DCHdigi],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice,EventDataSvc("EventDataSvc"),UniqueIDGenSvc("uidSvc")],
    OutputLevel=INFO,
)
