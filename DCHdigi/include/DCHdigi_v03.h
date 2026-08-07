/** ======= DCHdigi_v03 ==========
 * Gaudi Algorithm for DCH digitization
 *
 * @author Muhammad Saiel, Giovanni Francesco Tassielli, Nicola De Filippis
 * @date   2026-02
 *
 * <h4>Input collections and prerequisites</h4>
 * Processor requires a collection of SimTrackerHits <br>
 * This code uses DD4hep length natural unit (cm), but EDM4hep data is (usually) in mm. Please be careful with units. <br>
 *
 * <br> <h4>Output</h4>
 * Processor produces digitized hits and optional waveform collections. <br>
 *
 * @param DCH_simhits The name of input collection, type edm4hep::SimTrackerHitCollection <br>
 * (default name empty) <br>
 *
 * @param DCH_DigiCollection The name of output collection, type extension::SenseWireHitCollection <br>
 * (default name DCH_DigiCollection) <br>
 *
 * @param DCH_DigiSimAssociationCollection The name of output association collection, type extension::SenseWireHitSimTrackerHitLinkCollection <br>
 * (default name DCH_DigiSimAssociationCollection) <br>
 *
 * @param DCH_name DCH subdetector name <br>
 * (default value DCH_v2) <br>
 *
 * @param zResolution_mm Resolution (sigma for gaussian smearing) along the sense wire, in mm <br>
 * (default value 1 mm) <br>
 *
 * @param xyResolution_mm Resolution (sigma for gaussian smearing) perpendicular the sense wire, in mm <br>
 * (default value 0.1 mm) <br>
 *
 * @param XTFileName ROOT file containing XT relation input (used by DCHXT2DLUT) <br>
 * (default value par.root) <br>
 *
 * @param EnableFFTNoise Enable FFT-based noise injection (template from ROOT file) <br>
 * (default value true) <br>
 *
 * @param NoiseFileName ROOT file containing FFT noise template <br>
 * (default value par.root) <br>
 *
 * @param NoiseDirName Directory/key name for FFT noise objects in the ROOT file <br>
 * (default value fft_noise) <br>
 *
 * @param DoWaveformDigitization Digitize waveform into ADC samples (produces RawTimeSeries L/R) <br>
 * (default value true) <br>
 *
 * @param ADCSamplePeriod_ns ADC sampling period (ns) <br>
 * (default value 0.5 ns) <br>
 *
 * @param ADCPolarity Electronics polarity (+1 positive pulse, -1 negative pulse) <br>
 * (default value +1) <br>
 *
 * @param ADCBits ADC resolution bits <br>
 * (default value 12) <br>
 *
 * @param ADCLSB_mV ADC step (LSB) in mV <br>
 * (default value 0.5 mV) <br>
 *
 * @param ADCBaseline_mV Baseline in mV <br>
 * (default value 50.0 mV) <br>
 *
 * @param create_debug_histograms Optional flag to enable QA/debug histograms registered to THistSvc <br>
 * Note: output ROOT file name is configured in runDCHdigi.py via THistSvc.Output (e.g. DCHdigi_debug_plots.root). <br>
 * (default value false) <br>
 *
 * @param GeoSvcName Geometry service name <br>
 * (default value GeoSvc) <br>
 *
 * @param uidSvcName The name of the UniqueIDGenSvc instance (reproducible seeds) <br>
 * (default value uidSvc) <br>
 */

#pragma once

// Gaudi Transformer baseclass headers
#include "Gaudi/Property.h"
#include "Gaudi/Accumulators.h"
#include "k4FWCore/Transformer.h"
#include <GaudiKernel/SmartIF.h>
#include "GaudiKernel/ITHistSvc.h"

// Gaudi services
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// EDM4HEP
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/RawTimeSeriesCollection.h"

// EDM4HEP extension
#include "extension/SenseWireHitCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"

//add extension of the waveform
//#include "extension/SenseWireAnalogWaveformCollection.h"
//#include "extension/SenseWireDigiWaveformCollection.h"

// DD4hep
#include "DD4hep/Detector.h" // for dd4hep::VolumeManager
#include "DDSegmentation/BitFieldCoder.h"

// STL
#include <random>
#include <string>

// data extension for detector DCH_v2
#include "DDRec/DCH_info.h"

// ROOT headers
#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TProfile.h"
#include "TGraph.h"
#include "DCHXT2DLUT.h"
#include <memory>


//requirments for the xtrelTIME
//class XTRELTIME; //uncommnent when needed

class TF1;
class TFile;

/// constant to convert from mm (EDM4hep) to DD4hep (cm)

struct DCHdigi_v03 final
    : k4FWCore::MultiTransformer<
          std::tuple<extension::SenseWireHitCollection, 
	  extension::SenseWireHitSimTrackerHitLinkCollection,
	  
	  //waveforms
	  edm4hep::RawTimeSeriesCollection,	//Digitized waveform left
	  edm4hep::RawTimeSeriesCollection	//Digitized waveform right
	     >(
              const edm4hep::SimTrackerHitCollection&, 
	      const edm4hep::EventHeaderCollection&)> {

  DCHdigi_v03(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<
	  extension::SenseWireHitCollection, 
	  extension::SenseWireHitSimTrackerHitLinkCollection,
	  
	  //Digi waveform collections
	  edm4hep::RawTimeSeriesCollection,
	  edm4hep::RawTimeSeriesCollection>

  operator()(const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&) const override;

private:
  /// conversion factor mm to cm, static to the class to avoid clash with DD4hep
  static constexpr double MM_TO_CM = 0.1;

  //------------------------------------------------------------------
  //          machinery for geometry

  /// Geometry service name
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
  Gaudi::Property<std::string> m_uidSvcName{this, "uidSvcName", "uidSvc", "The name of the UniqueIDGenSvc instance"};

  /// Detector name
  Gaudi::Property<std::string> m_DCH_name{this, "DCH_name", "DCH_v2", "Name of the Drift Chamber detector"};

  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;

  mutable Gaudi::Accumulators::Counter<Gaudi::Accumulators::atomicity::full, unsigned int> m_event_counter;

  /// Decoder for the cellID
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

  /// Pointer to drift chamber data extension
  dd4hep::rec::DCH_info* dch_data = {nullptr};

  double m_halfChamberLength_mm = 0.0;	// gas half-length along z

  //------------------------------------------------------------------
  //          machinery for smearing the position

  /// along the sense wire position resolution in mm
  Gaudi::Property<float> m_z_resolution{
      this, "zResolution_mm", 1.0,
      "Spatial resolution in the z direction (from reading out the wires at both sides) in mm. Default 1 mm."};
  /// xy resolution in mm
  Gaudi::Property<float> m_xy_resolution{this, "xyResolution_mm", 0.1,
                                         "Spatial resolution in the xy direction in mm. Default 0.1 mm."};

  /// create seed using the uid
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  /// Create random engine, initialized with seed out of Event Header
  std::tuple<std::mt19937_64, TRandom3> CreateRandomEngines(const edm4hep::EventHeaderCollection& headers) const;
  //------------------------------------------------------------------
  //        ancillary functions

  bool IsFileGood(std::string& ifilename) const { return std::ifstream(ifilename).good(); }

  /// Print algorithm configuration
  void PrintConfiguration(std::ostream& io);

  /// Send error message to logger and then throw exception
  void ThrowException(std::string s) const;

  int CalculateLayerFromCellID(dd4hep::DDSegmentation::CellID id) const {
    // return m_decoder->get(id, "layer") + dch_data->nlayersPerSuperlayer * m_decoder->get(id, "superlayer") + 1;
    return dch_data->CalculateILayerFromCellIDFields(m_decoder->get(id, "layer"), m_decoder->get(id, "superlayer"));
  }

  int CalculateNphiFromCellID(dd4hep::DDSegmentation::CellID id) const { return m_decoder->get(id, "nphi"); }

  TVector3 Convert_EDM4hepVector_to_TVector3(const edm4hep::Vector3d& v, double scale) const {
    return {v[0] * scale, v[1] * scale, v[2] * scale};
  };
  edm4hep::Vector3d Convert_TVector3_to_EDM4hepVector(const TVector3& v, double scale) const {
    return {v.x() * scale, v.y() * scale, v.z() * scale};
  };

  bool IsParticleCreatedInsideDriftChamber(const edm4hep::MCParticle&) const;
   
  //Mean excitation energy
  Gaudi::Property<double> m_MeanExcEnergy_eV{this, "MeanExcitationEnergy_eV", 48.48, 
	  "Mean excitation Energy I in eV for the gas (Ref: materials_o1_v02.xml)"};

  //Gas density(default: He-isobutane mixture in 90:10)
  Gaudi::Property<double> m_GasDensity_g_cm3{this, "GasDensity_g_cm3", 3.984e-4,
	  "Gas density in g/cm3 used to convert MeV.cm2/g -> MeV/cm"};

  //Energy to produce single ion pair (default: He-isobutane mixture in 90:10)
  Gaudi::Property<double> m_W_eff_eV{this, "W_eff_eV", 35.0,
	  "Effective energy per cluster in eV (W_eff)"};

  //Mass of the MC particle
  //Gaudi::Property<double> m_MassForBB_GeV{this, "MassForBB_GeV", 0.105658,
//	  "Reference particle mass (GeV) for BetaGamma conversion (default muon mass)"};

  //function to calculate cluster density using Betagamma
  double get_dNcldx_per_cm(double betagamma, const edm4hep::MCParticle& mc) const;

  //// ROOT file containing the x-t relation
  Gaudi::Property<std::string> m_xtFileName{this, "XTFileName", "par.root",
	  "ROOT file with x-t relation used to convert radius to drift time" };

  /************************************************************
   * This block of code is useful for using r vs t
   ************************************************************
  mutable XTRELTIME* m_xtHelper{nullptr};
  double electronDriftTime(double r_cm, TRandom3& myRandom) const;
  double electronDriftVelocity_cm_per_us(double r_cm) const;
  double electricField_V_per_cm(double r_cm) const;
  *************************************************************/

  //pointer to the LUT
  std::unique_ptr<DCHXT2DLUT> m_xt2d;

  //Gaudi properties for the Gas Gain/ polya parameters:
  Gaudi::Property<double> m_GasGain{this, "GasGain", 2.0e5,
	  "Mean gas gain used in Polya distribution (arbitrary units)" };
  
  //Polya distribution shape parameter
  Gaudi::Property<double> m_PolyaTheta{this, "PolyaTheta", 0.5,
	  "Polya shape parameter theta (typicaly around 5 for drift tube)"};
  
  TF1* m_polya{nullptr};

  //Function to produce avalunche charge sample from polya
  double avalancheCharge(TRandom3& myRandom) const;

  //------------------------- pulse shaping block -------------------------//

  //Scale parameter
  Gaudi::Property<double> m_pulseAmplitudeScale{this, "PulseAmplitudeScale", 5.e-6,
          "Global Scale factor converting avalache charge to voltage"};
  
  //Pulse rise time
  Gaudi::Property<double> m_pulseRiseTime_ns{this, "PulseRiseTime_ns", 1.0,
          "Rise time of the single electron pulse (ns)"};

  //Pulse decay time
  Gaudi::Property<double> m_pulseFallTime_ns{this, "PulseFallTime_ns", 7.0,
          "Fall time of the single electron pulse (ns)"};

  //Pulse bin size
  Gaudi::Property<double> m_pulseBinSize_ns{this, "PulseBinSize_ns", 1.0,
          "Time bin size for the analog pulse (ns)"};

    // Time window used only for the debug single-electron pulse histogram
  Gaudi::Property<double> m_signalPulseWindow_ns{this, "SignalPulseWindow_ns", 500.0,
          "Displayed time window of the debug single-electron pulse histogram (ns)"};

  //Analog waveform window
  Gaudi::Property<double> m_waveformTimeWindow_ns{this, "WaveformTimeWindow_ns", 3000.0,
          "Time window of analog waveform (ns)"};

  //Analog waveform bin size
  Gaudi::Property<double> m_waveformBinSize_ns{this, "WaveformBinSize_ns", 1.0,
          "Time bin size for the analog waveform (ns)"};

  //Function used to produce Single Pulse
  double singleElectronPulse(double t_ns, double t0_ns, double q) const;

  // --------------- Signal propagation along sense wire  ----------------------//
  
  //Signal speed along the wire
  Gaudi::Property<double> m_signalSpeed_mm_per_ns{this, "signalSpeed_mm_per_ns", 200.0,
	  "Propagation speed along sense wire [mm/ns]. 200 mm/ns ~ 5 ns/m"};

  //Enable or disable the attenuation
  Gaudi::Property<bool> m_enableWireAttenuation{this, "EnableWireAttenuation", true,
          "Enable signal attenuation along the sense wire (charge-division effect)"};

  //Attenuation in signal
  Gaudi::Property<double> m_wireAttenuationLength_mm{this, "WireAttenuationLength_mm", 1.0e4,
	  "Attenuation length lambda [mm] for exp(-d/lambda). Very large disables attenuation"};

  // ---------- impedence mismatch + electronics transfer fucntion: -----------//

  //Front end gain of the electronics
  Gaudi::Property<double> m_frontEndGain{this, "FrontEndGain", 1.0,
	  "Overall analog gain applied before ADC (dimensionless)"};

  //Apply mismatch
  Gaudi::Property<bool> m_enableImpedanceMismatch{this, "EnableImpedanceMismatch", true,
	  "Enable impedance mismatch scaling (reflection/transmission)"};

  //Set the value of mismatch resistance
  Gaudi::Property<double> m_matchingRes_Ohm{this, "MatchingRes_Ohm", 330.0,
	  "Matching / termination resistance at preamp input [Ohm]"};

  //Set the value of the Tube Impedance
  Gaudi::Property<double> m_tubeImpedance_Ohm{this, "TubeImpedance_Ohm", 50.0,
	  "Transmission line / tube characteristic impedance [Ohm]"};

  //Apply Electronic transfer function
  Gaudi::Property<bool> m_enableElectronicsTF{this, "EnableElectronicsTF", true,
	  "Enable electronics transfer function (shaper impulse response)"};

  //Set the value of Time rise
  Gaudi::Property<double> m_elecTauRise_ns{this, "ElectronicsTauRise_ns", 3.0,
	  "Electronics shaping rise time constant (ns)"};

  //Set the value of Time fall1
  Gaudi::Property<double> m_elecTauFall1_ns{this, "ElectronicsTauFall1_ns", 9.0,
          "Electronics shaping fall time constant #1 (ns)"};

  //Set the value of Time fall2
  Gaudi::Property<double> m_elecTauFall2_ns{this, "ElectronicsTauFall2_ns", 25.0,
          "Electronics shaping fall time constant #2 (ns)"};

  //Set the fraction value of the Time fall1 and fall2
  Gaudi::Property<double> m_elecMixFraction{this, "ElectronicsMixFraction", 0.5,
	  "Mixing fraction between fall1 and fall2 (0..1)"};

  // Length (ns) of the electronics impulse response used in waveform convolution.
  Gaudi::Property<double> m_elecKernelLength_ns{this, "ElectronicsKernelLength_ns", 200.0,
	  "Length of the impulse response kernel used for discrete convolution (ns)"};

  //--------- FFT/IFFT colored noise (from par.root/fft_noise) --------------//
  
  // Enable addition of frequency-shaped (colored) electronic noise using FFT/IFFT.
  Gaudi::Property<bool> m_enableFFTNoise{this, "EnableFFTNoise", true,
	  "Enable adding colored noise using FFT magnitude from par.root"};
  
  // ROOT file containing the FFT-based noise frequency template.
  Gaudi::Property<std::string> m_noiseFileName{this, "NoiseFileName", "par.root",
	  "ROOT file that contains directory fft_noise with fft_freq/fft_mag/fft_amp"};
  
  // Directory name inside the ROOT file holding the FFT noise histograms.
  Gaudi::Property<std::string> m_noiseDirName{this, "NoiseDirName", "fft_noise",
	  "Directory name inside ROOT file containing fft_freq/fft_mag/fft_amp"};
  
  // Global scaling factor applied to the generated time-domain noise amplitude.
  Gaudi::Property<double> m_noiseScale{this, "NoiseScale", 1.0e-3,
	  "Overall scale factor for FFT noise (multiplies the generated time-noise)"};
  
  // Remove DC (baseline) component from the generated noise waveform.
  Gaudi::Property<bool> m_noiseRemoveDC{this, "NoiseRemoveDC", true,
	  "Remove DC offset from generated noise (recommended)"};
  // Cached FFT-noise template (loaded once in initialize)
  mutable bool m_fftNoiseReady{false};
  mutable std::vector<double> m_fftMag;
  mutable double m_fftAmpNorm{1.0};     // fft_amp (normalization from file)
  mutable double m_fftMaxFreq{0.0};     // max frequency from fft_freq (same units as 1/ns)
  mutable int    m_fftSize{0};          // size of m_fftMag

  // Enable conversion of the analog waveform into a digitized ADC waveform.
  Gaudi::Property<bool>   m_doWaveformDigitization{this, "DoWaveformDigitization", true,
  "Digitize analog waveform into ADC time bins"};

  // ADC sampling period for digitized waveform
  Gaudi::Property<double> m_adcSamplePeriod_ns{this, "ADCSamplePeriod_ns", 0.5,
	  "ADC sampling period (ns) for digitized waveform"};
  
  // Set electronics polarity: +1 keeps pulse sign, -1 inverts it.
  Gaudi::Property<int> m_adcPolarity{this, "ADCPolarity", +1,
  "+1: positive pulse, -1: negative pulse (electronics polarity)"};
  
  // Number of bits defining the ADC amplitude resolution.
  Gaudi::Property<int>    m_adcBits{this, "ADCBits",12, 
	  "ADC resolution bits"};
  
  // Voltage step size (mV) corresponding to one ADC count
  Gaudi::Property<double> m_adcLSB_mV{this, "ADCLSB_mV", 0.5, 
	  "ADC step (LSB) in mV"};
  
  // Baseline voltage (mV) added to the digitized waveform
  Gaudi::Property<double> m_adcBaseline_mV{this, "ADCBaseline_mV", 50.0, 
	  "Baseline in mV"};

  // Enable diagnostic histograms registered to THistSvc
  Gaudi::Property<bool> m_create_debug_histos{this, "create_debug_histograms", false,
                                              "Enable diagnostic histograms via THistSvc"};

  /// histogram to store distance from sim hit position to the sense wire
  SmartIF<ITHistSvc> m_histSvc;

  TH1D* hDpw{nullptr};

  /// histogram to store distance from digi-hit to the wire. Should be zero because digi-hit central position lies on
  /// the wire. This histogram is a consistency check, because the function used to calculate the distance to the wire
  /// is different from the function used to calculate the digi-hit central position from a sim-hit position
  TH1D* hDww{nullptr};

  /// histogram to store smearing along the wire
  TH1D* hSz{nullptr};

  /// histogram to store smearing perpendicular the wire
  TH1D* hSxy{nullptr};

  TH1F* hPathLength{nullptr};
  TH2F* hPLvsnC{nullptr};
  TH1F* hX{nullptr};
  TH1F* hY{nullptr};
  TH1F* hZ{nullptr};
  TH3F* hXYZ{nullptr};
  TH1F* heDep{nullptr};
  TH1F* hR{nullptr};
  TH2F* hRvsZ{nullptr};
  TH1F* hTotPathCell{nullptr};
  TH1F* hBetaGammaCell{nullptr};

  TH1F* hNcl{nullptr};
  TH1F* hNe{nullptr};
  TH1F* hClSpacing_mm{nullptr};
  
  TH1F* hTd{nullptr};
  TProfile* hXT{nullptr};

  TH1F* hAvalancheQ{nullptr};

  mutable TH1F* hSignalPulse;
  TH1F* hWaveform{nullptr};
  TH1F* hWaveformL{nullptr};
  TH1F* hWaveformR{nullptr};
  TH1F* hWaveformElecL{nullptr};
  TH1F* hWaveformElecR{nullptr};
  TH1F* hNoise{nullptr};
  TH1F* hWaveformDigiL{nullptr};
  TH1F* hWaveformDigiR{nullptr};

  mutable bool m_waveformFilled{false};

};

DECLARE_COMPONENT(DCHdigi_v03);


