#include "DCHdigi_v03.h"
#include "BetheBloch.h"
#include "DCHUtilities.h"
#include "DCHFFTNoise.h"
#include "DCHXT2DLUT.h"

#include "TMath.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TDirectory.h"
#include <numeric> 

// STL
#include <algorithm>
#include <iostream>
#include <sstream>
#include <random>
#include <cmath>
#include <unordered_map>
#include "extension/MutableSenseWireHit.h"

// Gaudi
#include <GaudiKernel/MsgStream.h>

//header for the waveform
//#include "extension/SenseWireAnalogWaveformCollection.h"
//#include "extension/SenseWireDigiWaveformCollection.h"
//#include "extension/MutableSenseWireAnalogWaveform.h"
//#include "extension/MutableSenseWireDigiWaveform.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"


///////////////////////////////////////////////////////////////////////////////////////
//////////////////////       DCHdigi_v03 constructor       ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// -- KeyValues("name of the variable that holds the name of the collection exposed in the python steering file",
// {"default name for the collection"}),
DCHdigi_v03::DCHdigi_v03(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {
		       	   //input collection
                           KeyValues("DCH_simhits", {""}),
                           KeyValues("HeaderName", {"EventHeader"}),
                       },
                       {
		       //output collection
		       KeyValues("DCH_DigiCollection", {"DCH_DigiCollection"}),
                        KeyValues("DCH_DigiSimAssociationCollection", {"DCH_DigiSimAssociationCollection"}),
			// hWaveform
			//KeyValues("DCH_AnalogWaveformCollection", {"DCH_AnalogWaveforms"}),
			
			// hWaveformL
			//KeyValues("DCH_AnalogWaveformLCollection", {"DCH_AnalogWaveformsL"}),
			// hWaveformR
			//KeyValues("DCH_AnalogWaveformRCollection", {"DCH_AnalogWaveformsR"}),
			
			// hWaveformDigiL
			KeyValues("DCH_DigiWaveformLCollection",   {"DCH_DigiWaveformsL"}),
			// hWaveformDigiR
			KeyValues("DCH_DigiWaveformRCollection",   {"DCH_DigiWaveformsR"})
			}) {
  m_geoSvc = serviceLocator()->service(m_geoSvcName);	// Access detector geometry service
  m_uidSvc = serviceLocator()->service(m_uidSvcName);	// Access reproducible unique-ID seed service
}

// Random number generator
std::tuple<std::mt19937_64, TRandom3> DCHdigi_v03::CreateRandomEngines(
                const edm4hep::EventHeaderCollection& headers) const {
  auto engine_seed = m_uidSvc->getUniqueID(headers, this->name());
  auto rng_engine = std::mt19937_64(engine_seed);
  auto random_seed = m_uidSvc->getUniqueID(headers, this->name()+"_1");
  auto myRandom = TRandom3(random_seed);
  // advance internal state to minimize possibility of creating correlations
  rng_engine.discard(10);
  for (int i = 0; i < 10; ++i)
    myRandom.Rndm();
  return {rng_engine, myRandom};
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       initialize       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode DCHdigi_v03::initialize() 
{
  if (!m_uidSvc)
    ThrowException("Unable to get UniqueIDGenSvc");

  if (0 > m_z_resolution.value())
    ThrowException("Z resolution input value can not be negative!");

  if (0 > m_xy_resolution.value())
    ThrowException("Radial (XY) resolution input value can not be negative!");

  // Retrieve the subdetector
  std::string DCH_name(m_DCH_name.value());
  if (0 == m_geoSvc->getDetector()->detectors().count(DCH_name)) {
    ThrowException("Detector <<" + DCH_name + ">> does not exist.");
  }
  dd4hep::DetElement DCH_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);

  //Event counter
  m_event_counter.reset();

  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////  retrieve data extension     //////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  this->dch_data = DCH_DE.extension<dd4hep::rec::DCH_info>();
  if (not dch_data->IsValid())
    ThrowException("No valid data extension was found for detector <<" + DCH_name + ">>.");

  ///////////////////////////////////////////////////////////////////////////////////
	
  // Retrieve the readout associated with the detector element (subdetector)
  dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector(DCH_name);
  if (not dch_sd.isValid())
    ThrowException("No valid Sensitive Detector was found for detector <<" + DCH_name + ">>.");

  dd4hep::Readout dch_readout = dch_sd.readout();
  // set the cellID decoder
  m_decoder = dch_readout.idSpec().decoder();

  // Read DCH half-length from XML constant (DCH_standalone_o1_v02.xml)
  dd4hep::Detector& det = dd4hep::Detector::getInstance();

  if (!det.constantAsDouble("DCH_gas_Lhalf")) {
    error() << "Missing DD4hep constant: DCH_gas_Lhalf" << endmsg;
    return StatusCode::FAILURE;
  }

  m_halfChamberLength_mm = det.constantAsDouble("DCH_gas_Lhalf") / MM_TO_CM; // returns in mm
  info() << "[DCHdigi] m_halfChamberLength_mm = " << m_halfChamberLength_mm << " mm" << endmsg;


  ///////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////  initialize x-t relation  //////////////////////////
  /////////////////////////////////////////////////////////////////////////////////// 
  // Load the 2D x-t lookup table used to convert transverse cluster position into drift time.
  m_xt2d = std::make_unique<DCHXT2DLUT>();	// Create the x-t LUT helper object.
  
  // Load mean and sigma
  const bool xt_file = m_xt2d->load(m_xtFileName.value(), "xt_mean", "xt_error");
  if (!xt_file) {
	  error() << "[XT2D] Failed to load/build 2D XT LUT from " << m_xtFileName.value()
		  << " (need TGraph2D xt_mean and xt_error)" << endmsg;
	  return StatusCode::FAILURE;
  }

  info() << "[XT2D] Built regular-grid XT LUT from " << m_xtFileName.value()<< endmsg;

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////        Initialize Noise    ////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  // Load FFT-noise template
  if (m_enableFFTNoise.value()) {
	  std::string err;
	  // Load FFT magnitude template + normalization + frequency
	  const bool noisefile = dchfft::loadFFTNoiseTemplate(
		    	  m_noiseFileName.value(),	// ROOT file containing the FFT-noise template
		    	  m_noiseDirName.value(),	

		    	  m_fftMag, 	//Magnitude
			  m_fftAmpNorm,
			  m_fftMaxFreq,
			  m_fftSize,
		    	  &err
			  );

	  if (!noisefile) {
	      	  error() << "[FFTNoise] " << err << endmsg;
	      	  return StatusCode::FAILURE;
	  }

	  // Mark noise system ready for waveform generation
	  m_fftNoiseReady = true;

	  info() << "[FFTNoise] Loaded fft_noise:"
	 	  << " size=" << m_fftSize
	 	  << " maxFreq=" << m_fftMaxFreq
	 	  << " ampNorm=" << m_fftAmpNorm << endmsg;
  }

  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////  initialize Polya gain  ///////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  // Polya PDF for gas gain:
  // P(Q) ~ ((1+theta)/(G * Gamma(theta+1))) *
  //        [Q*(1+theta)/G]^theta * exp( - Q*(1+theta)/G )
  //
  // Parameters:
  //   [0] normalization (kept at 1; ROOT will renormalize as needed)
  //   [1] mean gas gain G
  //   [2] theta (Polya shape parameter)
  //
  // Range (1e3, 1e7) is taken from the original GMCT implementation
  // and can be tuned later if needed.
  m_polya = new TF1("Polya",
		  "((1.+[2])/([1]*TMath::Gamma([2]+1))) * pow(x*(1.+[2])/[1],[2]) * exp(-x*(1.+[2])/[1])",
		  1.e3, 1.e7);
  m_polya->SetParameter(0, 1.0);
  m_polya->SetParameter(1, m_GasGain.value());
  m_polya->SetParameter(2, m_PolyaTheta.value());
  //-----------------------------------------------------------------------------------//

  std::stringstream ss;
  PrintConfiguration(ss);
  info() << ss.str().c_str() << endmsg;

  // Create histograms for debugging:
  if (m_create_debug_histos.value()) 
  {
    m_histSvc = service("THistSvc");
    if (!m_histSvc) {
      error() << "Unable to locate Histogram Service" << endmsg;
      return StatusCode::FAILURE;
    }
    hDpw = new TH1D("hDpw", "Distance from sim-hit to the wire, in cm", 100, 0, 1);
    if (m_histSvc->regHist("/rec/hDpw", hDpw).isFailure()) {
            error() << "Couldn't register hDpw histogram" << endmsg;
    }
    hDww = new TH1D(
        "hDww",
        "Distance from digi-hit to the wire, in cm. Should be zero",
        100, 0, 1);
    if (m_histSvc->regHist("/rec/hDww", hDww).isFailure()) {
            error() << "Couldn't register hDww histogram" << endmsg;
    }
    hSz = new TH1D("hSz", "Smearing along the wire, in cm", 100, 0, 5 * m_z_resolution.value());
    if (m_histSvc->regHist("/rec/hSz", hSz).isFailure()) {
            error() << "Couldn't register hSz histogram" << endmsg;
    }
    hSxy = new TH1D("hSxy", "Smearing perpendicular the wire, in cm", 100, 0, 5 * m_xy_resolution.value());
    if (m_histSvc->regHist("/rec/hSxy", hSxy).isFailure()) {
            error() << "Couldn't register hSxy histogram" << endmsg;
    }

  hPathLength = new TH1F("hPathLength", "Path Length per Hit;Path Length [mm];Entries", 10, 0, 15);
  if (m_histSvc->regHist("/rec/hPathLength", hPathLength).isFailure()) {
            error() << "Couldn't register hPathLength histogram" << endmsg;
    }
  hPLvsnC = new TH2F("hPLvsnC","Ionized Length vs nCell; PathLength [mm]; nCell", 20, 0, 15,  20, 0, 200);
  if (m_histSvc->regHist("/rec/hPLvsnC", hPLvsnC).isFailure()) {
            error() << "Couldn't register hPLvsnC histogram" << endmsg;
    }

  hX = new TH1F("hX", "Hit X positions;X (mm);Counts", 50, 0, 2000);
  if (m_histSvc->regHist("/rec/hX", hX).isFailure()) {
            error() << "Couldn't register hX histogram" << endmsg;
    }
  hY = new TH1F("hY", "Hit Y positions;Y (mm);Counts", 50, 0, 2000);
  if (m_histSvc->regHist("/rec/hY", hY).isFailure()) {
            error() << "Couldn't register hY histogram" << endmsg;
    }
  hZ = new TH1F("hZ", "Hit Z positions;Z (mm);Counts", 50, 0, 2000);
  if (m_histSvc->regHist("/rec/hZ", hZ).isFailure()) {
            error() << "Couldn't register hZ histogram" << endmsg;
    }
  hXYZ = new TH3F("hXYZ", "Hit positions;X (mm);Y (mm);Z (mm)", 50, 0, 2000, 50, 0, 2000, 50, 0, 2000);
  if (m_histSvc->regHist("/rec/hXYZ", hXYZ).isFailure()) {
            error() << "Couldn't register hXYZ histogram" << endmsg;
    }
  heDep = new TH1F("heDep", "Energy Deposit", 50, 0, 10);
  if (m_histSvc->regHist("/rec/heDep", heDep).isFailure()) {
            error() << "Couldn't register heDep histogram" << endmsg;
    }
  hR = new TH1F("hR", "Distance(mm)", 50, 0, 2000);
  if (m_histSvc->regHist("/rec/hR", hR).isFailure()) {
            error() << "Couldn't register hR histogram" << endmsg;
    }
  hRvsZ = new TH2F("hRvsZ","R vs z; z [mm]; R [mm]", 50, 0, 2400,  50, 0, 2400);
  if (m_histSvc->regHist("/rec/hRvsZ", hRvsZ).isFailure()) {
            error() << "Couldn't register hRvsZ histogram" << endmsg;
    }
  hNcl = new TH1F("hNcl", "n_{cl} per cm; n_{cl}; entries", 40, 0, 50);
  if (m_histSvc->regHist("/rec/hNcl", hNcl).isFailure()) {
            error() << "Couldn't register hNcl histogram" << endmsg;
    }
  hNe = new TH1F("hNe", "Cluster Size Distribution; Cluster size; entries", 10, 0, 10);
  if (m_histSvc->regHist("/rec/hNe", hNe).isFailure()) {
            error() << "Couldn't register hNe histogram" << endmsg;
    }
  hClSpacing_mm = new TH1F("hClSpacing_mm", "Inter-cluster spacing; spacing [mm]; entries", 50, 0, 6);
  if (m_histSvc->regHist("/rec/hClSpacing_mm", hClSpacing_mm).isFailure()) {
            error() << "Couldn't register hClSpacing_mm histogram" << endmsg;
    }
  //-------------------------------------------------------------------------------------------//

  hTd = new TH1F  ("hTd", "Drift Time; t (ns); Entries", 50, 0, 1000);
  if (m_histSvc->regHist("/rec/hTd", hTd).isFailure()) {
            error() << "Couldn't register hTd histogram" << endmsg;
    }
  hXT = new TProfile ("hXT", "x-t relation; r (cm); t (ns)", 50, 0, 1.0);
  if (m_histSvc->regHist("/rec/hXT", hXT).isFailure()) {
            error() << "Couldn't register hXT histogram" << endmsg;
    }
  
  hAvalancheQ = new TH1F("hAvalancheQ", "Avalanche per Charge; Q  (a.u.); Entries", 100, 0, 150.e4);
  if (m_histSvc->regHist("/rec/hAvalancheQ", hAvalancheQ).isFailure()) {
            error() << "Couldn't register hAvalancheQ histogram" << endmsg;
    }
  //-------------------------------------------------------------------------------------------//

  const double pulseTmin = 0.0;
  const double pulseTmax = m_signalPulseWindow_ns.value();
  const double dtPulse_ns = m_pulseBinSize_ns.value();
  const int nBinsSignalPulse = std::max(
      1, static_cast<int>(std::lround((pulseTmax - pulseTmin) / dtPulse_ns)));

  hSignalPulse = new TH1F("hSignalPulse",
                          "Single-electron pulse; t [ns]; amplitude [V]",
                          nBinsSignalPulse, pulseTmin, pulseTmax);
  if (m_histSvc->regHist("/rec/hSignalPulse", hSignalPulse).isFailure()) {
    error() << "Couldn't register hSignalPulse histogram" << endmsg;
  }

  const double wfTmin = 0.0;
  const double wfTmax = m_waveformTimeWindow_ns.value();
  
  const double dtAnalog_ns = m_waveformBinSize_ns.value();
  const double dtDigi_ns = m_adcSamplePeriod_ns.value();


  const int nBinsAnalogWf  = static_cast<int>(std::lround((wfTmax - wfTmin) / dtAnalog_ns));
  const int nBinsDigiWf = static_cast<int>(std::lround((wfTmax - wfTmin)/dtDigi_ns));

  hWaveform = new TH1F("hWaveform", "Analog Waveform; t [ns]; Amplitude", nBinsAnalogWf, wfTmin, wfTmax);
  if (m_histSvc->regHist("/rec/hWaveform", hWaveform).isFailure()) {
            error() << "Couldn't register hWaveform histogram" << endmsg;
    }
  hWaveformL = new TH1F("hWaveformL", "Waveform at wire end 1; t [ns]; Amplitude", nBinsAnalogWf, wfTmin, wfTmax);
  if (m_histSvc->regHist("/rec/hWaveformL", hWaveformL).isFailure()) {
            error() << "Couldn't register hDpw histogram" << endmsg;
    }
  hWaveformR = new TH1F("hWaveformR","Waveform at wire end 2; t [ns]; Amplitude",nBinsAnalogWf, wfTmin, wfTmax);
  if (m_histSvc->regHist("/rec/hWaveformR", hWaveformR).isFailure()) {
            error() << "Couldn't register hWaveformR histogram" << endmsg;
    }
  //-------------------------------------------------------------------------------------------//

  hWaveformElecL = new TH1F("hWaveformElecL", "Waveform after electronics (end 1); t [ns]; Amplitude",
		  nBinsAnalogWf, wfTmin, wfTmax);
  if (m_histSvc->regHist("/rec/hWaveformElecL", hWaveformElecL).isFailure()) {
            error() << "Couldn't register hWaveformElecL histogram" << endmsg;
    }
  hWaveformElecR = new TH1F("hWaveformElecR","Waveform after electronics (end 2); t [ns]; Amplitude",
                          nBinsAnalogWf, wfTmin, wfTmax);
  if (m_histSvc->regHist("/rec/hWaveformElecR", hWaveformElecR).isFailure()) {
            error() << "Couldn't register hWaveformElecR histogram" << endmsg;
    }
  //------------------------------------------------------------------------------------------//

  // Debug histogram for FFT noise (single example waveform)
  hNoise = new TH1F("hNoise", "FFT Noise; t [ns]; Amplitude",
		  nBinsAnalogWf, wfTmin, wfTmax);
  if (m_histSvc->regHist("/rec/hNoise", hNoise).isFailure()) {
	  error() << "Couldn't register hNoise histogram" << endmsg;
  }

  // Digitized (reconstructed) waveform in mV
  hWaveformDigiL = new TH1F("hWaveformDigiL","Digitized waveform at end 1; t [ns]; voltage [mV]",
                          nBinsDigiWf, wfTmin, wfTmax);
  if (m_histSvc->regHist("/rec/hWaveformDigiL", hWaveformDigiL).isFailure()) {
            error() << "Couldn't register hWaveformDigiL histogram" << endmsg;
    }
  hWaveformDigiR = new TH1F("hWaveformDigiR","Digitized waveform at end 2; t [ns]; voltage [mV]",
                          nBinsDigiWf, wfTmin, wfTmax);
  if (m_histSvc->regHist("/rec/hWaveformDigiR", hWaveformDigiR).isFailure()) {
            error() << "Couldn't register hWaveformDigiR histogram" << endmsg;
    }
  } //End of histogram definition
  //----------------------------------------------------------------------------------------//

  return StatusCode::SUCCESS;

} //end of Initialize

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       operator()       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

std::tuple<extension::SenseWireHitCollection,
	extension::SenseWireHitSimTrackerHitLinkCollection,
	
	//use raw time series collection:
	edm4hep::RawTimeSeriesCollection,
	edm4hep::RawTimeSeriesCollection>

DCHdigi_v03::operator()(const edm4hep::SimTrackerHitCollection& input_sim_hits,
                    const edm4hep::EventHeaderCollection&   headers) const 
{

  if (m_event_counter.value()%1 == 0)
  {
	  info() << "Processing Events: " << m_event_counter.value() <<endmsg;
  }

  /// initialize engines
  auto [rng_engine, myRandom] = this->CreateRandomEngines(headers);

  // Gaussian random number generator used for the smearing of the z position, in cm!
  std::normal_distribution<double> gauss_z_cm{0., m_z_resolution.value() * MM_TO_CM};
  // Gaussian random number generator used for the smearing of the xy position, in cm!
  std::normal_distribution<double> gauss_xy_cm{0., m_xy_resolution.value() * MM_TO_CM};

  // Create the collections we are going to return
  extension::SenseWireHitCollection output_digi_hits;
  extension::SenseWireHitSimTrackerHitLinkCollection output_digi_sim_association;

  //creating output collection for waveform Digi
  edm4hep::RawTimeSeriesCollection out_Digitized_waveformsL;
  edm4hep::RawTimeSeriesCollection out_Digitized_waveformsR;

  //Varible initailizations:
  edm4hep::Vector3d pos;
  edm4hep::MCParticle mcParticle;
  float pathLength;
  float distance;
  double energyDep;
  int ilayer;
  int nphi;

  // loop over hit collection
  int loop_index = 0;

  // Store only one waveform for debuging
  m_waveformFilled = false;
  // Local analog waveform definition used for ALL cells
  const double wfTmin = 0.0;
  const double wfTmax = m_waveformTimeWindow_ns.value();
  const double dtAnalog_ns = m_waveformBinSize_ns.value();

  const int nBinsAnalogWf =
      std::max(1, static_cast<int>(std::lround((wfTmax - wfTmin) / dtAnalog_ns)));

  std::vector<double> analogTimeAxis(nBinsAnalogWf, 0.0);
  for (int i = 0; i < nBinsAnalogWf; ++i) {
    analogTimeAxis[i] = wfTmin + (i + 0.5) * dtAnalog_ns;
  }

  // Count how many input SimTrackerHits belong to the same cell in this event
  std::unordered_map<uint64_t, int> simHitCountPerCell;
  
  // Per-cell accumulation of electrons
  struct CellAccum {
    std::vector<double> times_ns;	// Absolute electron arrival times
    std::vector<double> charges;	// Avalanche charge per electron
    std::vector<double> wireZ_mm;	// z position of each electron production point [mm]
    std::vector<uint16_t> electronsPerCluster;		// number of electron per cluster
    std::vector<edm4hep::SimTrackerHit> simHits;	// Keep all contributing sim hits

    // Store representative metadata once per cell
    bool hasMeta = false;
    double representativeSimTime_ns = std::numeric_limits<double>::max();
    edm4hep::Vector3d positionSW{};
    float distanceToWire = 0.f;
    float wireAzimuthalAngle = 0.f;
    float wireStereoAngle = 0.f;
    std::int32_t type = 0;
    std::int32_t quality = 0;
    float eDep = 0.f;
    float eDepError = 0.f;
  };
  std::unordered_map<uint64_t, CellAccum> cellAccumMap;

  // start on each simulated hit(GEANT4 STEP SEGMENT):
  for (const auto& input_sim_hit : input_sim_hits) 
  {
    loop_index++;

    pos = input_sim_hit.getPosition();
    energyDep = input_sim_hit.getEDep();
    pathLength = input_sim_hit.getPathLength();
    mcParticle = input_sim_hit.getParticle();

    dd4hep::DDSegmentation::CellID cellid = input_sim_hit.getCellID();

    simHitCountPerCell[static_cast<uint64_t>(cellid)]++;
    ilayer = this->CalculateLayerFromCellID(cellid);
    nphi = this->CalculateNphiFromCellID(cellid);
    const int hit_ilayer = ilayer;
    const int hit_nphi   = nphi;

    auto hit_position = Convert_EDM4hepVector_to_TVector3(input_sim_hit.getPosition(), MM_TO_CM);
    distance = std::sqrt(pos.x*pos.x + pos.y*pos.y);

    //Looking for only MC Particle Hit:
    const bool isPrimaryHit = 
	    mcParticle.getParents().empty() &&
	    !input_sim_hit.isProducedBySecondary();

    // Looking for secondray hits:
    bool isDeltaElectronHit =(std::abs(mcParticle.getPDG()) == 11) &&
	    input_sim_hit.isProducedBySecondary() &&
	    IsParticleCreatedInsideDriftChamber(mcParticle);
	    
    if (!isPrimaryHit) continue;
    if (isDeltaElectronHit) continue;
 
    //Print some variable
    debug() << "Hit #" << loop_index
	    << " layer=" << ilayer
	    << " cell=" << nphi
	    << " pathLength_mm=" << input_sim_hit.getPathLength()
	    << endmsg;

    //filling histograms
    if (m_create_debug_histos.value()) {

    	hPathLength->Fill(pathLength);
    	hPLvsnC->Fill(pathLength, nphi);
    	hX->Fill(pos.x);
    	hY->Fill(pos.y);
    	hZ->Fill(pos.z);
    	hXYZ->Fill(pos.x, pos.y, pos.z);
    	heDep->Fill(energyDep);
    	hR->Fill(distance);
    	hRvsZ->Fill(pos.z, distance);
    }
//=======================================================================================
//                                  Calculate N_cluster                                //
//=======================================================================================
    
    const double L_mm = input_sim_hit.getPathLength();
    //const double ave_length = L_mm/nphi;

    if (L_mm <=0) continue;
    const double l_cm = L_mm*MM_TO_CM;

    //compute beta gamma
    const double bg = compute_beta_gamma(mcParticle);

    //mean value of the number of eletrons per cluster
    static const double Ne_mean = mean_cluster_size_He();
    const double lambda = get_dNcldx_per_cm(bg, mcParticle)/Ne_mean;
    double mu = lambda * l_cm;
    
    //Generate Cluster modeled by Zero Trancated poisson:
    const int Ncl_step = ZT_poisson(mu, rng_engine);
    const double Ncl_per_cm = (l_cm > 0.) ? 
	    static_cast<double>(Ncl_step)/l_cm : 0;
    if (m_create_debug_histos.value()) {
	    hNcl->Fill(Ncl_per_cm);
    }

    //Generate cluster position inside the loop:
    auto cl_positions_mm = generate_cluster_positions_mm(L_mm, lambda, rng_engine);
    for (size_t i = 1; i < cl_positions_mm.size(); ++i) {
	    if (m_create_debug_histos.value()) {
		    hClSpacing_mm->Fill(cl_positions_mm[i] - cl_positions_mm[i-1]);
	    }
    }

    //Generate Cluster size inside the loop:
    int n_clusters_in_step = cl_positions_mm.size();
    int Ne;
    std::vector<int> electrons_per_cluster;

    //loop over the number of clusters
    for (int i = 0; i < n_clusters_in_step; ++i) {
	    Ne = sample_cluster_size(rng_engine);
	    electrons_per_cluster.push_back(Ne); 
    }

    //loop over the number of electrons per cluster
    for (auto ne : electrons_per_cluster) {
	    if (m_create_debug_histos.value()) {
		    hNe->Fill(ne);
	    }
    }
    //-------------------------------------------------------------------------

    //calculate hit position projection into the wire
    TVector3 hit_to_wire_vector = this->dch_data->Calculate_hitpos_to_wire_vector(
		    ilayer, nphi, hit_position);

    TVector3 hit_projection_on_the_wire = hit_position + hit_to_wire_vector;
    if (m_create_debug_histos.value()) {
      double distance_hit_wire = hit_to_wire_vector.Mag();
      hDpw->Fill(distance_hit_wire);
    }
    TVector3 wire_direction_ez = this->dch_data->Calculate_wire_vector_ez(ilayer, nphi);

    // -------------------------------------------------------------------------
    // smear position along the wire
    double smearing_z = gauss_z_cm(rng_engine);
    if (m_create_debug_histos.value())
      hSz->Fill(smearing_z);

    hit_projection_on_the_wire += smearing_z * (wire_direction_ez.Unit());
    if (m_create_debug_histos.value()) {
      // the distance from the hit projection and the wire should be zero
      TVector3 dummy_vector = this->dch_data->Calculate_hitpos_to_wire_vector(
		      ilayer, nphi, hit_projection_on_the_wire);
      hDww->Fill(dummy_vector.Mag());
    }

    // smear position perpendicular to the wire
    double smearing_xy = gauss_xy_cm(rng_engine);
    if (m_create_debug_histos.value())
      hSxy->Fill(smearing_xy);
    float distanceToWire_real = hit_to_wire_vector.Mag();

    // protect against negative values
    float distanceToWire_smeared = std::max(0.0, distanceToWire_real + smearing_xy);

    std::int32_t type = 0;
    std::int32_t quality = 0;
    float eDepError = 0;
    // length units back to mm
    auto positionSW = Convert_TVector3_to_EDM4hepVector(hit_projection_on_the_wire, 1. / MM_TO_CM);
    // auto  directionSW    = Convert_TVector3_to_EDM4hepVector(wire_direction_ez, 1. / MM_TO_CM);
    float distanceToWire = distanceToWire_smeared / MM_TO_CM;

    // The direction of the sense wires can be calculated as:
    //   RotationZ(WireAzimuthalAngle) * RotationX(stereoangle)
    // One point of the wire is for example the following:
    //   RotationZ(WireAzimuthalAngle) * Position(cell_rave_z0, 0 , 0)
    // variables aredefined below
    auto WireAzimuthalAngle = this->dch_data->Get_cell_phi_angle(ilayer, nphi);
    float WireStereoAngle = 0;
    {
      auto l = this->dch_data->database.at(ilayer);
      // radial middle point of the cell at Z=0
      auto cell_rave_z0 = 0.5 * (l.radius_fdw_z0 + l.radius_fuw_z0);
      // when building the twisted tube, the twist angle is defined as:
      //     cell_twistangle    = l.StereoSign() * DCH_i->twist_angle
      // which forces the stereoangle of the wire to have the oposite sign
      WireStereoAngle = (-1.) * l.StereoSign() * dch_data->stereoangle_z0(cell_rave_z0);
    }

//====================================================================================
//                           Drift Time using x-t realtion                          //
//====================================================================================
    // Count total number of electrons from all clusters in this sim hit
    const int total_electrons = std::accumulate(electrons_per_cluster.begin(),
                                            electrons_per_cluster.end(), 0);
    
    // Vectors to store per-electron drift time
    // and reserve memory to avoid the reallocation
    std::vector<double> electron_times_ns;
    electron_times_ns.reserve(total_electrons);	

    // Store per-electron radial distance
    std::vector<double> electron_r_cm;
    electron_r_cm.reserve(total_electrons);

    // Store per-electron z position along the wire for propagation
    std::vector<double> electron_wireZ_mm;
    electron_wireZ_mm.reserve(total_electrons);

    // Build track direction from SimTrackerHit momentum
    TVector3 track_dir(input_sim_hit.getMomentum().x,
                   input_sim_hit.getMomentum().y,
                   input_sim_hit.getMomentum().z);
    if (track_dir.Mag() > 0.0) {
	    track_dir = track_dir.Unit();
    }
    else {
	    track_dir.SetXYZ(0.0, 0.0, 0.0); // safety
    }

    // Start from hit position as reference point
    TVector3 step_start = hit_position;
    const double L_cm = L_mm*MM_TO_CM;

    // In EDM4hep, SimTrackerHit position is usually the midpoint of the step.
    // So the start point is midpoint - (L/2) * direction.
    if (track_dir.Mag() > 0.0 && L_cm > 0.0) {
	    step_start = hit_position - 0.5 * L_cm * track_dir;
    }

    // for 2D LUT, we need x,y
    TVector3 u = wire_direction_ez.Unit();	//unit vector along th Z direction
    TVector3 e1(1.0, 0.0, 0.0);			//unit vector along the global X direction
    e1 -= (e1.Dot(u)) * u;

    if (e1.Mag() < 1e-9) {
	    e1 = TVector3(0.0, 1.0, 0.0);
	    e1 -= (e1.Dot(u)) * u;
    }

    e1 = e1.Unit();
    TVector3 e2 = u.Cross(e1).Unit();		//unit vector along the Y
    //-------------------------------------------------------------------------
    
    double s_cm = 0.0;
    int Ne_cluster = 0;
    double r_cluster_cm = 0.0;
    double t_ns = 0.0;
    double x_cluster_cm = 0.0;
    double y_cluster_cm = 0.0;

    for (size_t icl = 0; icl < cl_positions_mm.size(); ++icl)
    {
	   s_cm = cl_positions_mm[icl]*MM_TO_CM;		//mm -> cm

	   // Read number of electrons in this cluster
	   Ne_cluster = electrons_per_cluster[icl];

	   // Cluster position along the reconstructed track
	   TVector3 cl_pos = step_start + s_cm * track_dir;
	   
	   // Compute vector from cluster position to the closest point on the wire
	   TVector3 cl_to_wire = this->dch_data->Calculate_hitpos_to_wire_vector(
			   hit_ilayer, hit_nphi, cl_pos);

	   // Compute projected cluster position on the wire
	   TVector3 cl_on_wire = cl_pos + cl_to_wire;

	   // Store z coordinate of projected point on the wire in mm
	   const double z_wire_cl_mm = cl_on_wire.Z() / MM_TO_CM;

	   // Perpendicular displacement to the wire
	   TVector3 v_perp = cl_to_wire - (cl_to_wire.Dot(u)) * u;

	   // TRUE 2D coordinates in cm (wire-centered transverse plane)
	   x_cluster_cm = v_perp.Dot(e1);
	   y_cluster_cm = v_perp.Dot(e2);

	   // Keep r only for existing debug histos (optional)
	   r_cluster_cm = v_perp.Mag();

	   for (int ie = 0; ie < Ne_cluster; ++ie) {
		   //drift time of each electron in a cluster
		   t_ns = m_xt2d->sampleTimeNs(x_cluster_cm, y_cluster_cm, myRandom);
		   
		   // Store drift time for the electrons
		   electron_times_ns.push_back(t_ns);

		   // Store electron radial distance
		   electron_r_cm.push_back(r_cluster_cm);

		   // Store electron z position for wire
		   electron_wireZ_mm.push_back(z_wire_cl_mm);

		   if (m_create_debug_histos.value()) {
			   hTd->Fill(t_ns);
		   	   hXT->Fill(r_cluster_cm, t_ns);
		   }
	   }
    } //end of the loop over cluster for drift time calculation:

//==========================================================================================
//                               Apply Polya distribution                                 //
//==========================================================================================
    // Store avalanche charge for each drifted electron
    std::vector<double> electron_charges;
    electron_charges.reserve(electron_times_ns.size());
    double Qtot = 0.0;

    //for each drifted electron, sample the avalunche uisng
    //Polya distribution.
    for (size_t i = 0; i<electron_times_ns.size(); ++i)
    {
	    double q = avalancheCharge(myRandom);

	    electron_charges.push_back(q);
	    Qtot +=q;

	    if (m_create_debug_histos.value()) hAvalancheQ->Fill(q);
    }

//=======================================================================================
// 			Accumulate electrons per cellID				       //
//======================================================================================= 			
    // Create per-cell accumulator for this cellID
    auto& cellAccum = cellAccumMap[static_cast<uint64_t>(cellid)];

    // Accumulate total deposited energy summed over all contributing sim hits in the cell
    cellAccum.eDep += static_cast<float>(input_sim_hit.getEDep());

    // Absolute time of this SimTrackerHit
    const double simTime_ns = input_sim_hit.getTime();

    // Store representative metadata only once for cell	   
    if (!cellAccum.hasMeta) {
	    cellAccum.hasMeta           = true;
	    cellAccum.representativeSimTime_ns = simTime_ns;
	    cellAccum.positionSW        = positionSW;
	    cellAccum.distanceToWire    = distanceToWire;
	    cellAccum.wireAzimuthalAngle = WireAzimuthalAngle;
	    cellAccum.wireStereoAngle   = WireStereoAngle;
	    cellAccum.type              = type;
	    cellAccum.quality           = quality;
	    cellAccum.eDep              = static_cast<float>(input_sim_hit.getEDep());
	    cellAccum.eDepError         = eDepError;
    }
    else {
	    // For now, keep a simple cumulative eDep at cell level
	    cellAccum.eDep += static_cast<float>(input_sim_hit.getEDep());
    }

    // store absolute electron times in the event time frame:
    for (double tdrift_ns : electron_times_ns) {
	    cellAccum.times_ns.push_back(simTime_ns + tdrift_ns);
    }

    cellAccum.charges.insert(cellAccum.charges.end(), 
		    electron_charges.begin(), electron_charges.end());
    cellAccum.wireZ_mm.insert(cellAccum.wireZ_mm.end(), 
		    electron_wireZ_mm.begin(), electron_wireZ_mm.end());
    
    // Store number of electrons per cluster
    for (int ne : electrons_per_cluster) {
	    if (ne < 0) continue;
	    cellAccum.electronsPerCluster.push_back(static_cast<uint16_t>(ne));
    }
    //---------------------------------------------------------------------------------------
  } // end loop over hit collection

//***********************************************************************************************
//***********************************************************************************************
//***********************************************************************************************

  // look for the first cell to save only one debug histogram 
  uint64_t debugCellID = 0;

  // Retrieve accumulated data per cell (all electrons)
  for (const auto& kv : cellAccumMap) 
  {
	  const auto& acc = kv.second;
            
	  if (acc.times_ns.empty()) continue;  
	  if (!acc.hasMeta) continue;

	  // Store the electrons charge, drift time, and 
	  // distance along the wire per cell:
	  const auto& electron_times_ns = acc.times_ns;	   
	  const auto& electron_charges  = acc.charges;	    
	  const auto& electron_wireZ_mm = acc.wireZ_mm;

	    
	  // Local waveform buffers for THIS cell only
	  std::vector<double> wfCenter(nBinsAnalogWf, 0.0);
      	  std::vector<double> wfLeft(nBinsAnalogWf, 0.0);
	  std::vector<double> wfRight(nBinsAnalogWf, 0.0);
	  std::vector<double> wfElecL(nBinsAnalogWf, 0.0);
	  std::vector<double> wfElecR(nBinsAnalogWf, 0.0);


	  const float WireStereoAngle = acc.wireStereoAngle;
	  // Find earliest electron arrival time in cell
	  auto it_min = std::min_element(acc.times_ns.begin(), acc.times_ns.end());
	  // Find earliest electron arrival time in cell
	  const double t_hit_ns = (*it_min < 0.0) ? 0.0 : *it_min;

	  // Select first valid cell to store debug waveform histograms
	  if (m_create_debug_histos.value() && !m_waveformFilled) {
      		  debugCellID = kv.first;
	  }

	  const bool isDebugCell = m_create_debug_histos.value() &&
		  !m_waveformFilled &&
		  (kv.first == debugCellID);
	  //const int nBinsWave = hWaveform ? hWaveform->GetNbinsX() : 0;
	    
//==========================================================================================
//               Build analog waveform by summing the individual Pulses
//==========================================================================================

	  if (!electron_times_ns.empty()) 
	  {
    		  if (isDebugCell && hSignalPulse &&
			  	  electron_times_ns.size() == electron_charges.size()) {

		  	  hSignalPulse->Reset();

			  // Find earliest electron arrival time in this cell
		  	  auto it_min_e = std::min_element(electron_times_ns.begin(), 
					  electron_times_ns.end());

			  // Get index of the earliest electron
		  	  const size_t idx_min = std::distance(electron_times_ns.begin(), 
					  it_min_e);

			  // Use charge of the earliest electron for the reference single-electron pulse
		  	  const double q = electron_charges[idx_min];
			  const double t0_ns = electron_times_ns[idx_min];

			  // Fill the debug pulse histogram bin by bin in time
		  	  for (int ibin = 1; ibin <= hSignalPulse->GetNbinsX(); ++ibin) {
		      		  const double t = hSignalPulse->GetBinCenter(ibin);
		      		  const double v = singleElectronPulse(t, t0_ns, q);
		      		  hSignalPulse->SetBinContent(ibin, v);
			  }
	      	  }

		  // Loop over all time bins of the analog waveform
		  for (int ibin = 0; ibin < nBinsAnalogWf; ++ibin) {
			  const double t = analogTimeAxis[ibin];
			  double V_t = 0.0;

			  // Sum contribution from every electron arriving
		  	  for (size_t i = 0; i < electron_times_ns.size(); ++i) {
		      		  const double t0_ns = electron_times_ns[i];
		      		  const double q     = electron_charges[i];

				  // Add this electron pulse contribution to the waveform sample
		      		  V_t += singleElectronPulse(t, t0_ns, q);
		  	  }
		  	  wfCenter[ibin] = V_t;
	      	  }

//=========================================================================================
//                   propagate waveform to both ends (delay + attenuation)
//=========================================================================================

		  // Stereo angle
	      	  const double cosStereo = std::cos(WireStereoAngle);
	      	  const double cosangle  = (std::abs(cosStereo) > 1e-6) ? std::abs(cosStereo) : 1.0;

		  // Signal propagation speed along the sense wire.
	      	  const double v_mm_per_ns = m_signalSpeed_mm_per_ns.value();
		  
		  // Attenuation length for signal transport along the wire.
	      	  const double lambda_mm   = m_wireAttenuationLength_mm.value();

		  // Sum the propagated signal seen at both wire ends time-bin by time-bin.
	      	  for (int ibin = 0; ibin < nBinsAnalogWf; ++ibin) {
		  	  const double t = analogTimeAxis[ibin];

		  	  double VL = 0.0;
		  	  double VR = 0.0;

			  // Each electron reaches the two ends with different delay and attenuation.
		  	  for (size_t i = 0; i < electron_times_ns.size(); ++i) {
				  // Electron production point along the wire.
		      		  const double z_wire_mm = electron_wireZ_mm[i];

				  // Propagation distance from the electron point to the left-right wire end.
		      		  double distR_mm = (m_halfChamberLength_mm - z_wire_mm) / cosangle;
		      		  double distL_mm = (m_halfChamberLength_mm + z_wire_mm) / cosangle;
		      		  if (distR_mm < 0.0) distR_mm = 0.0;
		      		  if (distL_mm < 0.0) distL_mm = 0.0;

				  // Convert wire-travel distance into propagation time.
		      		  const double tpropR_ns = (v_mm_per_ns > 0.0) ? distR_mm / v_mm_per_ns : 0.0;
		      		  const double tpropL_ns = (v_mm_per_ns > 0.0) ? distL_mm / v_mm_per_ns : 0.0;
				  
				  // Left, Right-end amplitude loss from wire attenuation.
				  const double attR = (m_enableWireAttenuation.value() && lambda_mm > 0.0)
					  ? std::exp(-distR_mm / lambda_mm)
					  : 1.0;
				  const double attL = (m_enableWireAttenuation.value() && lambda_mm > 0.0)
					  ? std::exp(-distL_mm / lambda_mm)
					  : 1.0;

				  // Local electron time inside the waveform window.
				  const double t0_ns = electron_times_ns[i];

				  // Left, Right-end signal = delayed and attenuated single-electron pulse.
				  VR += singleElectronPulse(t, t0_ns + tpropR_ns,
						  electron_charges[i] * attR);

				  VL += singleElectronPulse(t, t0_ns + tpropL_ns,
						  electron_charges[i] * attL);
			  }

			  // Store propagated waveform amplitudes at the two readout ends.
			  wfLeft[ibin]  = VL;
			  wfRight[ibin] = VR;
		  }

	      	  // Fill debug histograms only for one example cell
	      	  if (isDebugCell) {
		  	  hWaveform->Reset();
		  	  hWaveformL->Reset();
		  	  hWaveformR->Reset();
		  	  hWaveformElecL->Reset();
		  	  hWaveformElecR->Reset();

			  // Fill central, left-end and right-end waveform of total induced signal.
		  	  for (int ibin = 1; ibin <= nBinsAnalogWf; ++ibin) {
		      		  hWaveform->SetBinContent(ibin,  wfCenter[ibin - 1]);
		      		  hWaveformL->SetBinContent(ibin, wfLeft[ibin - 1]);
		      		  hWaveformR->SetBinContent(ibin, wfRight[ibin - 1]);
		  	  }
	      	  }
	  }
// ==========================================================================================//
//                          impedance mismatch (reflection/transmission)                     //
//===========================================================================================//
//This function return the electronics transfer function impulse shape
	//by using the rise and fall time
	
	  auto gmct_signalShape2 = [&](double t0_ns, double amp) -> double 
	  {
		const double tauFall1 = m_elecTauFall1_ns.value();
		const double tauRise  = m_elecTauRise_ns.value();
		const double mix      = m_elecMixFraction.value();
		const double tauFall2 = m_elecTauFall2_ns.value();

		if (t0_ns <= 0.0) return 0.0;
		if (tauFall1 <= 0.0 || tauFall2 <= 0.0 || tauRise <= 0.0) return 0.0;

		// Double-exponential decay describing front-end shaping tail.
		double sign = mix * std::exp(-t0_ns / tauFall1) / tauFall1;
		sign       += (1.0 - mix) * std::exp(-t0_ns / tauFall2) / tauFall2;
		// Rising edge shaping
		sign       *= (1.0 - std::exp(-t0_ns / tauRise));

		 // Normalization factor
		 sign       *= (tauFall1 + tauRise) * (tauRise + tauFall2);
		 sign       /= (tauFall1 * tauFall2 + tauRise * tauFall2
				 + tauFall1 * tauRise * mix - tauRise * mix * tauFall2);

		 // Return shaped signal scaled by input amplitude.
		 return sign * amp;
	  };
	  //------------------------------------------------------------------------------------//

	  // Apply impedance mismatch + electronics shaping via convolution.
	  auto apply_step10_electronics_vec = [&](const std::vector<double>& in, 
			  std::vector<double>& out)
		  {
			  // Initialize output waveform
			  out.assign(in.size(), 0.0);
			  const int nBins = static_cast<int>(in.size());
			  if (nBins <= 0) return;

			  const double dt_ns = dtAnalog_ns;
			  
			  // Start with front-end gain scaling.
			  double scale = m_frontEndGain.value();

			  // Impedance mismatch modifies transmitted signal amplitude
			  // (reflection losses).
			  if (m_enableImpedanceMismatch.value()) {
				  const double R = m_matchingRes_Ohm.value();
				  const double Z = m_tubeImpedance_Ohm.value();
				  const double denom = (R + Z);

				  // Reflection coefficient (R-Z)/(R+Z)
				  if (std::abs(denom) > 0.0) {
					  const double reflect = (R - Z) / denom;
					  scale *= (1.0 - reflect);
				  }
			  }

			  // Apply global gain + transmission scaling to waveform.
			  std::vector<double> x(nBins, 0.0);
			  for (int i = 0; i < nBins; ++i) {
				  x[i] = in[i] * scale;
			  }
			  // If TF disabled, just copy scaled waveform
			  if (!m_enableElectronicsTF.value()) {
				  out = x;
				  return;
			  }

			  // Build finite-length impulse response kernel of electronics.
			  const int nKernel = std::min( nBins, std::max(1, 
						  int(std::ceil(m_elecKernelLength_ns.value() / dt_ns)) + 1));
			  
			  // Discrete electronics response sampled in time.
			  std::vector<double> k(nKernel, 0.0);    
			  for (int ik = 0; ik < nKernel; ++ik) {
				  k[ik] = gmct_signalShape2(ik * dt_ns, 1.0);
			  }

			  // Convolve input waveform with electronics response (signal shaping).
			  for (int i = 0; i < nBins; ++i) {
				  double accConv = 0.0;
				  const int jmax = std::min(i, nKernel - 1);
				  for (int j = 0; j <= jmax; ++j) {
					  accConv += x[i - j] * k[j];
				  }

				  // Multiply by dt to approximate continuous convolution integral.
				  out[i] = accConv * dt_ns;
		      	  }
		  };


	  // Left-end and right-end waveform after transmission losses
	  // and electronics shaping.
	  apply_step10_electronics_vec(wfLeft,  wfElecL);
	  apply_step10_electronics_vec(wfRight, wfElecR);
//============================================================================================//
//                                add colored noise using FFT/IFFT                            //
//============================================================================================//

	  // Add electronic noise to the shaped waveform using an FFT noise template.
	  if (m_enableFFTNoise.value() && m_fftNoiseReady)
	  {
		  // Use the same samples size as the analog waveform for the noise.
	      	  const int nBins = nBinsAnalogWf;
		  // Analog sampling step sets the waveform Nyquist frequency for the noise generation.
	      	  const double dt_ns = dtAnalog_ns;

		  // Generate one colored-noise for the left channel
		  // and one for the right channel
	      	  auto noiseL = dchfft::makeFFTNoiseSamples(
			  	  nBins, dt_ns, m_fftMag, m_fftMaxFreq, myRandom, m_noiseRemoveDC.value());
	      	  auto noiseR = dchfft::makeFFTNoiseSamples(
			  	  nBins, dt_ns, m_fftMag, m_fftMaxFreq, myRandom, m_noiseRemoveDC.value());

		  // Rescale the generated noise to the requested detector/electronics noise amplitude.
	      	  const double scale =
		  	  (m_fftAmpNorm != 0.0) ? (m_noiseScale.value() / m_fftAmpNorm) : 0.0;

		  // Add the scaled colored noise sample-by-sample to the electronics-shaped waveform.
	      	  for (int i = 0; i < nBins; ++i) {
		  	  wfElecL[i] += scale * noiseL[i];
		  	  wfElecR[i] += scale * noiseR[i];
	      	  }
		  // Fill ONE noise histogram for debugging (only first cell/event)
		  if (m_create_debug_histos.value() && !m_waveformFilled) {
			  hNoise->Reset();
			  for (int ibin = 1; ibin <= nBins; ++ibin) {
				  hNoise->SetBinContent(ibin, scale * noiseL[ibin - 1]);
			  }
		  }
	  }
//================================================================================================//
//                     Analog(L/R) -> Digitized(L/R) using the SAME time bin:                     //
//================================================================================================//

	  // Convert analog waveforms to discrete ADC samples
	  std::vector<uint16_t> adcL, adcR;
	  std::vector<float> digiL_mV, digiR_mV;

	  if (m_doWaveformDigitization.value())
	  {
		  // Maximum ADC code defined by the number of bits
		  const int adcMax = (m_adcBits.value() >= 1) ? ((1 << m_adcBits.value()) - 1) : 0;
		  // Voltage step per ADC count
		  const double lsb_mV = m_adcLSB_mV.value();
		  // Baseline offset representing electronics pedestal.
		  const double base_mV = m_adcBaseline_mV.value();
		  // Electronics polarity defines signal sign convention
		  const int pol = (m_adcPolarity.value() >= 0) ? +1 : -1;

		  // ADC sampling period defines discrete time grid of digitization.
		  const double dtDigi_ns = m_adcSamplePeriod_ns.value();
		  // Number of digitized samples within the waveform time window.
		  const int nBinsDigi = std::max(
				  1, static_cast<int>(std::lround((wfTmax - wfTmin) / dtDigi_ns)));

		  // Reserve memory for ADC codes and reconstructed voltages.
	      	  adcL.reserve(nBinsDigi);
	      	  adcR.reserve(nBinsDigi);
	      	  digiL_mV.reserve(nBinsDigi);
	      	  digiR_mV.reserve(nBinsDigi);

		  // Loop over ADC sampling times.
	      	  for (int ibin = 0; ibin < nBinsDigi; ++ibin) {
		  	  const double t_ns = wfTmin + (ibin + 0.5) * dtDigi_ns;

		  	  int idxA = static_cast<int>(std::floor((t_ns - wfTmin) / dtAnalog_ns));
		  	  if (idxA < 0) idxA = 0;
		  	  if (idxA >= nBinsAnalogWf) idxA = nBinsAnalogWf - 1;

			  // Analog voltage at this sampling point
		  	  const double analogL_V = wfElecL[idxA];
		  	  const double analogR_V = wfElecR[idxA];

			  // Convert analog voltage to mV, apply polarity and add baseline offset.
		  	  const double vinL_mV = base_mV + pol * (1000.0 * analogL_V);
		  	  const double vinR_mV = base_mV + pol * (1000.0 * analogR_V);

			  // Quantize input voltage into ADC integer code.
		  	  int codeL = static_cast<int>(vinL_mV / lsb_mV);
		  	  int codeR = static_cast<int>(vinR_mV / lsb_mV);

			  // Saturate ADC codes within allowed dynamic range.
		  	  if (codeL < 0) codeL = 0;
		  	  if (codeL > adcMax) codeL = adcMax;
		  	  if (codeR < 0) codeR = 0;
		  	  if (codeR > adcMax) codeR = adcMax;

			  // Store digitized ADC counts.
		  	  adcL.push_back(static_cast<uint16_t>(codeL));
		  	  adcR.push_back(static_cast<uint16_t>(codeR));
		  	  digiL_mV.push_back(static_cast<float>(codeL * lsb_mV));
		  	  digiR_mV.push_back(static_cast<float>(codeR * lsb_mV));
		  }
	  }
	  //---------------------------------------------------------------------------------//

	  // Store one representative waveform per event
	  // to monitor the full signal chain without histogram overlap.
	  if (isDebugCell)

	  {
		  if (hWaveform) {
			  hWaveform->Reset();
			  for (int ibin = 1; ibin <= nBinsAnalogWf; ++ibin) {
				  hWaveform->SetBinContent(ibin, wfCenter[ibin - 1]);
			  }
		  }
		  if (hWaveformL) {
			  hWaveformL->Reset();
			  for (int ibin = 1; ibin <= nBinsAnalogWf; ++ibin) {
				  hWaveformL->SetBinContent(ibin, wfLeft[ibin - 1]);
			  }
		  }
		  if (hWaveformR) {
			  hWaveformR->Reset();
			  for (int ibin = 1; ibin <= nBinsAnalogWf; ++ibin) {
				  hWaveformR->SetBinContent(ibin, wfRight[ibin - 1]);
			  }
		  }
		  if (hWaveformElecL) {
			  hWaveformElecL->Reset();
			  for (int ibin = 1; ibin <= nBinsAnalogWf; ++ibin) {
				  hWaveformElecL->SetBinContent(ibin, wfElecL[ibin - 1]);
			  }
		  }
		  if (hWaveformElecR) {
			  hWaveformElecR->Reset();
			  for (int ibin = 1; ibin <= nBinsAnalogWf; ++ibin) {
				  hWaveformElecR->SetBinContent(ibin, wfElecR[ibin - 1]);
			  }
		  }
		  if (hWaveformDigiL && !digiL_mV.empty() &&
				  (int)digiL_mV.size() == hWaveformDigiL->GetNbinsX()) {
			  hWaveformDigiL->Reset();
			  for (int ibin = 1; ibin <= hWaveformDigiL->GetNbinsX(); ++ibin) {
				  hWaveformDigiL->SetBinContent(ibin, digiL_mV[ibin - 1]);
			  }
		  }
		  if (hWaveformDigiR && !digiR_mV.empty() &&
				  (int)digiR_mV.size() == hWaveformDigiR->GetNbinsX()) {
			  hWaveformDigiR->Reset();
			  for (int ibin = 1; ibin <= hWaveformDigiR->GetNbinsX(); ++ibin) {
				  hWaveformDigiR->SetBinContent(ibin, digiR_mV[ibin - 1]);
			  }
		  }
		  m_waveformFilled = true;
	  } // End of debugging	  
	  //--------------------------------------------------------------------------------------//

            //extension::MutableSenseWireHit oDCHdigihit;
            auto  oDCHdigihit = output_digi_hits.create();

            oDCHdigihit.setCellID(kv.first);
            oDCHdigihit.setType(acc.type);
            oDCHdigihit.setQuality(acc.quality);

            oDCHdigihit.setTime(t_hit_ns);
            oDCHdigihit.setEDep(acc.eDep);
            oDCHdigihit.setEDepError(acc.eDepError);
            oDCHdigihit.setPosition(acc.positionSW);

            oDCHdigihit.setPositionAlongWireError(m_z_resolution);
            oDCHdigihit.setWireAzimuthalAngle(acc.wireAzimuthalAngle);
            oDCHdigihit.setWireStereoAngle(acc.wireStereoAngle);
            oDCHdigihit.setDistanceToWire(acc.distanceToWire);
            oDCHdigihit.setDistanceToWireError(m_xy_resolution);

	    for (auto ne : acc.electronsPerCluster) {
		    oDCHdigihit.addToNElectrons(ne);
	    }


	    for (const auto& simHit : acc.simHits) {
                    auto link = output_digi_sim_association.create();
                    link.setFrom(oDCHdigihit);
                    link.setTo(simHit);
            }

    	const float t0_ns = static_cast<float>(t_hit_ns);

    	// left side digi waveform:
    	if (m_doWaveformDigitization.value() && !adcL.empty()) {

		auto DigiwfL = out_Digitized_waveformsL.create();
        	DigiwfL.setCellID(kv.first);
        	//DigiwfL.setSide(1);                // convention: 1 = Left
        	DigiwfL.setTime(t0_ns);
        	DigiwfL.setInterval(m_adcSamplePeriod_ns.value());
		DigiwfL.setCharge(m_adcLSB_mV.value());

		for (uint16_t code : adcL) {
			DigiwfL.addToAdcCounts(code);
		}

		//DigiwfL.setHit(oDCHdigihit);
    	}

    	// left side digi waveform:
    	if (m_doWaveformDigitization.value() && !adcR.empty()) {

        	auto DigiwfR = out_Digitized_waveformsR.create();
        	DigiwfR.setCellID(kv.first);
        	//DigiwfR.setSide(0);                // convention: 0 = Right
        	DigiwfR.setTime(t0_ns);
        	DigiwfR.setInterval(m_adcSamplePeriod_ns.value());
		DigiwfR.setCharge(m_adcLSB_mV.value());

        	for (uint16_t code : adcR) {
			DigiwfR.addToAdcCounts(code);
		}

    	}
    
  } // End of the Cell loop:
//===============================================================================================//
//                       END OF L0 DIGITIZED WAVEFORM FOR IDEA DCH                               //
//===============================================================================================//

  m_event_counter+=1;

/////////////////////////////////////////////////////////////////
// return the digitized hits, truth links, and left/right waveform collections for this event
  return std::make_tuple<extension::SenseWireHitCollection, 
	 extension::SenseWireHitSimTrackerHitLinkCollection,
	 
	 edm4hep::RawTimeSeriesCollection,
	 edm4hep::RawTimeSeriesCollection
		 >(
      std::move(output_digi_hits), 
      std::move(output_digi_sim_association),

      std::move(out_Digitized_waveformsL),
      std::move(out_Digitized_waveformsR)
      );

} // End of the operator Block

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       finalize       //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

StatusCode DCHdigi_v03::finalize() 
{

  if (m_polya) {
	  delete m_polya;
	  m_polya = nullptr;
  }


  return StatusCode::SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       ThrowException       ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void DCHdigi_v03::ThrowException(std::string s) const 
{
  error() << s.c_str() << endmsg;
  throw std::runtime_error(s);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       PrintConfiguration       ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void DCHdigi_v03::PrintConfiguration(std::ostream& io) 
{
  io << "DCHdigi will use the following components:\n";
  io << "\tGeometry Service: " << m_geoSvcName.value().c_str() << "\n";
  io << "\tUID Service: " << m_uidSvcName.value().c_str() << "\n";
  io << "\tDetector name: " << m_DCH_name.value().c_str() << "\n";
  io << "\t\t|--Volume bitfield: " << m_decoder->fieldDescription().c_str() << "\n";
  io << "\t\t|--Number of layers: " << dch_data->database.size() << "\n";
  io << "\tResolution along the wire (mm): " << m_z_resolution.value() << "\n";
  io << "\tResolution perp. to the wire (mm): " << m_xy_resolution.value() << "\n";
    io << "\tCreate debug histograms (THistSvc): " << (m_create_debug_histos.value() ? "true" : "false") << "\n";
  return;
}

bool DCHdigi_v03::IsParticleCreatedInsideDriftChamber(const edm4hep::MCParticle& thisParticle) const {
  auto vertex = thisParticle.getVertex(); // in mm
  auto vertexRsquared = vertex[0] * vertex[0] + vertex[1] * vertex[1];
  auto vertexZabs = std::fabs(vertex[2]);
  float DCH_halflengh = dch_data->Lhalf / dd4hep::mm;                // 2000;
  float DCH_rin_squared = std::pow(dch_data->rin / dd4hep::mm, 2);   // 350 * 350;
  float DCH_rout_squared = std::pow(dch_data->rout / dd4hep::mm, 2); // 2000 * 2000;
  return (vertexZabs < DCH_halflengh) && (vertexRsquared > DCH_rin_squared) && (vertexRsquared < DCH_rout_squared);
}


