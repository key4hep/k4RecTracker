#pragma once

// ROOT headers
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDirectory.h"

// GAUDI
#include "Gaudi/Property.h"
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

// K4FWCORE
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"

// EDM4HEP
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

// DD4HEP
#include "DD4hep/Detector.h"  // for dd4hep::VolumeManager
#include "DD4hep/DetType.h"
#include "DD4hep/Shapes.h"
#include "DDRec/DetectorData.h" // for detector extensions
#include "DDRec/Vector3D.h"
#include "DDRec/Vector2D.h"
#include "DDRec/SurfaceManager.h"

#include "DDSegmentation/BitFieldCoder.h"

#include <vector>

/** @class VTXdigitizerDetailed
 *
 *  Algorithm for creating digitized (meaning 'reconstructed' for now) vertex detector hits (edm4hep::TrackerHitPlane) from Geant4 hits (edm4hep::SimTrackerHit).
 *  It considers the production, drift and collection of charges in pixel/strip sensors. The time is considered with a simple smear for now 
 *
 *  @author Jessy Daniel, Brieuc Francois, Armin Ilg
 *  @date   2025-02
 *
 */

class VTXdigitizerDetailed : public Gaudi::Algorithm {
public:
  explicit VTXdigitizerDetailed(const std::string&, ISvcLocator*);
  virtual ~VTXdigitizerDetailed();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Execute.
   *   @return status code
   */
  virtual StatusCode execute(const EventContext&) const final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

  // Typedef public si besoin
  typedef std::map<int, std::map<int, float>> hit_map_type;
  
private:
  // Input sim vertex hit collection name
  mutable DataHandle<edm4hep::SimTrackerHitCollection> m_input_sim_hits{"inputSimHits", Gaudi::DataHandle::Reader, this};
  // Output digitized vertex hit collection name
  mutable DataHandle<edm4hep::TrackerHitPlaneCollection> m_output_digi_hits{"outputDigiHits", Gaudi::DataHandle::Writer, this};
  // Output link between sim hits and digitized hits
  mutable DataHandle<edm4hep::TrackerHitSimTrackerHitLinkCollection> m_output_sim_digi_link{"outputSimDigiAssociation", Gaudi::DataHandle::Writer, this};

  // Threshold 
  bool Apply_Threshold(double& ChargeInPixel) const;
  // Detector name
  Gaudi::Property<std::string> m_detectorName{this, "detectorName", "Vertex", "Name of the detector (default: Vertex)"};
  // Detector readout names
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "VertexBarrelCollection", "Name of the detector readout"};
  // Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  // Decoder for the cellID
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  // Volume manager to get the physical cell sensitive volume
  dd4hep::VolumeManager m_volman;

  // Normal Vector direction in sensor local frame (may differ according to geometry definition within k4geo)
  StringProperty m_LocalNormalVectorDir{this, "LocalNormalVectorDir", "", "Normal Vector direction in sensor local frame (may differ according to geometry definition within k4geo)"};

  // List of sensor thickness per layer in millimeter
  std::vector<float> m_sensorThickness;

  // 2D dim
  std::vector<float> m_sensorWidth;
  std::vector<float> m_sensorLength;

  // Pour ZDiskPetalsData (endcap)
  std::vector<float> m_sensorWidthInner;
  std::vector<float> m_sensorWidthOuter;

  // t resolution in ns
  Gaudi::Property<std::vector<float>> m_t_resolution{this, "tResolution", {0.1}, "Time resolutions per layer [ns]"};

  // Threshold in electron 

  Gaudi::Property<float> m_Threshold{this, "Threshold", 0.0 , "Charge Threshold in e (default: 0.0)"};

  //Threshold smearing in electron

  Gaudi::Property<float> m_ThresholdSmearing{this, "ThresholdSmearing", 0.0, "Sigma of Threshold Gaussian Smearing in e (default: 0.0)"};

  // Surface manager used to project hits onto sensitive surface with forceHitsOntoSurface argument
  mutable const dd4hep::rec::SurfaceMap* _map;

  // Option to force hits onto sensitive surface
  BooleanProperty m_forceHitsOntoSurface{this, "forceHitsOntoSurface", false, "Project hits onto the surface in case they are not yet on the surface (default: false"};

  // Decephiring CellID / Surface Not Found for CellID Issue
  Gaudi::Property<size_t> m_cellIDBits{this, "CellIDBits", 64, "Number of bits to use for the cellID of the hits"}; 
  // Mask to use for the cellID of the hits
  std::uint64_t m_mask{static_cast<std::uint64_t>(-1)}; 
  // Map to store the surface for each cellID
  const dd4hep::rec::SurfaceMap* surfaceMap = nullptr;


  // Tangent of sensor's Lorentz angle (default is 0.1)
  FloatProperty m_tanLorentzAnglePerTesla{this, "tanLorentzAnglePerTesla", {0.1}, "Tangent of sensor's Lorentz angle per Tesla (default is 0.1)"};

  FloatProperty m_Sigma50{this, "Sigma50", 0.00151, "Charge diffusion in millimeters for 50 micron Si. Default = 0.00151mm taken from CMS"}; 

  float m_Dist50; //=0.050  // Define 50microns in mm

  float m_ClusterWidth; //=3.0  // Define the area of integration of the charge in number of sigma around the central position
  
  // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;

  // Gaussian random number generator used for time smearing
  std::vector<Rndm::Numbers> m_gauss_t_vec;

  Rndm::Numbers m_gauss_threshold; // pour faire le smearing du threshold

  // Define a class for 3D ionization points and energy
  /**
   * Internal use only.
   */
  class ChargeDepositUnit {
  public:
    ChargeDepositUnit() : _charge(0), _position(0, 0, 0) {}
    ChargeDepositUnit(float charge, float x, float y, float z) : _charge(charge), _position(x, y, z) {}
    ChargeDepositUnit(float charge, dd4hep::rec::Vector3D position) : _charge(charge), _position(position) {}
    float x() const { return _position.x(); }
    float y() const { return _position.y(); }
    float z() const { return _position.z(); }
    float charge() const {return _charge; }

  private:
    float _charge;
    dd4hep::rec::Vector3D _position;
  }; // End ChargeDepositUnit class definition

  // Define class to store signals on collecion surface
  /**
   * Internal useonly.
   */
  class SignalPoint {
  public:
    SignalPoint() : _pos(0, 0), _time(0), _amplitude(0), _sigma_x(1.), _sigma_y(1.), _hitp(nullptr) {}

    SignalPoint(float x, float y, float sigma_x, float sigma_y, float t, float a = 1.0)
      : _pos(x, y), _time(t), _amplitude(a), _sigma_x(sigma_x), _sigma_y(sigma_y), _hitp(nullptr) {}

    SignalPoint(float x, float y, float sigma_x, float sigma_y, float t, const edm4hep::SimTrackerHit& hit, float a = 1.0)
      : _pos(x, y), _time(t), _amplitude(a), _sigma_x(sigma_x), _sigma_y(sigma_y), _hitp(&hit) {}

    const dd4hep::rec::Vector2D position() const { return _pos; }
    float x() const { return _pos.u(); }
    float y() const { return _pos.v(); }
    float sigma_x() const { return _sigma_x; }
    float sigma_y() const { return _sigma_y; } 
    float time() const { return _time; }
    float amplitude() const { return _amplitude; }
    const edm4hep::SimTrackerHit& hit() { return *_hitp; }
    SignalPoint& set_amplitude(float amp) {
      _amplitude = amp;
      return *this;
    }
    
  private:
    dd4hep::rec::Vector2D _pos;
    float _time;
    float _amplitude;
    float _sigma_x;
    float _sigma_y;
    const edm4hep::SimTrackerHit* _hitp;
  }; // End SignalPoint class definition

  //typedef std::map<int,std::map<int,float,std::less<int>>,std::less<int>> hit_map_type;
    
private:
  // Additional member functions
  // Private methods
  template<typename T> void getSensorThickness();
  template<typename T> void getSensorSize();
  
  void primary_ionization(const edm4hep::SimTrackerHit& hit,
			  std::vector<ChargeDepositUnit>& ionizationPoints) const;

  void drift(const edm4hep::SimTrackerHit& hit,
	     const std::vector<ChargeDepositUnit>& ionizationPoints,
	     std::vector<SignalPoint>& collectionPoints) const;

  dd4hep::rec::Vector3D DriftDirection(const edm4hep::SimTrackerHit& hit) const;
  
  void get_charge_per_pixel(const edm4hep::SimTrackerHit& hit,
			    const std::vector<SignalPoint>& collectionPoints,
			    hit_map_type& hit_map) const;

  void clampCloudToSensorBounds(
     float& CloudMinX, float& CloudMaxX,
     float& CloudMinY, float& CloudMaxY,
     float CloudCenterX, float CloudCenterY,
     const edm4hep::SimTrackerHit& hit) const;

  void generate_output(const edm4hep::SimTrackerHit hit, edm4hep::TrackerHitPlaneCollection* output_digi_hits, edm4hep::TrackerHitSimTrackerHitLinkCollection* output_sim_digi_link_col, const hit_map_type& hit_map) const;

  void SetProperDirectFrame(TGeoHMatrix& sensorTransformMatrix) const;

private:
  // Additional Debug information

  BooleanProperty m_DebugHistos{this, "DebugHistos", false, "If set to true, create a file containing debug histograms. Should then also set DebugFileName to set the name of the debug root file."};
  Gaudi::Property<std::string> m_DebugFileName{this,"DebugFileName", "", "Name of the file containing debug histograms."};

  TH1D* hErrorX; // Histogram to store the distance in X between the true hit and digitized one in mm
  TH1D* hErrorY; // Histogram to store the distance in Y between the true hit and digitized one in mm
  TH1D* hErrorZ; // Histogram to store the distance in Z between the true hit and digitized one in mm
  TH1D* hError;  // Histogram to store the distance between the true hit and digitized one in mm

  TH1D* hChargeAboveThreshold; // Histogram to store the pixel charge after Threshold 
  TH1D* hChargeBeforeThreshold; // Histogram to store the pixel charge before Threshold
  TH1D* hChargePerClusterOrDigis; // Histogram to store charge per Digis ie Cluster 

  TH1D* hActivePixelCountBeforeThreshold;  // Histogram to store the number of active pixels per Cluster before Threshold
  TH1D* hActivePixelCountAfterThreshold; // Histogram to store the number of active pixels per Cluster after Threshold
  

  TH1D* hXDriftDueToMagField; // Histogram to store the X drift due to magnetic field in mm 
  TH1D* hYDriftDueToMagField; // Histogram to store the Y drift due to magnetic field in mm
  
  TH1D* hDigisPerLayer; // Histo to get hits per layer for occupancy studies via cellID
  TH1D* hActivePixelPerlayer; // histo occupancy
  

  void Create_outputROOTfile_for_debugHistograms() const;
};
