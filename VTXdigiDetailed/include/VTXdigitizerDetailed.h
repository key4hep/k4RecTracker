#pragma once

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
#include "edm4hep/TrackCollection.h"
#if __has_include("edm4hep/TrackerHitPlaneCollection.h")
#include "edm4hep/TrackerHitPlaneCollection.h"
#else
#include "edm4hep/TrackerHitCollection.h"
namespace edm4hep {
  using TrackerHitPlaneCollection = edm4hep::TrackerHitCollection;
}  // namespace edm4hep
#endif
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

// DD4HEP
#include "DD4hep/Detector.h"  // for dd4hep::VolumeManager
#include "DD4hep/DetType.h"
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

private:
  // Input sim vertex hit collection name
  mutable DataHandle<edm4hep::SimTrackerHitCollection> m_input_sim_hits{"inputSimHits", Gaudi::DataHandle::Reader, this};
  // Output digitized vertex hit collection name
  mutable DataHandle<edm4hep::TrackerHitPlaneCollection> m_output_digi_hits{"outputDigiHits", Gaudi::DataHandle::Writer, this};
  // Output link between sim hits and digitized hits
  mutable DataHandle<edm4hep::TrackerHitSimTrackerHitLinkCollection> m_output_sim_digi_link{"outputSimDigiAssociation", Gaudi::DataHandle::Writer, this};

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

  // t resolution in ns
  Gaudi::Property<std::vector<float>> m_t_resolution{this, "tResolution", {0.1}, "Time resolutions per layer [ns]"};

  // Surface manager used to project hits onto sensitive surface with forceHitsOntoSurface argument
  mutable const dd4hep::rec::SurfaceMap* _map;

  // Option to force hits onto sensitive surface
  BooleanProperty m_forceHitsOntoSurface{this, "forceHitsOntoSurface", false, "Project hits onto the surface in case they are not yet on the surface (default: false"};

  // Tangent of sensor's Lorentz angle (default is 0.1)
  FloatProperty m_tanLorentzAnglePerTesla{this, "tanLorentzAnglePerTesla", {0.1}, "Tangent of sensor's Lorentz angle per Tesla (default is 0.1)"};

  FloatProperty m_Sigma50{this, "Sigma50", 0.00151, "Charge diffusion in millimeters for 50 micron Si. Default = 0.00151mm taken from CMS"}; 

  float m_Dist50; //=0.050  // Define 50microns in mm

  float m_ClusterWidth; //=3.0  // Define the area of integration of the charge in number of sigma around the central position
  
  // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;

  // Gaussian random number generator used for time smearing
  std::vector<Rndm::Numbers> m_gauss_t_vec;

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

  typedef std::map<int,std::map<int,float,std::less<int>>,std::less<int>> hit_map_type;
    
private:
  // Additional member functions
  // Private methods
  template<typename T> void getSensorThickness();
  
  void primary_ionization(const edm4hep::SimTrackerHit& hit,
			  std::vector<ChargeDepositUnit>& ionizationPoints) const;

  void drift(const edm4hep::SimTrackerHit& hit,
	     const std::vector<ChargeDepositUnit>& ionizationPoints,
	     std::vector<SignalPoint>& collectionPoints) const;

  dd4hep::rec::Vector3D DriftDirection(const edm4hep::SimTrackerHit& hit) const;
  
  void get_charge_per_pixel(const edm4hep::SimTrackerHit& hit,
			    const std::vector<SignalPoint>& collectionPoints,
			    hit_map_type& hit_map) const;

  void generate_output(const edm4hep::SimTrackerHit hit, edm4hep::TrackerHitPlaneCollection* output_digi_hits, edm4hep::TrackerHitSimTrackerHitLinkCollection* output_sim_digi_link_col, const hit_map_type& hit_map) const;

  void SetProperDirectFrame(TGeoHMatrix& sensorTransformMatrix) const;
  
};
