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
#include "DDRec/Vector3D.h"
#include "DDRec/Vector2D.h"
#include "DDRec/SurfaceManager.h"

#include "DDSegmentation/BitFieldCoder.h"

#include <vector>

/** @class VTXdigitizerDetailed
 *
 *  Algorithm for creating digitized (meaning 'reconstructed' for now) vertex detector hits (edm4hep::TrackerHitPlane) from Geant4 hits (edm4hep::SimTrackerHit).
 *  
 *  @author Brieuc Francois, Armin Ilg, Jessy Daniel
 *  @date   2025-01
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

  // x resolution in um
  Gaudi::Property<std::vector<float>> m_x_resolution{this, "xResolution", {0.1}, "Spatial resolutions in the x direction per layer [um] (r-phi direction in barrel, z direction in disks)"};

  // y resolution in um
  Gaudi::Property<std::vector<float>> m_y_resolution{this, "yResolution", {0.1}, "Spatial resolutions in the y direction per layer [um] (r direction in barrel, r-phi direction in disks)"};

  // t resolution in ns
  Gaudi::Property<std::vector<float>> m_t_resolution{this, "tResolution", {0.1}, "Time resolutions per layer [ns]"};

  // Surface manager used to project hits onto sensitive surface with forceHitsOntoSurface argument
  mutable const dd4hep::rec::SurfaceMap* _map;

  // Option to force hits onto sensitive surface
  BooleanProperty m_forceHitsOntoSurface{this, "forceHitsOntoSurface", false, "Project hits onto the surface in case they are not yet on the surface (default: false"};

  // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;

  // Gaussian random number generator used for smearing
  std::vector<Rndm::Numbers> m_gauss_x_vec;
  std::vector<Rndm::Numbers> m_gauss_y_vec;
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
    
private:
  // Additioal member functions
  // Private methods
  void primary_ionization(const edm4hep::SimTrackerHit& hit,
			  std::vector<ChargeDepositUnit>& ionizationPoints) const;
  void drift(const edm4hep::SimTrackerHit& hit,
	     const std::vector<ChargeDepositUnit>& ionizationPoints,
	     std::vector<SignalPoint>& collectionPoints) const;

  void induce_signal(const edm4hep::SimTrackerHit& hit,
		     const std::vector<SignalPoint>& collectionPoints,
		     std::vector<dd4hep::DDSegmentation::CellID,float>& signal_map) const;
  
};
