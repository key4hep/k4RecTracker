#pragma once

// ROOT headers
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <memory> // for std::unique_ptr
#include <optional>

// GAUDI
#include "Gaudi/Property.h"

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IMetaDataSvc.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/SmartIF.h"
#include "k4FWCore/MetadataUtils.h"
#include <k4FWCore/Transformer.h>

// K4FWCORE
#include "k4Interface/IGeoSvc.h"

// EDM4HEP
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

// DD4HEP
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"    // for dd4hep::VolumeManager
#include "DDRec/DetectorData.h" // for detector extensions
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Vector2D.h"
#include "DDRec/Vector3D.h"

#include "DDSegmentation/BitFieldCoder.h"

#include <vector>

/** @class VTXdigitizerDetailed
 *
 *  Algorithm for creating digitized (meaning 'reconstructed' for now) vertex detector hits (edm4hep::TrackerHitPlane)
 * from Geant4 hits (edm4hep::SimTrackerHit). It considers the production, drift and collection of charges in
 * pixel/strip sensors. The time is considered with a simple smear for now
 *
 *  @author Jessy Daniel, Adrien Sabard
 *  @date   2025-09
 *
 */

class VTXdigitizerDetailed final // Gaudi::Functional::
    : public k4FWCore::MultiTransformer<
          std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection>(
              const edm4hep::SimTrackerHitCollection&)> {
public:
  explicit VTXdigitizerDetailed(const std::string& aname, ISvcLocator* asvcLoc);
  ~VTXdigitizerDetailed() override = default;

  /**  Initialize.
   *   @return status code
   */
  StatusCode initialize() override;

  /**  Operator (Execute).
   *   @return status code
   */
  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection>
  operator()(const edm4hep::SimTrackerHitCollection& inputSimHits) const override;
  /**  Finalize.
   *   @return status code
   */
  StatusCode finalize() override;

private:
  // Detector name
  Gaudi::Property<std::string> m_detectorName{this, "detectorName", "Vertex", "Name of the detector (default: Vertex)"};
  // List of pixels sizes along Phi angle in mm per layer (for strips, this is the pitch)
  Gaudi::Property<std::vector<float>> m_PixSizeX{
      this,
      "PixSizePhi",
      {0.020},
      "List of pixels sizes along Phi angle in mm per layer (for strips, this is the pitch)"};
  // List of pixels sizes along Theta angle in mm per layer (for strips, this is the pitch) (Theta is along BeamPipe for
  // Barrel and along ring radius for Endcaps)
  Gaudi::Property<std::vector<float>> m_PixSizeY{
      this,
      "PixSizeTheta",
      {0.020},
      "List of pixels sizes along Theta angle in mm per layer (for strips, this is the pitch)"};
  // Length of ionization steps in millimeters
  Gaudi::Property<float> m_IonizationGranularity{this, "IonizationGranularity", 0.001,
                                                 "Length of ionization steps in millimeters (default is 1 micron)"};
  // Tangent of sensor's Lorentz angle (default is 0.1)
  Gaudi::Property<float> m_tanLorentzAnglePerTesla{this, "tanLorentzAnglePerTesla", 0.1f,
                                                   "Tangent of sensor's Lorentz angle per Tesla (default is 0.1)"};
  // Charge diffusion
  Gaudi::Property<float> m_Sigma50{
      this, "Sigma50", 0.00151, "Charge diffusion in millimeters for 50 micron Si. Default = 0.00151mm taken from CMS"};
  // t resolution in ns
  Gaudi::Property<std::vector<float>> m_t_resolution{this, "tResolution", {0.1}, "Time resolutions per layer [ns]"};
  // Threshold in electron
  Gaudi::Property<float> m_Threshold{this, "Threshold", 0.0, "Charge Threshold in e (default: 0.0)"};
  // Threshold smearing in electron
  Gaudi::Property<float> m_ThresholdSmearing{this, "ThresholdSmearing", 0.0,
                                             "Sigma of Threshold Gaussian Smearing in e (default: 0.0)"};

  // List of sensor thickness per layer in millimeter
  std::vector<float> m_sensorThickness;
  // Normal Vector direction in sensor local frame (normal to sensor surface) (initialised during the first hit in
  // primary_ionization -> should be the same in the full sub-detector)
  mutable std::string m_LocalNormalVectorDir = "";
  // Define 50microns in mm
  float m_Dist50 = 0.050;
  // Define the area of integration of the charge in number of sigma around the central position
  float m_ClusterWidth = 3.0;

  // Declare the geometry service (initialised in the constructor)
  SmartIF<IGeoSvc> m_geoSvc;
  // Decoder for the cellID
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  // Volume manager to get the physical cell sensitive volume
  dd4hep::VolumeManager m_volman;
  // Surface map from SurfaceManager
  const dd4hep::rec::SurfaceMap* surfaceMap;

  // Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;

  // Gaussian random number generator used for time smearing
  std::vector<Rndm::Numbers> m_gauss_t_vec;

  Rndm::Numbers m_gauss_threshold; // for threshold's smearing

  // Create a class to map the charge per pixel (x,y) in electrons
  using hit_map_type = std::map<int, std::map<int, float, std::less<int>>, std::less<int>>;

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
    float charge() const { return _charge; }

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

    SignalPoint(float x, float y, float sigma_x, float sigma_y, float t, const edm4hep::SimTrackerHit& hit,
                float a = 1.0)
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

private:
  // Additional member functions
  // Private methods
  template <typename T>
  void getSensorThickness();

  void GetNormalVectorLocal(const edm4hep::SimTrackerHit& input_sim_hit) const;

  std::vector<ChargeDepositUnit> primary_ionization(const edm4hep::SimTrackerHit& hit) const;

  std::vector<SignalPoint> drift(const edm4hep::SimTrackerHit& hit,
                                 const std::vector<ChargeDepositUnit>& ionizationPoints) const;

  dd4hep::rec::Vector3D DriftDirection(const edm4hep::SimTrackerHit& hit) const;

  hit_map_type get_charge_per_pixel(const edm4hep::SimTrackerHit& hit,
                                    const std::vector<SignalPoint>& collectionPoints) const;

  bool Apply_Threshold(double& ChargeInPixel) const;

  void generate_output(const edm4hep::SimTrackerHit& hit, edm4hep::TrackerHitPlaneCollection& output_digi_hits,
                       edm4hep::TrackerHitSimTrackerHitLinkCollection& output_sim_digi_link_col,
                       const hit_map_type& hit_map) const;

  void SetProperDirectFrame(TGeoHMatrix& sensorTransformMatrix) const;

  void SetLocalPos_In_ProperDirectFrame(float& x, float& y, float& z) const;

private:
  // Additional Debug information

  Gaudi::Property<bool> m_DebugHistos{this, "DebugHistos", false,
                                      "If set to true, create a file containing debug histograms. Should then also set "
                                      "DebugFileName to set the name of the debug root file."};
  Gaudi::Property<std::string> m_DebugFileName{this, "DebugFileName", "",
                                               "Name of the file containing debug histograms."};

  TH1D* hErrorX; // Histogram to store the distance in X between the true hit and digitized one in mm
  TH1D* hErrorY; // Histogram to store the distance in Y between the true hit and digitized one in mm
  TH1D* hErrorZ; // Histogram to store the distance in Z between the true hit and digitized one in mm
  TH1D* hError;  // Histogram to store the distance between the true hit and digitized one in mm

  TH1D* hChargeAboveThreshold;  // Histogram to store the pixel charge after Threshold
  TH1D* hChargeBeforeThreshold; // Histogram to store the pixel charge before Threshold
  TH1D* hChargePerCluster;      // Histogram to store charge per Cluster

  TH1D* hActivePixelCountBeforeThreshold; // Histogram to store the number of active pixels per Cluster before Threshold
  TH1D* hActivePixelCountAfterThreshold;  // Histogram to store the number of active pixels per Cluster after Threshold

  TH1D* hXDriftDueToMagField; // Histogram to store the X drift due to magnetic field in mm
  TH1D* hYDriftDueToMagField; // Histogram to store the Y drift due to magnetic field in mm

  TH1D* hClusterPerLayer;     // Histogram to store the number of cluster per layer
  TH1D* hActivePixelPerlayer; // Histogram to store the amount of active pixels per layer

  void Create_outputROOTfile_for_debugHistograms() const;
};