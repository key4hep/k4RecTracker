/** ======= DCHdigi_v01 ==========
 * Gaudi Algorithm for DCH digitization
 *
 *
 * @author Alvaro Tolosa-Delgado, Brieuc Francois
 * @date   2024-08
 *
 * <h4>Input collections and prerequisites</h4>
 * Processor requires a collection of SimTrackerHits <br>
 * This code uses DD4hep length natural unit (cm), but EDM4hep data is (usually) in mm. Please be careful with units.  <br>
 * <h4>Output</h4>
 * Processor produces collection of Digitized hits of Drift Chamber v2<br>
 * @param DCH_simhits The name of input collection, type edm4hep::SimTrackerHitCollection <br>
 * (default name empty) <br>
 * @param DCH_DigiCollection The name of out collection, type extension::SenseWireHitCollection <br>
 * (default name DCH_DigiCollection) <br>
 * @param DCH_name DCH subdetector name <br>
 * (default value DCH_v2) <br>
 * @param calculate_dndx Optional flag to calcualte dNdx information <br>
 * (default value false) <br>
 * @param fileDataAlg File needed for calculating cluster count and size <br>
 * (default value /eos/.../DataAlgFORGEANT.root) <br>
 * @param zResolution_mm Resolution (sigma for gaussian smearing) along the sense wire, in mm <br>
 * (default value 1 mm) <br>
 * @param xyResolution_mm Resolution (sigma for gaussian smearing) perpendicular the sense wire, in mm <br>
 * (default value 0.1 mm) <br>
 * @param create_debug_histograms Optional flag to create debug histograms <br>
 * (default value false) <br>
 * @param GeoSvcName Geometry service name <br>
 * (default value GeoSvc) <br>
 * @param uidSvcName The name of the UniqueIDGenSvc instance, used to create seed for each event/run, ensuring reproducibility. <br>
 * (default value uidSvc) <br>
 * <br>
 */

#ifndef DCHDIGI_V01_H
#define DCHDIGI_V01_H

// Gaudi Transformer baseclass headers
#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

// Gaudi services
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// EDM4HEP
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/ParticleIDData.h"
#include "edm4hep/SimTrackerHitCollection.h"

// EDM4HEP extension
#include "extension/SenseWireHitCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"

// DD4hep
#include "DD4hep/Detector.h"  // for dd4hep::VolumeManager
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

// Class developed by Walaa for the CLS
#include "AlgData.h"

/// constant to convert from mm (EDM4hep) to DD4hep (cm)

struct DCHdigi_v01 final
    : k4FWCore::MultiTransformer<
          std::tuple<extension::SenseWireHitCollection, extension::SenseWireHitSimTrackerHitLinkCollection>(
              const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&)> {
  DCHdigi_v01(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<extension::SenseWireHitCollection, extension::SenseWireHitSimTrackerHitLinkCollection>
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

  /// Decoder for the cellID
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

  /// Pointer to drift chamber data extension
  dd4hep::rec::DCH_info* dch_data = {nullptr};

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

  //------------------------------------------------------------------
  //        ancillary functions

  bool IsFileGood(std::string& ifilename) const { return std::ifstream(ifilename).good(); }

  /// Print algorithm configuration
  void PrintConfiguration(std::ostream& io);

  /// Send error message to logger and then throw exception
  void ThrowException(std::string s) const;

  int CalculateLayerFromCellID(dd4hep::DDSegmentation::CellID id) const {
    //return m_decoder->get(id, "layer") + dch_data->nlayersPerSuperlayer * m_decoder->get(id, "superlayer") + 1;
    return dch_data->CalculateILayerFromCellIDFields( m_decoder->get(id, "layer"), m_decoder->get(id, "superlayer") );
  }

  int CalculateNphiFromCellID(dd4hep::DDSegmentation::CellID id) const { return m_decoder->get(id, "nphi"); }

  TVector3 Convert_EDM4hepVector_to_TVector3(const edm4hep::Vector3d& v, double scale) const {
    return {v[0] * scale, v[1] * scale, v[2] * scale};
  };
  edm4hep::Vector3d Convert_TVector3_to_EDM4hepVector(const TVector3& v, double scale) const {
    return {v.x() * scale, v.y() * scale, v.z() * scale};
  };

  //------------------------------------------------------------------
  //        cluster calculation, developed by Walaa
  
  /// Flag to create to calculate cluster counting information
  Gaudi::Property<bool> m_calculate_dndx{this, "calculate_dndx", false,
                                              "Calculate number of clusters and electron per cluster"};
                                              
  /// file with distributions to be sampled
  Gaudi::Property<std::string> m_fileDataAlg{
      this, "fileDataAlg", "/eos/project/f/fccsw-web/www/filesForSimDigiReco/IDEA/DataAlgFORGEANT.root",
      "ROOT file with cluster size distributions"};

  /// pointer to wrapper class, which contains the cluster size and number distributions
  AlgData* flData;

  /// code developed by Walaa for calculating number of clusters and cluster size of each one
  std::pair<uint32_t, std::vector<int>> CalculateClusters(const edm4hep::SimTrackerHit& input_sim_hit, TRandom3 * myRandom) const;

  bool IsParticleCreatedInsideDriftChamber(const edm4hep::MCParticle &) const ;


  //------------------------------------------------------------------
  //        debug information

  /// Flag to create output file with debug histgrams
  Gaudi::Property<bool> m_create_debug_histos{this, "create_debug_histograms", false,
                                              "Create output file with histograms for debugging"};

  /// name for the file that will contain the histograms for debugging
  Gaudi::Property<std::string> m_out_debug_filename{this, "out_debug_filename", "dch_digi_alg_debug.root",
                                                    "name for the file that will contain the histograms for debugging"};
  /// histogram to store distance from sim hit position to the sense wire
  TH1D* hDpw;

  /// histogram to store distance from digi-hit to the wire. Should be zero because digi-hit central position lies on the wire.
  /// This histogram is a consistency check, because the function used to calculate the distance to the wire is different from
  /// the function used to calculate the digi-hit central position from a sim-hit position
  TH1D* hDww;

  /// histogram to store smearing along the wire
  TH1D* hSz;

  /// histogram to store smearing perpendicular the wire
  TH1D* hSxy;

  /// Create ROOT file for debug histograms
  /// Does not change ROOT directory
  void Create_outputROOTfile_for_debugHistograms();
};

DECLARE_COMPONENT(DCHdigi_v01);

#endif
