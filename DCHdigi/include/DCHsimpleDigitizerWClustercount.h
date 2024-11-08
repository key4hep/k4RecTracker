#pragma once

// GAUDI
#include "Gaudi/Property.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

// K4FWCORE
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"

// EDM4HEP
#include "edm4hep/SimTrackerHitCollection.h"

// EDM4HEP extension
#include "extension/DriftChamberDigiCollection.h"

// DD4HEP
#include "DD4hep/Detector.h"  // for dd4hep::VolumeManager
#include "DDSegmentation/BitFieldCoder.h"
#include "TF1Convolution.h"
/** @class DCHsimpleDigitizerWClustercount
 *
 *  Algorithm for creating digitized drift chamber hits (extension::DriftChamberDigi) from edm4hep::SimTrackerHit.
 *  You have to specify the expected resolution in z and in xy (distance to the wire). The smearing is applied in the wire reference frame.
 *  
 *  @author Brieuc Francois
 *  @date   2023-03
 *
 */
class AlgData; //Added by Walaa for the CLS 

class DCHsimpleDigitizerWClustercount : public GaudiAlgorithm {
public:
  explicit DCHsimpleDigitizerWClustercount(const std::string&, ISvcLocator*);
  virtual ~DCHsimpleDigitizerWClustercount();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Execute.
   *   @return status code
   */
  virtual StatusCode execute() final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:
  // Input sim tracker hit collection name
  DataHandle<edm4hep::SimTrackerHitCollection> m_input_sim_hits{"inputSimHits", Gaudi::DataHandle::Reader, this};
  // Output digitized tracker hit collection name
  DataHandle<extension::DriftChamberDigiCollection> m_output_digi_hits{"outputDigiHits", Gaudi::DataHandle::Writer, this};

  // Detector readout name
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "CDCHHits", "Name of the detector readout"};
  // Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  // Decoder for the cellID
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  // Volume manager to get the physical cell sensitive volume
  dd4hep::VolumeManager m_volman;

  // z position resolution in mm
  FloatProperty m_z_resolution{this, "zResolution", 1.0,
                               "Spatial resolution in the z direction (from reading out the wires at both sides) [mm]"};
  // xy resolution in mm
  FloatProperty m_xy_resolution{this, "xyResolution", 0.1, "Spatial resolution in the xy direction [mm]"};

  // Random Number Service
  IRndmGenSvc* m_randSvc;
  // Gaussian random number generator used for the smearing of the z position
  Rndm::Numbers m_gauss_z;
  // Gaussian random number generator used for the smearing of the xy position
  Rndm::Numbers m_gauss_xy;
  
    // Parameters needed for the CL Algo //Added by Walaa///////
    ///from Createclusters.cpp////
  Int_t fNClusters;
  double me = 0.511;
  float rnd;
  float Eloss= 0.0; float DeltaE_per_track= 0.0; float EtotCell= 0.0;
  int NCl,NCl1,NClp;
  float ClSz,ClSzP;
  float Rt=0.87;
  float EIzp=15.8;
  float ExECl1;
  float ExECl;
  float cut= 1000;//controlla
  float prc=0.83;
  float EIzs=25.6;
  float ratio=prc/EIzs;
  int NCltot=0;
  float Edep=0.0;
  float DeltaE=0.0;
  int count;
  bool loop;
  float prvDiff;
  float tmpDiff[10];
  float vecExECl[10];
  float mintmpDiff;
  int iloop;
  float ExECl1totRec;
  float SgmaxExECl;
  int Ncl;  
  float rndCorr;
  float maxExECl=0.0; float ExSgm=0.0; 
  float totExECl;
//  int ParentID;
  const int nhEp=10;//10
  float hEpcut[10]={100,200,300,400,500,600,700,800,900,1000};
  int minE=1000;
  int maxE=10000;
  int binE=1000;
  int nhE=(maxE-minE)/binE;
  TString parClass;
  int Ncltot=0;
  int NEltot=0;
  float Etot=0.0;
  float maxEx0len,ExSgmlen;
  int choice= 2;
  int Nev=500000;
  Int_t val;
  float EtCut=0.7;
 ///////end from Createclusters.cpp////
  
///////from Createclusters.h/////////
  TF1Convolution **cnvl =new TF1Convolution*[100];//100 for 1 cm gas box //Added by Walaa
  TF1 **fcnvl=new TF1 *[100];  //Added by Walaa
  float fnorm[100]; //Added by Walaa
  double vecpar[5]; //Added by Walaa
  std::vector <float> vecEtr,vecEtrTot;
  float sumVecTot,meanVecTot,sumVec,meanVec;
/////////end from Createclusters.h/////////
   int NCld;
   TString parName;
   double bg, Momentum, mass;
   AlgData *flData;
   std::vector<float> vecExtraD;

   double Ffrac, Fmpv1, Fsgm1, Fmpv2, Fsgm2;
   std::vector<double> CorrpMean,CorrpSgm,Corrdgmean,Corrdgsgm,Corrdglfrac,Corrdglmpvl,Corrdglsgml,Corrdglmeang,Corrdglsgmg;
   float Ekdelta;
   float maxEx0, maxExSlp, ExSgmlep, ExSgmhad;
   float MPVEx, SgmEx, MeanEx1, SgmEx1, frac, Slp, CorrSlp, CorrInt;
   Double_t LengthTrack, Etot_per_track;
   //////end parameters

};
