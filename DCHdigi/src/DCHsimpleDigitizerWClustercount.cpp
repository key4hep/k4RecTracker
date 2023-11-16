#include "DCHsimpleDigitizerWClustercount.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDRec/Vector3D.h"

// ROOT
#include "Math/Cylindrical3D.h"
#include "AlgData.h" 
#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "edm4hep/ParticleIDData.h"
#include "TRandom.h"
DECLARE_COMPONENT(DCHsimpleDigitizerWClustercount)

DCHsimpleDigitizerWClustercount::DCHsimpleDigitizerWClustercount(const std::string& aName, ISvcLocator* aSvcLoc)
    : GaudiAlgorithm(aName, aSvcLoc), m_geoSvc("GeoSvc", "DCHsimpleDigitizerWClustercount") {
  declareProperty("inputSimHits", m_input_sim_hits, "Input sim tracker hit collection name");
  declareProperty("outputDigiHits", m_output_digi_hits, "Output digitized tracker hit collection name");
}

DCHsimpleDigitizerWClustercount::~DCHsimpleDigitizerWClustercount() {}

StatusCode DCHsimpleDigitizerWClustercount::initialize() {
  // Initialize random services
  if (service("RndmGenSvc", m_randSvc).isFailure()) {
    error() << "Couldn't get RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_gauss_z.initialize(m_randSvc, Rndm::Gauss(0., m_z_resolution)).isFailure()) {
    error() << "Couldn't initialize RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_gauss_xy.initialize(m_randSvc, Rndm::Gauss(0., m_xy_resolution)).isFailure()) {
    error() << "Couldn't initialize RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }

  // check if readout exists
  if (m_geoSvc->lcdd()->readouts().find(m_readoutName) == m_geoSvc->lcdd()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  // set the cellID decoder
  m_decoder = m_geoSvc->lcdd()->readout(m_readoutName).idSpec().decoder();
  // retrieve the volume manager
  m_volman = m_geoSvc->lcdd()->volumeManager();

////Added by Walaa for the CLS///////////
  TString DataAlg= "/afs/cern.ch/work/w/welmeten/public/IDEA_Key4HEP_Jun22/Cluster_sim_code/DataAlgFORGEANT.root";
  debug() << "root file opened"<<endmsg;
  flData = new AlgData();
  flData->Load_file(DataAlg.Data());
  flData->Load_interp();
////////////////////////////////////////
  return StatusCode::SUCCESS;
}

StatusCode DCHsimpleDigitizerWClustercount::execute() {
  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  debug() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Digitize the sim hits
  extension::DriftChamberDigiCollection* output_digi_hits = m_output_digi_hits.createAndPut();
  for (const auto& input_sim_hit : *input_sim_hits) {
    auto output_digi_hit = output_digi_hits->create();
    // smear the hit position: need to go in the wire local frame to smear in the direction aligned/perpendicular with the wire for z/distanceToWire, taking e.g. stereo angle into account
    // retrieve the cell detElement
    dd4hep::DDSegmentation::CellID cellID         = input_sim_hit.getCellID();
    auto                           cellDetElement = m_volman.lookupDetElement(cellID);
    // retrieve the wire (in DD4hep 1.23 there is no easy way to access the volume daughters we have to pass by detElements, in later versions volumes can be used)
    const std::string& wireDetElementName =
        Form("superLayer_%d_layer_%d_phi_%d_wire", m_decoder->get(cellID, "superLayer"),
             m_decoder->get(cellID, "layer"), m_decoder->get(cellID, "phi"));
    dd4hep::DetElement wireDetElement = cellDetElement.child(wireDetElementName);
    // get the transformation matrix used to place the wire
    const auto& wireTransformMatrix = wireDetElement.nominal().worldTransformation();
    // Retrieve global position in mm and apply unit transformation (translation matrix is stored in cm)
    double simHitGlobalPosition[3] = {input_sim_hit.getPosition().x * dd4hep::mm,
                                      input_sim_hit.getPosition().y * dd4hep::mm,
                                      input_sim_hit.getPosition().z * dd4hep::mm};
    double simHitLocalPosition[3]  = {0, 0, 0};
    // get the simHit coordinate in cm in the wire reference frame to be able to apply smearing of radius perpendicular to the wire
    wireTransformMatrix.MasterToLocal(simHitGlobalPosition, simHitLocalPosition);
    debug() << "Cell ID string: " << m_decoder->valueString(cellID) << endmsg;
    debug() << "Global simHit x " << simHitGlobalPosition[0] << " --> Local simHit x " << simHitLocalPosition[0]
            << " in cm" << endmsg;
    debug() << "Global simHit y " << simHitGlobalPosition[1] << " --> Local simHit y " << simHitLocalPosition[1]
            << " in cm" << endmsg;
    debug() << "Global simHit z " << simHitGlobalPosition[2] << " --> Local simHit z " << simHitLocalPosition[2]
            << " in cm" << endmsg;
    // build a vector to easily apply smearing of distance to the wire
    dd4hep::rec::Vector3D simHitLocalPositionVector(simHitLocalPosition[0], simHitLocalPosition[1],
                                                    simHitLocalPosition[2]);
    // get the smeared distance to the wire (cylindrical coordinate as the smearing should be perpendicular to the wire)
    debug() << "Original distance to wire: " << simHitLocalPositionVector.rho() << endmsg;
    double smearedDistanceToWire = simHitLocalPositionVector.rho() + m_gauss_xy.shoot() * dd4hep::mm;
    while(smearedDistanceToWire < 0){
      debug() << "Negative smearedDistanceToWire (" << smearedDistanceToWire << ") shooting another random number" << endmsg;
      smearedDistanceToWire = simHitLocalPositionVector.rho() + m_gauss_xy.shoot() * dd4hep::mm;
    }
    debug() << "Smeared distance to wire: " << smearedDistanceToWire << endmsg;
    // smear the z position (in local coordinate the z axis is aligned with the wire i.e. it take the stereo angle into account);
    double smearedZ = simHitLocalPositionVector.z() + m_gauss_z.shoot() * dd4hep::mm;
    // fill the output DriftChamberDigi (making sure we are back in mm)
    output_digi_hit.setCellID(cellID);
    output_digi_hit.setDistanceToWire(smearedDistanceToWire / dd4hep::mm);
    output_digi_hit.setZPositionAlongWire(smearedZ / dd4hep::mm);

 //_________________SET NECESSARY PARAMETERS FOR THE CLS ALGORITHM-----WALAA_________________//
  
     Momentum = sqrt((input_sim_hit.getMomentum().x * input_sim_hit.getMomentum().x) + (input_sim_hit.getMomentum().y * input_sim_hit.getMomentum().y) + (input_sim_hit.getMomentum().z * input_sim_hit.getMomentum().z)) ;
     
    NCld =0;
    int ParentID = input_sim_hit.getCellID();
    int pdgid = static_cast<int>(input_sim_hit.getTime());
    
    if(ParentID==1&&pdgid==11){
        NCld++;
		Ekdelta=(TMath::Sqrt(Momentum*Momentum+me*me)-me)*1e6;
		vecExtraD.push_back(Ekdelta);
        }
    debug() << "NCld= "<<NCld<<"vecExtraD= "<<vecExtraD<<endmsg;
    if (NCld >0)    debug() << "There is delta rays"<<endmsg;
    TParticlePDG* partRoot=TDatabasePDG::Instance()->GetParticle(pdgid);
	if(partRoot!=nullptr){
	parName=partRoot->GetName();
	parClass=partRoot->ParticleClass();
	debug() << "parName "<<parName<<" parClass= "<<parClass<<endmsg;	
		mass=partRoot->Mass()*1e+3; //MeV
		bg=Momentum/mass;
	} 
	
 Ffrac=flData->get_Ffrac(bg);//0.9;
 Fmpv1=flData->get_Fmpv1(bg);//43;
 Fsgm1=flData->get_Fsgm1(bg);//9.64;
 Fmpv2=flData->get_Fmpv2(bg);//311.3;
 Fsgm2=flData->get_Fsgm2(bg);//35;
 CorrpMean=flData->get_ClSzCorrpmean(bg);
 CorrpSgm=flData->get_ClSzCorrpsgm(bg);
 Corrdgmean=flData->get_ClSzCorrdgmean(bg);
 Corrdgsgm=flData->get_ClSzCorrdgsgm(bg);
 Corrdglfrac=flData->get_ClSzCorrdglfrac(bg);
 Corrdglmpvl=flData->get_ClSzCorrdglmpvl(bg);
 Corrdglsgml=flData->get_ClSzCorrdglsgml(bg);
 Corrdglmeang=flData->get_ClSzCorrdglmeang(bg);
 Corrdglsgmg=flData->get_ClSzCorrdglsgmg(bg);
 maxEx0=flData->get_maxEx0(bg);//379.4 electron 100 gev
 maxExSlp=flData->get_maxExSlp();
 if(parClass=="Lepton"){
 ExSgmlep=flData->get_ExSgmlep();
 }else{
  ExSgmhad=flData->get_ExSgmhad();
 }

 MPVEx=flData->get_MPVExtra(bg);//103.2
 SgmEx=flData->get_SgmExtra(bg);//28.5;//
 MeanEx1=flData->get_MeanExtra1(bg);//13.8;//
 SgmEx1=flData->get_SgmExtra1(bg);//9.72;//
 frac=flData->get_FracExtra1(bg);//0.13;//
 Slp=flData->get_SlopeExtra1(bg);//7.95;//
debug() <<  " maxEx0 "<< maxEx0<< " ExSgmhad "<<ExSgmhad<<" MPVEx "<<MPVEx<<" SgmEx "<<SgmEx<<endmsg;
 CorrSlp=flData->get_ClSzCorrSlp(bg);
 CorrInt=flData->get_ClSzCorrInt(bg);

  vecpar[0]=Ffrac;
  vecpar[1]=Fmpv1;
  vecpar[2]=Fsgm1;
  vecpar[3]=Fmpv2;
  vecpar[4]=Fsgm2;
  
double Tmax = (2.0*me*pow(bg,2)/(1+(2.0*(1+pow(bg,2))*me/mass)+pow(me/mass,2)))*1e+6;
float maxEcut=cut;
if(Tmax<maxEcut){maxEcut=Tmax;}
//cout << "Tmax ="<<Tmax<<"maxEcut= "<<maxEcut<<"cut= "<<maxEcut<<endl;
/*int nCv=20;
fcnvl[0] = new TF1("fcnvl-0","[0]*TMath::Landau(x,[1],[2],true)+(1.0-[0])*TMath::Landau(x,[3],[4],true)");
fcnvl[0]->SetRange(0,maxEcut);
fcnvl[0]->SetParameters(vecpar[0],vecpar[1],vecpar[2],vecpar[3],vecpar[4]);
cnvl[0]=nullptr;*/


	TF1 *land= new TF1("land","landaun");
	land->SetParameter(0,1);
	land->SetParameter(1,MPVEx);
	land->SetParameter(2,SgmEx);
	land->SetRange(0,maxEcut);

	TF1 *exGauss=new TF1("exGauss","[0]*([1]*TMath::Gaus(x,[2],[3],true)+(1.0-[1])*TMath::Exp(-x/[4])/[4])");
	exGauss->SetParameter(0,1);
	exGauss->SetParameter(1,frac);
	exGauss->SetParameter(2,MeanEx1);
	exGauss->SetParameter(3,SgmEx1);
	exGauss->SetParameter(4,Slp);
	exGauss->SetRange(0,90);


//	int MaxEv = Entries > Nev ? Nev : Entries;
//	cout<<" MAXEV"<<MaxEv<<endl;
	sumVecTot=meanVecTot=sumVec=meanVec=0.0;
	vecEtr.clear();
	vecEtrTot.clear();
  //_________________End SET NECESSARY PARAMETERS FOR ALGORITHM_________________//
  
    NCl=0;
	NCl1=0;
	NClp=0;
	totExECl=0.0;
	ExECl=0.0;
    LengthTrack =  input_sim_hit.getPathLength ();
    Etot_per_track = input_sim_hit.getEDep (); 
	LengthTrack*=0.1;//from mm to cm
	Eloss=Etot_per_track*1.e09; //from GeV to eV
	maxEx0len=maxEx0*LengthTrack;    //maxEx0 is a parameter
	if(parClass=="Lepton"){
		ExSgmlen=ExSgmlep*TMath::Sqrt(LengthTrack);
	}else{
		ExSgmlen=ExSgmhad*TMath::Sqrt(LengthTrack);
	}
	maxExECl=(Eloss-maxEx0len+gRandom->Gaus(0,ExSgmlen))/maxExSlp;
 //  cout<<"maxExECl=   "<<maxExECl<<endl;
//	SgmaxExECl=ExSgm/maxExSlp;
	if(maxExECl<EIzs) {maxExECl=0.0;} //EIzs const = 25.6
		
	//else if(choice==2){
	if(choice==2){
		//	cout<<"start new"<<"Eloss= "<<Eloss<<"EIzp=  "<<EIzp<< "EIzs= "<< EIzs<<"maxExECl=  "<<maxExECl<<"totExECl= "<<totExECl<<endl;
			debug() << "Eloss= "<<Eloss<<"EIzs= "<<EIzs<<"EIzp= "<<EIzp<<"maxExECl= "<<maxExECl<<"totExECl"<<totExECl<<endmsg;

			while(Eloss>(EIzp+EIzs)&&maxExECl>totExECl){
	//		cout<<"start again"<<endl;
	//		cout<<"Eloss= "<<Eloss<<"EIzp=  "<<EIzp<< "EIzs= "<< EIzs<<"maxExECl=  "<<maxExECl<<"totExECl= "<<totExECl<<endl;
				ExECl=land->GetRandom(0,maxEcut);
//				cout<<"ExECl= "<<ExECl<<endl;
//				ExECl=fcnvl[0]->GetRandom();
			//	hExtraCl->Fill(ExECl);
			
				if(ExECl>EIzs){
					Eloss-=EIzp;
					if(ExECl>(maxExECl-totExECl)) {
						ExECl=maxExECl-totExECl;
						//							cout<< " count "<<count<<endl;
					}
					if(ExECl>Eloss) ExECl=Eloss;
					totExECl+=ExECl;
	//				cout<<"totExECl= " <<totExECl<<endl;
					Eloss-=ExECl;
	//				cout<<"Eloss= " <<Eloss<<endl;
				
					NClp++;
	//				cout<<"ExECl= "<<ExECl<<"totExECl= "<<totExECl<<"Eloss= "<<Eloss<<"NClp= " <<NClp<<endl;
					//CLSZ
					float tmpCorr=0.0;
					for(int i=0;i<nhEp;++i){
						if(ExECl>=(i==0 ? 0 :hEpcut[i-1])&&ExECl<hEpcut[i]){
							tmpCorr=gRandom->Gaus(CorrpMean[i],CorrpSgm[i]);
						}
					}
					ClSz=TMath::Nint(ExECl*CorrSlp+CorrInt-tmpCorr);
					if(ClSz<2) {ClSz=2;}
					output_digi_hit.setClusterSize(ClSz);
			//		hClSzRec->Fill(ClSz);
			//		hClSzRecND->Fill(ClSz);
					NEltot+=ClSz;
         //           fClusterCharge.push_back(ClSz);
			//	    cout<<"ClSz=  "<<"fClusterCharge= "<<fClusterCharge<<endl;
				}
//				cout<<"try1= "<<endl;
			}
  //             cout<<"try2= "<<endl;
  //             cout<<"Eloss= "<<Eloss<<"EIzp=  "<<EIzp<< "EIzs= "<< EIzs<<"maxExECl=  "<<maxExECl<<"totExECl= "<<totExECl<<endl;

		}	
	debug() << "Eloss= "<<Eloss<<"EIzp= "<<EIzp<<endmsg;

	while(Eloss>=EIzp){

			Eloss-=EIzp;
			ExECl1=exGauss->GetRandom();
		//	hExtraCl1->Fill(ExECl1);
			if(ExECl1>Eloss) {ExECl1=Eloss;}
			ExECl1totRec+=ExECl1;
			NCl1++;
			//				EIzp=hEIzp->GetRandom();
			Eloss-=ExECl1;
			ClSz=1;
			output_digi_hit.setClusterSize(ClSz);
		//	hClSzRec->Fill(ClSz);
		//	hClSzRecND->Fill(ClSz);
		 //   fClusterCharge.push_back(ClSz);
			NEltot+=ClSz;

		}		

	float tmpCl;
		for(int i=0;i<vecExtraD.size();++i){
			int tmphE=(vecExtraD[i]-minE)/binE;
			if(tmphE>=nhE) tmphE=nhE-1;
//			cout<<" tmphE "<<tmphE<<endl;
			if(tmphE==nhE-1){
				rndCorr=gRandom->Uniform(0,1);
				if(rndCorr<Corrdglfrac[tmphE]){
					tmpCl=gRandom->Gaus(Corrdglmeang[tmphE],Corrdglsgmg[tmphE]);
				}else{
					tmpCl=gRandom->Landau(Corrdglmpvl[tmphE],Corrdglsgml[tmphE]);
				}
			}else{
				tmpCl=gRandom->Gaus(Corrdgmean[tmphE],Corrdgsgm[tmphE]);
			}

			ClSz=TMath::Nint(vecExtraD[i]*CorrSlp+CorrInt-tmpCl);
			if(ClSz<2) {ClSz=2;}
			output_digi_hit.setClusterSize(ClSz);
	//		hClSzRec->Fill(ClSz);
	//		hClSzDRec->Fill(ClSz);
	  //      fClusterCharge.push_back(ClSz);
		 	NEltot+=ClSz;

		}
		
	  Ncl=NCl1+NClp+NCld;
	debug() <<"Ncl= "<< Ncl<<" NCl1= "<<NCl1<<"NClp= "<<NClp<<"NCld= "<<NCld<<endmsg;
	debug() <<"NCld= "<< NCld<<" vecExtraD.size() "<<vecExtraD.size()<< endmsg;
	  if (NCld !=vecExtraD.size()) debug() <<"There is a bug"<< endmsg;
  debug() << "Output Digi Hit collection size: " << output_digi_hits->size() << endmsg;
  
 output_digi_hit.setClusterCount(Ncl);
 
/* for (int icl=0;icl<Ncl;icl++) {
 output_digi_hit.setClusterSize(NEltot);
 }*/
// output_digi_hit.setClusterSize(NEltot);


  }
 

/////////////////////////////////////////////////////////////////
  return StatusCode::SUCCESS;
}

StatusCode DCHsimpleDigitizerWClustercount::finalize() { return StatusCode::SUCCESS; }
