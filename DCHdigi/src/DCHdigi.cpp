#include "DCHdigi.h"

// STL
#include <iostream>
#include <sstream>

#include "extension/MutableDriftChamberDigiV2.h"

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       DCHdigi constructor       ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// -- KeyValues("name of the variable that holds the name of the collection exposed in the python steering file", {"default name for the collection"}),
DCHdigi::DCHdigi(const std::string& name, ISvcLocator* svcLoc)
: MultiTransformer(name, svcLoc,
    {
        KeyValues("DCH_simhits", {""}),
        KeyValues("HeaderName", {"EventHeader"}),
    },
    {
        KeyValues("DCH_DigiCollection", {"DCH_DigiCollection"}),
        KeyValues("DCH_DigiSimAssociationCollection", {"DCH_DigiSimAssociationCollection"})
    }
  )
{
    m_geoSvc = serviceLocator()->service(m_geoSvcName);
    m_uidSvc = serviceLocator()->service(m_uidSvcName);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       initialize       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode DCHdigi::initialize() {

    if (!m_uidSvc)
    ThrowException( "Unable to get UniqueIDGenSvc" );

    m_gauss_z_cm  = std::normal_distribution<double>(0., m_z_resolution.value()*MM_TO_CM );
    m_gauss_xy_cm = std::normal_distribution<double>(0., m_xy_resolution.value()*MM_TO_CM);

    //-----------------
    // Retrieve the subdetector
    std::string DCH_name(m_DCH_name.value());
    if ( 0 == m_geoSvc->getDetector()->detectors().count(DCH_name) )
    {
        ThrowException( "Detector <<" + DCH_name + ">> does not exist." );
    }

    dd4hep::DetElement DCH_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);

    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  retrieve data extension     //////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    this->dch_data = DCH_DE.extension<dd4hep::rec::DCH_info>();

    ///////////////////////////////////////////////////////////////////////////////////

    //-----------------
    // Retrieve the readout associated with the detector element (subdetector)
    dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector(DCH_name);
    if(not dch_sd.isValid() )
    ThrowException("No valid Sensitive Detector was found for detector <<" + DCH_name + ">>.");

    dd4hep::Readout dch_readout = dch_sd.readout();
    // set the cellID decoder
    m_decoder = dch_readout.idSpec().decoder();

    ///////////////////////////////////////////////////////////////////////////////////
    /////////////////////////  initialize Walaa's code for CLS  ///////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
	debug() << Form("Opening %s...",m_fileDataAlg.value().c_str())<<endmsg;
    if( not IsFileGood(m_fileDataAlg.value()) ) ThrowException("File <<" + m_fileDataAlg.value() + ">> not found.");
	flData = new AlgData();
	flData->Load_file(m_fileDataAlg.value().c_str());
	flData->Load_interp();
	debug() << Form("Opening %s... Done",m_fileDataAlg.value().c_str())<<endmsg;
    ///////////////////////////////////////////////////////////////////////////////////

    std::stringstream ss;
    PrintConfiguration(ss);
    info() << ss.str().c_str() <<endmsg;
    if( m_create_debug_histos.value() )
    {
        hDpw = new TH1D("hDpw", "Distance hit to the wire, in cm", 100,0,1);
        hDpw->SetDirectory(0);
        hDww = new TH1D("hDww", "Distance hit projection to the wire, in cm. Should be zero", 100,0,1);
        hDww->SetDirectory(0);
        hSz  = new TH1D("hSz", "Smearing along the wire, in cm", 100,0,5*m_z_resolution.value());
        hSz->SetDirectory(0);
        hSxy = new TH1D("hSxy", "Smearing perpendicular the wire, in cm", 100,0,5*m_xy_resolution.value());
        hSxy->SetDirectory(0);
    }
    return StatusCode::SUCCESS;
  }



///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       operator()       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
std::tuple<colltype_out,colltype_out2>
DCHdigi::operator()(const colltype_in& input_sim_hits,
    const edm4hep::EventHeaderCollection&  headers) const {

    // initialize seed for random engine
    this->PrepareRandomEngine( headers );

    debug() << "Input Sim Hit collection size: " << input_sim_hits.size() << endmsg;

    // Create the collections we are going to return
    colltype_out output_digi_hits;
    colltype_out2 output_digi_sim_association;

    //loop over hit collection
    for (const auto& input_sim_hit : input_sim_hits)
    {
        dd4hep::DDSegmentation::CellID cellid = input_sim_hit.getCellID();
        int ilayer = this->CalculateLayerFromCellID(cellid );
        int nphi   = this->CalculateNphiFromCellID(cellid );
        auto hit_position =  Convert_EDM4hepVector_to_TVector3( input_sim_hit.getPosition(), MM_TO_CM );

        // -------------------------------------------------------------------------
        //      calculate hit position projection into the wire
        TVector3 hit_to_wire_vector = this->Calculate_hitpos_to_wire_vector(ilayer, nphi,hit_position);
        TVector3 hit_projection_on_the_wire = hit_position + hit_to_wire_vector;
        if( m_create_debug_histos.value() )
        {
            double distance_hit_wire = hit_to_wire_vector.Mag();
            hDpw->Fill(distance_hit_wire);
        }
        TVector3 wire_direction_ez = this->Calculate_wire_vector_ez(ilayer, nphi);

        // -------------------------------------------------------------------------
        //       smear the position

        //       smear position along the wire
        double smearing_z = m_gauss_z_cm( m_engine );
        if( m_create_debug_histos.value() ) hSz->Fill( smearing_z );

        hit_projection_on_the_wire +=  smearing_z*(wire_direction_ez.Unit());
        if( m_create_debug_histos.value() )
        {
            // the distance from the hit projection and the wire should be zero
            TVector3 dummy_vector = this->Calculate_hitpos_to_wire_vector(ilayer, nphi,hit_projection_on_the_wire);
            hDww->Fill( dummy_vector.Mag() );
        }

        //       smear position perpendicular to the wire
        double smearing_xy = m_gauss_xy_cm( m_engine );
        if( m_create_debug_histos.value() ) hSxy->Fill( smearing_xy );
        float distanceToWire_real = hit_to_wire_vector.Mag();

        // protect against negative values
        float distanceToWire_smeared = std::max(0.0, distanceToWire_real +smearing_xy );

        std::int32_t type = 0;
        std::int32_t quality = 0;
        float eDepError =0;
        // length units back to mm
        auto positionSW   = Convert_TVector3_to_EDM4hepVector(hit_projection_on_the_wire, 1./MM_TO_CM );
        auto directionSW  = Convert_TVector3_to_EDM4hepVector(wire_direction_ez         , 1./MM_TO_CM );
        float distanceToWire = distanceToWire_smeared/MM_TO_CM;

        auto [clusterCount,clusterSize] = CalculateClusters(input_sim_hit);

        extension::MutableDriftChamberDigiV2 oDCHdigihit;
        oDCHdigihit.setCellID(input_sim_hit.getCellID());
        oDCHdigihit.setType(type);
        oDCHdigihit.setQuality(quality);
        oDCHdigihit.setTime(input_sim_hit.getTime());
        oDCHdigihit.setEDep(input_sim_hit.getEDep());
        oDCHdigihit.setEDepError(eDepError);
        oDCHdigihit.setPosition(positionSW);
        oDCHdigihit.setDirectionSW(directionSW);
        oDCHdigihit.setDistanceToWire(distanceToWire);
        oDCHdigihit.setClusterCount(clusterCount);
        oDCHdigihit.setClusterSize(clusterSize);

        output_digi_hits.push_back(oDCHdigihit);

        extension::MutableMCRecoDriftChamberDigiV2Association oDCHsimdigi_association;
        oDCHsimdigi_association.setDigi( output_digi_hits.at( output_digi_hits.size()-1 )  );
        oDCHsimdigi_association.setSim( input_sim_hit  );

        }// end loop over hit collection

    /////////////////////////////////////////////////////////////////
    return std::make_tuple<colltype_out,colltype_out2>(std::move(output_digi_hits),std::move(output_digi_sim_association));
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       finalize       //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode DCHdigi::finalize()
{
    if( m_create_debug_histos.value() )
    {
        std::unique_ptr<TFile> ofile{TFile::Open ( "dch_digi_alg_debug.root", "recreate" ) };
        ofile->cd();
        hDpw->Write();
        hDww->Write();
        hSxy->Write();
        hSz->Write();
    }

    return StatusCode::SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       ThrowException       ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void DCHdigi::ThrowException(std::string s) const {
    error() << s.c_str()  << endmsg;
    throw std::runtime_error(s);
  }

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       PrintConfiguration       ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void DCHdigi::PrintConfiguration(std::ostream& io)
{
    io << "DCHdigi will use the following components:\n";
    io << "\tGeometry Service: "                  << m_geoSvcName.value().c_str()           << "\n";
    io << "\tUID Service: "                       << m_uidSvcName.value().c_str()           << "\n";
    io << "\tDetector name: "                     << m_DCH_name.value().c_str()             << "\n";
    io << "\t\t|--Volume bitfield: "              << m_decoder->fieldDescription().c_str()  << "\n";
    io << "\t\t|--Number of layers: "             << dch_data->database.size()              << "\n";
    io << "\tCluster distributions taken from: "  << m_fileDataAlg.value().c_str()          << "\n";
    io << "\tResolution along the wire (mm): "    << m_z_resolution.value()                 << "\n";
    io << "\tResolution perp. to the wire (mm): " << m_xy_resolution.value()                << "\n";
    return;
}

void DCHdigi::PrepareRandomEngine(const edm4hep::EventHeaderCollection&  headers) const
{
    uint32_t evt_n = headers[0].getEventNumber();
    uint32_t run_n = headers[0].getRunNumber();
    size_t seed = m_uidSvc->getUniqueID(evt_n, run_n, this->name() );
    m_engine.seed(seed);
    // test random engine...
    m_engine.discard(10);
    myRandom.SetSeed(seed+42);
}

///////////////////////////////////////////////////////////////////////////////////////
/////       Ancillary functions for calculating the distance to the wire       ////////
///////////////////////////////////////////////////////////////////////////////////////

double DCHdigi::Calculate_phi_rot_equivalent_to_hit_to_wire_distance(int ilayer, double hit_to_wire_distance) const
{
    auto & l = this->dch_data->database.at(ilayer);
    double rz0 = l.radius_sw_z0;
    return 2*atan( (hit_to_wire_distance/2.)/rz0 );
}



TVector3 DCHdigi::Calculate_wire_vector_ez(int ilayer, int nphi) const
{
    auto & l = this->dch_data->database.at(ilayer);

    // See original paper Hoshina et al, Computer Physics Communications 153 (2003) 3
    // eq. 2.9, for the definition of ez, vector along the wire

    // initialize some variables
    int stereosign = l.StereoSign();
    double rz0 = l.radius_sw_z0;
    double dphi = dch_data->twist_angle;
    // kappa is the same as in eq. 2.9
    double kappa = (1./dch_data->Lhalf)*tan(dphi/2);

    //--- calculating wire position
    // the points p1 and p2 correspond to the ends of the wire

    // point 1
    // double x1 = rz0; // m
    // double y1 = 0.; // m
    // double z1 = 0.; // m
    double x1 = rz0;                                        // m
    double y1 = -stereosign*rz0*kappa*dch_data->Lhalf;      // m
    double z1 = -dch_data->Lhalf;                           // m

    TVector3 p1 (x1,y1,z1);


    // point 2
    double x2 = rz0;                                        // m
    double y2 = stereosign*rz0*kappa*dch_data->Lhalf;       // m
    double z2 = dch_data->Lhalf;                            // m

    TVector3 p2 (x2,y2,z2);

    // calculate phi rotation of whole twisted tube, ie, rotation at z=0
    double phi_z0 = Calculate_wire_phi_z0(ilayer,nphi);
    p1.RotateZ(phi_z0);
    p2.RotateZ(phi_z0);

    //--- end calculating wire position

    return (p2-p1).Unit();

}

TVector3 DCHdigi::Calculate_wire_z0_point(int ilayer, int nphi) const
{
    auto & l = this->dch_data->database.at(ilayer);
    double rz0 = l.radius_sw_z0;
    TVector3 p1 (rz0,0,0);
    double phi_z0 = Calculate_wire_phi_z0(ilayer,nphi);
    p1.RotateZ(phi_z0);
    return p1;
}

// calculate phi rotation of whole twisted tube, ie, rotation at z=0
double DCHdigi::Calculate_wire_phi_z0(int ilayer, int nphi) const
{
    auto & l = this->dch_data->database.at(ilayer);
    int ncells = l.nwires/2;
    double phistep = TMath::TwoPi()/ncells;
    double phi_z0 = (nphi + 0.25*(l.layer%2))*phistep;
    return phi_z0;
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Calculate vector from hit position to wire   /////////////////
///////////////////////////////////////////////////////////////////////////////////////
TVector3 DCHdigi::Calculate_hitpos_to_wire_vector(int ilayer, int nphi, const TVector3 & hit_position /*in cm*/) const
{
    // Solution distance from a point to a line given here:
    // https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Vector_formulation
    TVector3 n = this->Calculate_wire_vector_ez(ilayer, nphi);
    TVector3 a = this->Calculate_wire_z0_point (ilayer, nphi);
    // Remember using cm as natural units of DD4hep consistently!
    // TVector3 p {hit_position.x()*MM_TO_CM,hit_position.y()*MM_TO_CM,hit_position.z()*MM_TO_CM};

    TVector3 a_minus_p = a - hit_position;
    double a_minus_p_dot_n = a_minus_p.Dot( n );
    TVector3 scaled_n = a_minus_p_dot_n * n;
    //hit_to_wire_vector = a_minus_p - scaled_n;
    return (a_minus_p - scaled_n);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       CalculateNClusters       ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
std::pair<uint32_t,uint32_t> DCHdigi::CalculateClusters(const edm4hep::SimTrackerHit & input_sim_hit) const
{
    std::pair<uint32_t,uint32_t> return_values = {0,0};
    uint32_t & Ncl  = return_values.first;
    uint32_t & ClSz = return_values.second;
    //_________________SET NECESSARY PARAMETERS FOR THE CLS ALGORITHM-----WALAA_________________//

    auto thisParticle = input_sim_hit.getParticle();
    // Parameters needed for the CL Algo //Added by Walaa///////
    ///from Createclusters.cpp////
    // int fNClusters;
  double me = 0.511;
  // float rnd;
  float Eloss= 0.0;
  // float DeltaE_per_track= 0.0;
  // float EtotCell= 0.0;
  int NCl;
  int NCl1,NClp;
  // float ClSz,ClSzP;
  // float Rt=0.87;
  float EIzp=15.8;
  float ExECl1;
  float ExECl;
  float cut= 1000;//controlla
  // float prc=0.83;
  float EIzs=25.6;
  // float ratio=prc/EIzs;
  // int NCltot=0;
  // float Edep=0.0;
  // float DeltaE=0.0;
  // int count;
  // bool loop;
  // float prvDiff;
  // float tmpDiff[10];
  // float vecExECl[10];
  // float mintmpDiff;
  // int iloop;
  float ExECl1totRec;
  // float SgmaxExECl;
  // int Ncl;
  float rndCorr;
  float maxExECl=0.0;
  // float ExSgm=0.0;
  float totExECl;
//  int ParentID;
  const int nhEp=10;//10
  float hEpcut[10]={100,200,300,400,500,600,700,800,900,1000};
  int minE=1000;
  int maxE=10000;
  int binE=1000;
  int nhE=(maxE-minE)/binE;
  TString parClass;
  // int Ncltot=0;
  int NEltot=0;
  // float Etot=0.0;
  float maxEx0len,ExSgmlen;
  int choice= 2;
  // int Nev=500000;
  // int val;
  // float EtCut=0.7;
 ///////end from Createclusters.cpp////

///////from Createclusters.h/////////
  // TF1Convolution **cnvl =new TF1Convolution*[100];//100 for 1 cm gas box //Added by Walaa
  // TF1 **fcnvl=new TF1 *[100];  //Added by Walaa
  // float fnorm[100]; //Added by Walaa
  double vecpar[5]; //Added by Walaa
  std::vector <float> vecEtr,vecEtrTot;
  float sumVecTot,meanVecTot,sumVec,meanVec;
/////////end from Createclusters.h/////////

  // get rid of spurious compiler error, Werror=unused-but-set-variable
  (void)NCl;
  (void)vecpar;
  (void)sumVecTot;

   int NCld;
   TString parName;
   double bg, Momentum, mass=-1e300;
   std::vector<float> vecExtraD;

   double Ffrac, Fmpv1, Fsgm1, Fmpv2, Fsgm2;
   std::vector<double> CorrpMean,CorrpSgm,Corrdgmean,Corrdgsgm,Corrdglfrac,Corrdglmpvl,Corrdglsgml,Corrdglmeang,Corrdglsgmg;
   float Ekdelta;
   float maxEx0, maxExSlp, ExSgmlep, ExSgmhad;
   float MPVEx, SgmEx, MeanEx1, SgmEx1, frac, Slp, CorrSlp, CorrInt;
   double LengthTrack, Etot_per_track;
   //////end parameters
		bool IsSecondaryWithinDCH = false;
		{
			auto vertex = thisParticle.getVertex(); // in mm
			auto vertexRsquared = vertex[0]*vertex[0] + vertex[1]*vertex[1];
			auto vertexZabs     = std::fabs(vertex[2]);
			float DCH_halflengh = 2000;
			float DCH_rin_squared  = 350*350;
			float DCH_rout_squared = 2000*2000;
			IsSecondaryWithinDCH = (vertexZabs < DCH_halflengh) && (vertexRsquared>DCH_rin_squared) && (vertexRsquared<DCH_rout_squared);
		}
		Momentum = sqrt((input_sim_hit.getMomentum().x * input_sim_hit.getMomentum().x) + (input_sim_hit.getMomentum().y * input_sim_hit.getMomentum().y) + (input_sim_hit.getMomentum().z * input_sim_hit.getMomentum().z)) ;

		NCld =0;
		int pdgid = thisParticle.getPDG();
		if( IsSecondaryWithinDCH && 11 == pdgid){
			NCld++;
			Ekdelta=(TMath::Sqrt(Momentum*Momentum+me*me)-me)*1e6;
			vecExtraD.push_back(Ekdelta);
		}
		debug() << "NCld= "<<NCld<<"vecExtraD= "<<vecExtraD<<endmsg;
		if (NCld >0)    debug() << "There is delta rays"<<endmsg;

		mass=(thisParticle.getMass()/1000.);// mass in GeV, required in MeV
		bg = Momentum/ mass;

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
		// lepton PDG id goes from 11 to 16 (antiparticles have negative id)
		bool IsLepton = (11 <= abs(pdgid) ) && (16 >= abs(pdgid) );
		if( IsLepton ){
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
		 * fcnvl[0] = new TF1("fcnvl-0","[0]*TMath::Landau(x,[1],[2],true)+(1.0-[0])*TMath::Landau(x,[3],[4],true)");
		 * fcnvl[0]->SetRange(0,maxEcut);
		 * fcnvl[0]->SetParameters(vecpar[0],vecpar[1],vecpar[2],vecpar[3],vecpar[4]);
		 * cnvl[0]=nullptr;*/


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
		if(IsLepton){
			ExSgmlen=ExSgmlep*TMath::Sqrt(LengthTrack);
		}else{
			ExSgmlen=ExSgmhad*TMath::Sqrt(LengthTrack);
		}
		maxExECl=(Eloss-maxEx0len+this->myRandom.Gaus(0,ExSgmlen))/maxExSlp;
		//  cout<<"maxExECl=   "<<maxExECl<<endl;
		//	SgmaxExECl=ExSgm/maxExSlp;
		if(maxExECl<EIzs) {maxExECl=0.0;} //EIzs const = 25.6

		//else if(choice==2){
		if(choice==2){
			//	cout<<"start new"<<"Eloss= "<<Eloss<<"EIzp=  "<<EIzp<< "EIzs= "<< EIzs<<"maxExECl=  "<<maxExECl<<"totExECl= "<<totExECl<<endl;
			debug() << "Eloss= "<<Eloss<<"EIzs= "<<EIzs<<"EIzp= "<<EIzp<<"maxExECl= "<<maxExECl<<"totExECl"<<totExECl<<endmsg;

			// avoid possibility of infinity while loop
			int while1counter=0;
			while(Eloss>(EIzp+EIzs)&&maxExECl>totExECl && while1counter<1e6){
				while1counter++;
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
							tmpCorr=this->myRandom.Gaus(CorrpMean[i],CorrpSgm[i]);
						}
					}
					ClSz=TMath::Nint(ExECl*CorrSlp+CorrInt-tmpCorr);
					if(ClSz<2) {ClSz=2;}
					// TODO Alvaro: check what is ClSz???
					// output_digi_hit.setClusterSize(ClSz);
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

		// avoid possibility of infinity while loop
		int while2counter=0;
		while(Eloss>=EIzp && while2counter<1e6){
			while2counter++;

			Eloss-=EIzp;
			ExECl1=exGauss->GetRandom();
			//	hExtraCl1->Fill(ExECl1);
			if(ExECl1>Eloss) {ExECl1=Eloss;}
			ExECl1totRec+=ExECl1;
			NCl1++;
			//				EIzp=hEIzp->GetRandom();
			Eloss-=ExECl1;
			ClSz=1;
			// TODO Alvaro: check what is ClSz???
			// output_digi_hit.setClusterSize(ClSz);
			//	hClSzRec->Fill(ClSz);
			//	hClSzRecND->Fill(ClSz);
			//   fClusterCharge.push_back(ClSz);
			NEltot+=ClSz;

		}

		float tmpCl;
		for(unsigned int i=0;i<vecExtraD.size();++i){
			int tmphE=(vecExtraD[i]-minE)/binE;
			if(tmphE>=nhE) tmphE=nhE-1;
			//			cout<<" tmphE "<<tmphE<<endl;
			if(tmphE==nhE-1){
				rndCorr=this->myRandom.Uniform(0,1);
				if(rndCorr<Corrdglfrac[tmphE]){
					tmpCl=this->myRandom.Gaus(Corrdglmeang[tmphE],Corrdglsgmg[tmphE]);
				}else{
					tmpCl=this->myRandom.Landau(Corrdglmpvl[tmphE],Corrdglsgml[tmphE]);
				}
			}else{
				tmpCl=this->myRandom.Gaus(Corrdgmean[tmphE],Corrdgsgm[tmphE]);
			}

			ClSz=TMath::Nint(vecExtraD[i]*CorrSlp+CorrInt-tmpCl);
			if(ClSz<2) {ClSz=2;}
			// TODO Alvaro: check what is ClSz???
			// output_digi_hit.setClusterSize(ClSz);
			//		hClSzRec->Fill(ClSz);
			//		hClSzDRec->Fill(ClSz);
			//      fClusterCharge.push_back(ClSz);
			NEltot+=ClSz;

		}

		Ncl=NCl1+NClp+NCld;
		debug() <<"Ncl= "<< Ncl<<" NCl1= "<<NCl1<<"NClp= "<<NClp<<"NCld= "<<NCld<<endmsg;
		debug() <<"NCld= "<< NCld<<" vecExtraD.size() "<<vecExtraD.size()<< endmsg;
		int vecExtraD_size = vecExtraD.size();
		if (NCld != vecExtraD_size) debug() <<"There is a bug"<< endmsg;
		//debug() << "Output Digi Hit collection size: " << output_digi_hits->size() << endmsg;

		// output_digi_hit.setClusterCount(Ncl);

		/* for (int icl=0;icl<Ncl;icl++) {
		 * output_digi_hit.setClusterSize(NEltot);
	}*/
		// output_digi_hit.setClusterSize(NEltot);
    return return_values;
}
