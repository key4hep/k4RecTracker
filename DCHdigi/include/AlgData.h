#include <iostream>

#include "TTree.h"
#include "TChain.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include <TSpline.h>
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

using namespace std;


class AlgData{
public:
	AlgData();
	~AlgData();

public:
	void read_file(TString dataAlg);
	TF1* read_graph(TString dataAlg, TString cvName, TString fitName);


	TF1* get_fit();
	TFormula* get_formula();
	double get_Fitvalue(double betagamma)  const;
	ROOT::Math::Interpolator* Interpvalue(vector<double> bg, vector<double> Ydata,int Intertype);
	void Load_file(TString dataAlg);
	void Load_interp();
	double get_MPVExtra (double betagamma)   const;
	double get_SgmExtra (double betagamma)   const;
	double get_MeanExtra1 (double betagamma) const;
	double get_SgmExtra1 (double betagamma)  const;
	double get_SlopeExtra1 (double betagamma)const;
	double get_FracExtra1(double betagamma) const;
	double get_FfracExtra(double betagamma) ;
	double get_maxEx0(double betagamma);
	double get_maxExSlp();
	double get_ExSgmlep();
	double get_ExSgmhad();
	double get_Ffrac(double betagamma);
	double get_Fmpv1(double betagamma);
	double get_Fsgm1(double betagamma);
	double get_Fmpv2(double betagamma);
	double get_Fsgm2(double betagamma);
	double get_ClSzCorrInt(double betagamma);
	double get_ClSzCorrSlp(double betagamma);
	vector<double> get_ClSzCorrpmean(double betagamma);
	vector<double> get_ClSzCorrpsgm(double betagamma);

	vector<double> get_ClSzCorrdgmean(double betagamma);
	vector<double> get_ClSzCorrdgsgm(double betagamma);

	vector<double> get_ClSzCorrdglfrac(double betagamma);
	vector<double> get_ClSzCorrdglmpvl(double betagamma);
	vector<double> get_ClSzCorrdglsgml(double betagamma);
	vector<double> get_ClSzCorrdglmeang(double betagamma);
	vector<double> get_ClSzCorrdglsgmg(double betagamma);


private:
	TString DataAlg;
	TChain *data1;
	TChain *datalep;
	TChain *datahad;
	TFile *file;

	double bgT,bgTlep,bgThad,tmpmaxEx0Tot,tmpErrmaxEx0Tot,tmpmaxExSlpTot,tmpErrmaxExSlpTot,tmpExSgmTotlep,tmpExSgmTothad;
	double tmptFfracTot,tmptFerrfracTot,tmptFmpv1Tot,tmptFerrmpv1Tot,tmptFsgm1Tot,tmptFerrsgm1Tot,tmptFmpv2Tot,tmptFerrmpv2Tot,tmptFsgm2Tot,tmptFerrsgm2Tot;
	double tmpCorrIntTot,tmpErrCorrIntTot,tmpCorrSlpTot,tmpErrCorrSlpTot;

	vector <double> *tmpCorrmeanSliceTot=nullptr;
	vector <double> *tmpErrCorrmeanSliceTot=nullptr;
	vector <double> *tmpCorrsgmSliceTot=nullptr;
	vector <double> *tmpErrCorrsgmSliceTot=nullptr;

	vector <double> *tmpCorrdglfracSliceTot=nullptr;
	vector <double> *tmpErrCorrdglfracSliceTot=nullptr;
	vector <double> *tmpCorrdglmeangSliceTot=nullptr;
	vector <double> *tmpErrCorrdglmeangSliceTot=nullptr;
	vector <double> *tmpCorrdglsgmgSliceTot=nullptr;
	vector <double> *tmpErrCorrdglsgmgSliceTot=nullptr;
	vector <double> *tmpCorrdglmpvSliceTot=nullptr;
	vector <double> *tmpErrCorrdglmpvSliceTot=nullptr;
	vector <double> *tmpCorrdglsgmlSliceTot=nullptr;
	vector <double> *tmpErrCorrdglsgmlSliceTot=nullptr;

	vector <double> *tmpCorrdgmeanSliceTot=nullptr;
	vector <double> *tmpErrCorrdgmeanSliceTot=nullptr;
	vector <double> *tmpCorrdgsgmSliceTot=nullptr;
	vector <double> *tmpErrCorrdgsgmSliceTot=nullptr;

	//____________________________________________________//
	vector <double> bgv,bglep,bghad;
	vector <double> maxEx0Tot,ErrmaxEx0Tot,maxExSlpTot,ErrmaxExSlpTot,ExSgmTotlep,ExSgmTothad;
	vector <double>  tFfracTot, tFerrfracTot, tFmpv1Tot, tFerrmpv1Tot, tFsgm1Tot, tFerrsgm1Tot, tFmpv2Tot, tFerrmpv2Tot, tFsgm2Tot, tFerrsgm2Tot;
	vector <double> CorrmeanSliceTot[20],ErrCorrmeanSliceTot[20],CorrsgmSliceTot[20],ErrCorrsgmSliceTot[20];
	vector <double> CorrIntTot,ErrCorrIntTot,CorrSlpTot,ErrCorrSlpTot;
	vector <double> CorrdglfracSliceTot[20],ErrCorrdglfracSliceTot[20],CorrdglmeangSliceTot[20],ErrCorrdglmeangSliceTot[20],CorrdglsgmgSliceTot[20],ErrCorrdglsgmgSliceTot[20],CorrdglmpvSliceTot[20],ErrCorrdglmpvSliceTot[20],CorrdglsgmlSliceTot[20],ErrCorrdglsgmlSliceTot[20];
	vector <double> CorrdgmeanSliceTot[20],ErrCorrdgmeanSliceTot[20],CorrdgsgmSliceTot[20],ErrCorrdgsgmSliceTot[20];

	//__________________________________________//
	TF1 *Ft;
	TF1 *FitMPVExtra;
	TF1 *FitSgmExtra;
	TF1 *FitMeanExtra1;
	TF1 *FitSgmExtra1;
	TF1 *FitSlopeExtra1;
	TF1 *FitFracExtra1;

	TFormula *FtFormula;

	vector<TString> Corrmean;
	vector<TString> Corrsgm;

	vector<TString> Corrdgmean;
	vector<TString> Corrdgsgm;

	vector<TString> Corrdglfrac;
	vector<TString> Corrdglmeang;
	vector<TString> Corrdglsgmg;
	vector<TString> Corrdglmpv;
	vector<TString> Corrdglsgml;

	vector<double> ClSzCorrpmean;
	vector<double> ClSzCorrpsgm;
	vector<double> ClSzCorrdgmean;
	vector<double> ClSzCorrdgsgm;
	vector<double> ClSzCorrdglfrac;
	vector<double> ClSzCorrdglmpvl;
	vector<double> ClSzCorrdglsgml;
	vector<double> ClSzCorrdglmeang;
	vector<double> ClSzCorrdglsgmg;

	double maxExSlp;
	double ExSgmlep;
	double ExSgmhad;
	void Calc_maxExSlp();
	void Calc_ExSgmlep();
	void Calc_ExSgmhad();



	//___________________________________________//

	map<TString,ROOT::Math::Interpolator*> itpm;
	ROOT::Math::Interpolator* itp1;

};

	//_______________________________________________//

AlgData::AlgData(){

	DataAlg="";
	data1=nullptr;
	datalep=nullptr;
	datahad=nullptr;
	file=nullptr;


	tmpCorrmeanSliceTot=nullptr;
	tmpErrCorrmeanSliceTot=nullptr;
	tmpCorrsgmSliceTot=nullptr;
	tmpErrCorrsgmSliceTot=nullptr;

	tmpCorrdglfracSliceTot=nullptr;
	tmpErrCorrdglfracSliceTot=nullptr;
	tmpCorrdglmeangSliceTot=nullptr;
	tmpErrCorrdglmeangSliceTot=nullptr;
	tmpCorrdglsgmgSliceTot=nullptr;
	tmpErrCorrdglsgmgSliceTot=nullptr;
	tmpCorrdglmpvSliceTot=nullptr;
	tmpErrCorrdglmpvSliceTot=nullptr;
	tmpCorrdglsgmlSliceTot=nullptr;
	tmpErrCorrdglsgmlSliceTot=nullptr;

	tmpCorrdgmeanSliceTot=nullptr;
	tmpErrCorrdgmeanSliceTot=nullptr;
	tmpCorrdgsgmSliceTot=nullptr;
	tmpErrCorrdgsgmSliceTot=nullptr;

	Ft=nullptr;
	FtFormula=nullptr;
	FitMPVExtra=nullptr;
	FitSgmExtra=nullptr;
	FitMeanExtra1=nullptr;
	FitSgmExtra1=nullptr;
	FitSlopeExtra1=nullptr;
	FitFracExtra1=nullptr;


}

	AlgData::~AlgData(){
		if(file->IsOpen())file->Close();
	}


	void AlgData::read_file(TString dataAlg){

		DataAlg=dataAlg;

		cout<<" ---reading of---"<<DataAlg<<endl;

		data1 = new TChain("DataAlg");
		datalep =new TChain("DataAlglep");
		datahad =new TChain("DataAlghad");

		if (dataAlg.Contains(".root")) {

			data1->Add(dataAlg.Data());
			datalep->Add(dataAlg.Data());
			datahad->Add(dataAlg.Data());

		}

		data1->SetBranchAddress("bgT",&bgT);
		data1->SetBranchAddress("tmpmaxEx0Tot",&tmpmaxEx0Tot);
		data1->SetBranchAddress("tmpErrmaxEx0Tot",&tmpErrmaxEx0Tot);
		data1->SetBranchAddress("tmpmaxExSlpTot",&tmpmaxExSlpTot);
		data1->SetBranchAddress("tmpErrmaxExSlpTot",&tmpErrmaxExSlpTot);
		data1->SetBranchAddress("tmptFfracTot",&tmptFfracTot);
		data1->SetBranchAddress("tFerrfracTot",&tmptFerrfracTot);
		data1->SetBranchAddress("tmptFerrfracTot",&tmptFmpv1Tot);
		data1->SetBranchAddress("tmptFerrmpv1Tot",&tmptFerrmpv1Tot);
		data1->SetBranchAddress("tmptFerrmpv1Tot",&tmptFerrmpv1Tot);
		data1->SetBranchAddress("tmptFsgm1Tot",&tmptFsgm1Tot);
		data1->SetBranchAddress("tmptFerrsgm1Tot",&tmptFerrsgm1Tot);
		data1->SetBranchAddress("tmptFmpv2Tot",&tmptFmpv2Tot);
		data1->SetBranchAddress("tmptFerrmpv2Tot",&tmptFerrmpv2Tot);
		data1->SetBranchAddress("tmptFsgm2Tot",&tmptFsgm2Tot);
		data1->SetBranchAddress("tmptFerrsgm2Tot",&tmptFerrsgm2Tot);

		data1->SetBranchAddress("tmpCorrIntTot",&tmpCorrIntTot);
		data1->SetBranchAddress("tmpErrCorrIntTot",&tmpErrCorrIntTot);
		data1->SetBranchAddress("tmpCorrSlpTot",&tmpCorrSlpTot);
		data1->SetBranchAddress("tmpErrCorrSlpTot",&tmpErrCorrSlpTot);

		data1->SetBranchAddress("tmpCorrmeanSliceTot",&tmpCorrmeanSliceTot);
		data1->SetBranchAddress("tmpErrCorrmeanSliceTot",&tmpErrCorrmeanSliceTot);
		data1->SetBranchAddress("tmpCorrsgmSliceTot",&tmpCorrsgmSliceTot);
		data1->SetBranchAddress("tmpErrCorrsgmSliceTot",&tmpErrCorrsgmSliceTot);

		data1->SetBranchAddress("tmpCorrdgmeanSliceTot",&tmpCorrdgmeanSliceTot);
		data1->SetBranchAddress("tmpErrCorrdgmeanSliceTot",&tmpErrCorrdgmeanSliceTot);
		data1->SetBranchAddress("tmpCorrdgsgmSliceTot",&tmpCorrdgsgmSliceTot);
		data1->SetBranchAddress("tmpErrCorrdgsgmSliceTot",&tmpErrCorrdgsgmSliceTot);

		data1->SetBranchAddress("tmpCorrdglfracSliceTot",&tmpCorrdglfracSliceTot);
		data1->SetBranchAddress("tmpErrCorrdglfracSliceTot",&tmpErrCorrdglfracSliceTot);
		data1->SetBranchAddress("tmpCorrdglmeangSliceTot",&tmpCorrdglmeangSliceTot);
		data1->SetBranchAddress("tmpErrCorrdglmeangSliceTot",&tmpErrCorrdglmeangSliceTot);
		data1->SetBranchAddress("tmpCorrdglsgmgSliceTot",&tmpCorrdglsgmgSliceTot);
		data1->SetBranchAddress("tmpErrCorrdglsgmgSliceTot",&tmpErrCorrdglsgmgSliceTot);
		data1->SetBranchAddress("tmpCorrdglmpvSliceTot",&tmpCorrdglmpvSliceTot);
		data1->SetBranchAddress("tmpErrCorrdglmpvSliceTot",&tmpErrCorrdglmpvSliceTot);
		data1->SetBranchAddress("tmpCorrdglsgmlSliceTot",&tmpCorrdglsgmlSliceTot);
		data1->SetBranchAddress("tmpErrCorrdglsgmlSliceTot",&tmpErrCorrdglsgmlSliceTot);

		datalep->SetBranchAddress("bgTlep",&bgTlep);
		datalep->SetBranchAddress("tmpExSgmTotlep",&tmpExSgmTotlep);

		datahad->SetBranchAddress("bgThad",&bgThad);
		datahad->SetBranchAddress("tmpExSgmTothad",&tmpExSgmTothad);


		for (int i=0;i<data1->GetEntries();i++) {
			data1->GetEntry(i);
			bgv.push_back(bgT);
//			cout<<bgT<<endl;
			maxEx0Tot.push_back(tmpmaxEx0Tot);
			ErrmaxEx0Tot.push_back(tmpErrmaxEx0Tot);
			maxExSlpTot.push_back(tmpmaxExSlpTot);
			ErrmaxExSlpTot.push_back(tmpErrmaxExSlpTot);

			tFfracTot.push_back(tmptFfracTot);
			tFerrfracTot.push_back(tmptFerrfracTot);
			tFmpv1Tot.push_back(tmptFmpv1Tot);
			tFerrmpv1Tot.push_back(tmptFerrmpv1Tot);
			tFsgm1Tot.push_back(tmptFsgm1Tot);
			tFerrsgm1Tot.push_back(tmptFerrsgm1Tot);
			tFmpv2Tot.push_back(tmptFmpv2Tot);
			tFerrmpv2Tot.push_back(tmptFerrmpv2Tot);
			tFsgm2Tot.push_back(tmptFsgm2Tot);
			tFerrsgm2Tot.push_back(tmptFerrsgm2Tot);

			CorrIntTot.push_back(tmpCorrIntTot);
			ErrCorrIntTot.push_back(tmpErrCorrIntTot);
			CorrSlpTot.push_back(tmpCorrSlpTot);
			ErrCorrSlpTot.push_back(tmpErrCorrSlpTot);
			for(unsigned int j=0;j<tmpCorrmeanSliceTot->size();++j){
				CorrmeanSliceTot[j].push_back(tmpCorrmeanSliceTot->at(j));
				ErrCorrmeanSliceTot[j].push_back(tmpErrCorrmeanSliceTot->at(j));
				CorrsgmSliceTot[j].push_back(tmpCorrsgmSliceTot->at(j));
				ErrCorrsgmSliceTot[j].push_back(tmpErrCorrsgmSliceTot->at(j));
			}

			for( int n=0;n<tmpCorrdgmeanSliceTot->size();++n){
				CorrdgmeanSliceTot[n].push_back(tmpCorrdgmeanSliceTot->at(n));
				ErrCorrdgmeanSliceTot[n].push_back(tmpErrCorrdgmeanSliceTot->at(n));
				CorrdgsgmSliceTot[n].push_back(tmpCorrdgsgmSliceTot->at(n));
				ErrCorrdgsgmSliceTot[n].push_back(tmpErrCorrdgsgmSliceTot->at(n));
			}

			for(unsigned int k=0;k<tmpCorrdglfracSliceTot->size();++k){
				CorrdglfracSliceTot[k].push_back(tmpCorrdglfracSliceTot->at(k));
				ErrCorrdglfracSliceTot[k].push_back(tmpErrCorrdglfracSliceTot->at(k));
				CorrdglmeangSliceTot[k].push_back(tmpCorrdglmeangSliceTot->at(k));
				ErrCorrdglmeangSliceTot[k].push_back(tmpErrCorrdglmeangSliceTot->at(k));
				CorrdglsgmgSliceTot[k].push_back(tmpCorrdglsgmgSliceTot->at(k));
				ErrCorrdglsgmgSliceTot[k].push_back(tmpErrCorrdglsgmgSliceTot->at(k));
				CorrdglmpvSliceTot[k].push_back(tmpCorrdglmpvSliceTot->at(k));
				ErrCorrdglmpvSliceTot[k].push_back(tmpErrCorrdglmpvSliceTot->at(k));
				CorrdglsgmlSliceTot[k].push_back(tmpCorrdglsgmlSliceTot->at(k));
				ErrCorrdglsgmlSliceTot[k].push_back(tmpErrCorrdglsgmlSliceTot->at(k));
//				cout<<tmpCorrdglfracSliceTot->at(k)<<endl;
			}

		}

		for (int i=0;i<datalep->GetEntries();i++) {
			datalep->GetEntry(i);
			bglep.push_back(bgTlep);
			ExSgmTotlep.push_back(tmpExSgmTotlep);
//			cout<< " bgTlep "<<bgTlep<<" tmpExSgmTotlep " <<tmpExSgmTotlep<<endl;
		}

		for (int i=0;i<datahad->GetEntries();i++) {
			datahad->GetEntry(i);
			bghad.push_back(bgThad);
			ExSgmTothad.push_back(tmpExSgmTothad);
		}

	}



	TF1* AlgData::read_graph(TString dataAlg, TString cvName, TString fitName){
		DataAlg=dataAlg;
		file = TFile::Open(dataAlg.Data(),"read");
		TCanvas *cv=(TCanvas *)file->Get(cvName.Data());
		TGraph *gr=(TGraph *)cv->GetListOfPrimitives()->FindObject("Graph");
		TF1* ft=(TF1*)gr->GetListOfFunctions()->FindObject(fitName.Data());
		TFormula *ftFormula=(TFormula *)ft->GetFormula();
		FtFormula=ftFormula;
		return ft;
	}

	void AlgData::Load_file(TString dataAlg){
		read_file(dataAlg.Data());
		FitMPVExtra=read_graph(dataAlg.Data(), "cMPVExtrabgTot", "fit_MPV");
		FitSgmExtra=read_graph(dataAlg.Data(), "cSgmExtrabgTot", "fit_sgmEx");
		FitMeanExtra1=read_graph(dataAlg.Data(), "cMeanExtra1bgTot", "expeff");
		FitSgmExtra1=read_graph(dataAlg.Data(), "cSgmExtra1bgTot", "expeffNeg");
		FitSlopeExtra1=read_graph(dataAlg.Data(), "cSlopeExtra1bgTot", "fit_slp");
		FitFracExtra1=read_graph(dataAlg.Data(), "cfracbgTot", "fit_frEx");
		Calc_ExSgmhad();
		Calc_ExSgmlep();
		Calc_maxExSlp();
	}

	TF1* AlgData::get_fit(){
		return Ft;
	}

	TFormula * AlgData::get_formula(){
		return FtFormula;
	}

	double AlgData::get_MPVExtra (double betagamma)  const{
		return FitMPVExtra->Eval(betagamma);
	}

	double AlgData::get_SgmExtra(double betagamma) const{
		return FitSgmExtra->Eval(betagamma);
	}

	double AlgData::get_MeanExtra1(double betagamma) const{
		return FitMeanExtra1->Eval(betagamma);
	}

	double AlgData::get_SgmExtra1(double betagamma) const{
		return FitSgmExtra1->Eval(betagamma);
	}

	double AlgData::get_FracExtra1(double betagamma) const{
		return FitFracExtra1->Eval(betagamma);
	}

	double AlgData::get_SlopeExtra1(double betagamma) const{
		return FitSlopeExtra1->Eval(betagamma);
	}


	ROOT::Math::Interpolator* AlgData::Interpvalue(vector<double>bg,vector<double>Ydata,int Intertype) {
		//tipo di interpolazione 0=linear;1=pol;2=cspline;4=akima
		return itp1=new ROOT::Math::Interpolator(bg,Ydata,(ROOT::Math::Interpolation::Type)Intertype);
	}

	void  AlgData::Load_interp(){
		itpm["maxEx0Tot"]=Interpvalue(bgv,maxEx0Tot,4);
		itpm["tFfracTot"]=Interpvalue(bgv,tFfracTot,4);
		itpm["tFmpv1Tot"]=Interpvalue(bgv,tFmpv1Tot,4);
		itpm["tFsgm1Tot"]=Interpvalue(bgv,tFsgm1Tot,4);
		itpm["tFmpv2Tot"]=Interpvalue(bgv,tFmpv2Tot,4);
		itpm["tFsgm2Tot"]=Interpvalue(bgv,tFsgm2Tot,4);

		itpm["CorrIntTot"]=Interpvalue(bgv,CorrIntTot,4);
		itpm["CorrSlpTot"]=Interpvalue(bgv,CorrSlpTot,4);
		for(unsigned int i=0;i<tmpCorrmeanSliceTot->size();++i){
			TString nameCm="CorrmeanSliceTot";
			TString nameCs="CorrsgmSliceTot";
			Corrmean.push_back(nameCm+i);
			Corrsgm.push_back(nameCs+i);
			itpm[Corrmean.at(i)]=Interpvalue(bgv, CorrmeanSliceTot[i], 4);
			itpm[Corrsgm.at(i)]=Interpvalue(bgv, CorrsgmSliceTot[i], 4);

		}

		for(unsigned int j=0;j<tmpCorrdglfracSliceTot->size();++j){
			TString nameCdf="Corrdglfrac";
			TString nameCdm="Corrdglmeang";
			TString nameCds="Corrdglsgmg";
			TString nameCdmpv="Corrdglmpv";
			TString nameCdsgl="Corrdglsgmgl";
			Corrdglfrac.push_back(nameCdf+j);
			Corrdglmeang.push_back(nameCdm+j);
			Corrdglsgmg.push_back(nameCds+j);
			Corrdglmpv.push_back(nameCdmpv+j);
			Corrdglsgml.push_back(nameCdsgl+j);
			itpm[Corrdglfrac.at(j)]=Interpvalue(bgv, CorrdglfracSliceTot[j], 4);
			itpm[Corrdglmeang.at(j)]=Interpvalue(bgv, CorrdglmeangSliceTot[j], 4);
			itpm[Corrdglsgmg.at(j)]=Interpvalue(bgv, CorrdglsgmgSliceTot[j], 4);
			itpm[Corrdglmpv.at(j)]=Interpvalue(bgv, CorrdglmpvSliceTot[j], 4);
			itpm[Corrdglsgml.at(j)]=Interpvalue(bgv, CorrdglsgmlSliceTot[j], 4);
		}

		for(unsigned int n=0;n<tmpCorrdgmeanSliceTot->size();++n){
			TString nameCdmg="Corrdgmean";
			TString nameCdsgg="Corrdgsgm";
			Corrdgmean.push_back(nameCdmg+n);
			Corrdgsgm.push_back(nameCdsgg+n);
			itpm[Corrdgmean.at(n)]=Interpvalue(bgv, CorrdgmeanSliceTot[n], 4);
			itpm[Corrdgsgm.at(n)]=Interpvalue(bgv, CorrdgsgmSliceTot[n], 4);
		}
	}

	double AlgData::get_maxEx0(double betagamma) {
		return itpm["maxEx0Tot"]->Eval(betagamma);
	}
	double AlgData::get_Ffrac(double betagamma) {
		return itpm["tFfracTot"]->Eval(betagamma);
	}
	double AlgData::get_Fmpv1(double betagamma) {
		return itpm["tFmpv1Tot"]->Eval(betagamma);
	}
	double AlgData::get_Fsgm1(double betagamma) {
		return itpm["tFsgm1Tot"]->Eval(betagamma);
	}
	double AlgData::get_Fmpv2(double betagamma) {
		return itpm["tFmpv2Tot"]->Eval(betagamma);
	}
	double AlgData::get_Fsgm2(double betagamma) {
		return itpm["tFsgm2Tot"]->Eval(betagamma);
	}
	double AlgData::get_ClSzCorrInt(double betagamma){
		return itpm["CorrIntTot"]->Eval(betagamma);
	}
	double AlgData::get_ClSzCorrSlp(double betagamma){
		return itpm["CorrSlpTot"]->Eval(betagamma);
	}

	vector<double> AlgData::get_ClSzCorrpmean(double betagamma){
		ClSzCorrpmean.clear();
//		cout<<"size "<<tmpCorrmeanSliceTot->size()<<endl;
		for(unsigned int i=0;i<tmpCorrmeanSliceTot->size();++i){
			ClSzCorrpmean.push_back(itpm[Corrmean.at(i)]->Eval(betagamma));
		}

		return ClSzCorrpmean;
	}

	vector<double> AlgData::get_ClSzCorrpsgm(double betagamma){
		ClSzCorrpsgm.clear();
		for(unsigned int i=0;i<tmpCorrmeanSliceTot->size();++i){
			ClSzCorrpsgm.push_back(itpm[Corrsgm.at(i)]->Eval(betagamma));
		}
		return ClSzCorrpsgm;
	}

	vector<double> AlgData::get_ClSzCorrdglfrac(double betagamma){
		ClSzCorrdglfrac.clear();
		for(unsigned int i=0;i<tmpCorrdglfracSliceTot->size();++i){
			ClSzCorrdglfrac.push_back(itpm[Corrdglfrac.at(i)]->Eval(betagamma));
		}
		return ClSzCorrdglfrac;
	}

	vector<double> AlgData::get_ClSzCorrdglmeang(double betagamma){
		ClSzCorrdglmeang.clear();
		for(unsigned int i=0;i<tmpCorrdglfracSliceTot->size();++i){
			ClSzCorrdglmeang.push_back(itpm[Corrdglmeang.at(i)]->Eval(betagamma));
		}
		return ClSzCorrdglmeang;
	}

	vector<double> AlgData::get_ClSzCorrdglsgmg(double betagamma){
		ClSzCorrdglsgmg.clear();
		for(unsigned int i=0;i<tmpCorrdglfracSliceTot->size();++i){
			ClSzCorrdglsgmg.push_back(itpm[Corrdglsgmg.at(i)]->Eval(betagamma));
		}
		return ClSzCorrdglsgmg;
	}

	vector<double> AlgData::get_ClSzCorrdglmpvl(double betagamma){
		ClSzCorrdglmpvl.clear();
		for(unsigned int i=0;i<tmpCorrdglfracSliceTot->size();++i){
			ClSzCorrdglmpvl.push_back(itpm[Corrdglmpv.at(i)]->Eval(betagamma));
		}
		return ClSzCorrdglmpvl;
	}

	vector<double> AlgData::get_ClSzCorrdglsgml(double betagamma){
		ClSzCorrdglsgml.clear();
		for(unsigned int i=0;i<tmpCorrdglfracSliceTot->size();++i){
			ClSzCorrdglsgml.push_back(itpm[Corrdglsgml.at(i)]->Eval(betagamma));
		}
		return ClSzCorrdglsgml;
	}

	vector<double> AlgData::get_ClSzCorrdgmean(double betagamma){
		ClSzCorrdgmean.clear();
		for(unsigned int i=0;i<tmpCorrdgmeanSliceTot->size();++i){
			ClSzCorrdgmean.push_back(itpm[Corrdgmean.at(i)]->Eval(betagamma));
		}
		return ClSzCorrdgmean;
	}

	vector<double> AlgData::get_ClSzCorrdgsgm(double betagamma){
		ClSzCorrdgsgm.clear();
		for(unsigned int i=0;i<tmpCorrdgmeanSliceTot->size();++i){
			ClSzCorrdgsgm.push_back(itpm[Corrdgsgm.at(i)]->Eval(betagamma));
		}
		return ClSzCorrdgsgm;
	}


	double AlgData::get_maxExSlp(){
		return maxExSlp;
	}

	void AlgData::Calc_maxExSlp(){
		double sum=0.0;
		for(unsigned int i=0;i<maxExSlpTot.size();++i){
			sum+=maxExSlpTot[i];
		}
		 maxExSlp=sum/(double)maxExSlpTot.size();
	}


	double AlgData::get_ExSgmlep(){
		return ExSgmlep;
	}

	void AlgData::Calc_ExSgmlep(){
		double sum=0.0;
		for(unsigned int i=0;i<ExSgmTotlep.size();++i){
			sum+=ExSgmTotlep[i];
		}
		ExSgmlep=sum/(double)ExSgmTotlep.size();
	}

	double AlgData::get_ExSgmhad(){
		return ExSgmhad;
	}

	void AlgData::Calc_ExSgmhad(){
			double sum=0.0;
			for(unsigned int i=0;i<ExSgmTothad.size();++i){
				sum+=ExSgmTothad[i];
			}
			ExSgmhad=sum/(double)ExSgmTothad.size();
		}

