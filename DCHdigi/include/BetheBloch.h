#pragma once
#include "TMath.h"
#include <cmath>

namespace BB
{
//Atomic number of material
double Z=2;

//Relative atomic mass
double A=4.0026;

//Atomic number of particle
double z_particle=1;

//Electron rest mass MeV
double me = 0.511;

//Electron charge in eV
double e = 1;

double K = 0.307075;

//double I=41.8*1e-6;
//variables for delta
double c=11.1393;
double x0=2.2017;
double x1=3.6122;
double a=0.13443;
double k1=5.8347;
const double k2=2.0*TMath::Log(10);

//double getDelta(double p, double m){
//	double delta;
//	double x= TMath::Log10(p/m);
//
//	cout<<x<<"\t"<<p<<endl;
//	if(x>=x1) delta=k2*x-c;
//	if(x<x1&&x>=x0) delta=k2*x-c+a*pow((x1-x),k1);
//	if(x<x0) delta=0.0;
//
//	return delta;
//}

double bethe_bloch( double *xp,double *par){
	double m=par[0];
	double I=par[1];
	double p=xp[0]*1e+3;
	double x = TMath::Log10(p/m);
	double delta = 0.0;

	double beta=p/TMath::Sqrt(p*p+m*m);
	double gam=1/TMath::Sqrt(1-beta*beta);

	if(x>=x1) delta=k2*x-c;
	if(x<x1&&x>=x0) delta=k2*x-c+a*pow((x1-x),k1);
	if(x<x0) delta=0.0;

	//Tmax formula
	double Tmax = 2.0*me*pow(gam*beta,2)/(1+(2.0*gam*me/m)+pow(me/m,2));
	//Bethe-bloch formula
	double func = K*pow(z_particle*e,2)*(Z/A)/pow(beta,2)*(0.5*TMath::Log(2.0*me*pow(beta*gam,2)*Tmax/(I*I))-pow(beta,2)-0.5*delta);
	return func;
}
}
