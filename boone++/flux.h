// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005  

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/
#include <vector>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <algorithm>
#include "basefunc.h"



//###################################################################################
// 		Begining of actual functions for flux calculations 
//##################################################################################

double SWpionflux (double Ppar) {
	// Self explanatory

	double p = Ppar+PFUDGE;	
	
	double c1 = 220.7;
	double c2 = 1.080;
	double c3 = 1.000;
	double c4 = 1.978;
	double c5 = 1.32;
	double c6 = 5.572;
	double c7 = 0.0868;
	double c8 = 9.686;
	double protonEN = 8.9;


	double SWang = SWNORM*c1*pow(p,c2)*(1 - p/(protonEN - 1))*exp(-c3*pow(p,c4)/pow(protonEN,c5) - c6*(TFUDGE)*(p - c7*protonEN*pow(cos(TFUDGE),c8)) );
	if (SWang<0){SWang=0;}
	return (SWang);
}



double EK(double p) {
	// Think this is just Energy of Kaon given p. Needed for peters function below.
	return sqrt(p*p+MKao*MKao);
}
double Kflux_FeynmanScaling(double pK, double th) {
	double beta = 8.89/(0.938+pow(8.89*8.89+0.938*0.938,0.5));
	double gamma = pow(1-beta*beta,-0.5);

	double p = pK + 1.8; //the offset should be fitted.
	double alpha = 2.0; // this should be fitted.

	double C1 = 11.7;
	double C2 = 0.88;
	double C3 = 4.77;
	double C4 = 1.51;
	double C5 = 2.21;
	double C6 = 2.17;
	double C7 = 1.51;

	double x = gamma*(p*cos(th) - beta*EK(p))/alpha;
	double f = C2*p*sin(th)+C3*pow(fabs(x),C4)+C5*pow(fabs(p*sin(th)),2.0)+C7*pow(fabs(x*p*sin(th)),C6);
	double flux = 2*M_PI*(pow(p,2.0)/EK(p))*C1*(1-fabs(x))*exp(-f);

	double eps = 1E-7;
	if((p<=eps)||(flux<eps)){ return 0.0; }
	else { return 1.5e-11*0.022*flux; }
}

double Kflux_Asaka (double pK) { 	
	//Testing only

	double p0 = 2.1, a0 = 2.4e-13, al = 1.285, bl = 0.6289, ah = 793.5,  bh = 3.230;
	double ans =  a0*(al*pow(pK,bl)*ah*pow(pK + p0,-bh ) )/(al*pow(pK,bl) + ah*pow(pK + p0,-bh) );
	if (pK>10){ans=0;}
	return ans; 
}

 


std::vector<double> kinconregions (double upbnd,  double lowbnd , double Es, double Ms, double Mmes) {

	double low1 =(-sqrt(Es*Es-Ms*Ms)*(-4.0*Mdau*Mdau+4.0*Mmes*Mmes+4.0*Ms*Ms)*lowbnd -8.0*Es*sqrt(0.25*pow(Mdau,4)-Es*Es*Mmes*Mmes-0.5*Mdau*Mdau*Mmes*Mmes+0.25*pow(Mmes,4)-0.5*Mdau*Mdau*Ms*Ms+0.5*Mmes*Mmes*Ms*Ms+0.25*pow(Ms,4)+Es*Es*Mmes*Mmes*lowbnd*lowbnd-Mmes*Mmes*Ms*Ms*lowbnd*lowbnd) )/  (-8.0*Es*Es + 8.0* Es*Es*lowbnd*lowbnd - 8.0*Ms*Ms*lowbnd)  ;
	
	double low2 = (-sqrt(Es*Es-Ms*Ms)*(-4.0*Mdau*Mdau+4.0*Mmes*Mmes+4.0*Ms*Ms)*lowbnd +8.0*Es*sqrt(0.25*pow(Mdau,4)-Es*Es*Mmes*Mmes-0.5*Mdau*Mdau*Mmes*Mmes+0.25*pow(Mmes,4)-0.5*Mdau*Mdau*Ms*Ms+0.5*Mmes*Mmes*Ms*Ms+0.25*pow(Ms,4)+Es*Es*Mmes*Mmes*lowbnd*lowbnd-Mmes*Mmes*Ms*Ms*lowbnd*lowbnd) )/  (-8.0*Es*Es + 8.0* Es*Es*lowbnd*lowbnd - 8.0*Ms*Ms*lowbnd)  ;

	double up1 = (-sqrt(Es*Es-Ms*Ms)*(-4.0*Mdau*Mdau+4.0*Mmes*Mmes+4.0*Ms*Ms)*upbnd -8.0*Es*sqrt(0.25*pow(Mdau,4)-Es*Es*Mmes*Mmes-0.5*Mdau*Mdau*Mmes*Mmes+0.25*pow(Mmes,4)-0.5*Mdau*Mdau*Ms*Ms+0.5*Mmes*Mmes*Ms*Ms+0.25*pow(Ms,4)+Es*Es*Mmes*Mmes*upbnd*upbnd-Mmes*Mmes*Ms*Ms*upbnd*upbnd) )/  (-8.0*Es*Es + 8.0* Es*Es*upbnd*upbnd - 8.0*Ms*Ms*upbnd)  ;
	
	double up2 =  (-sqrt(Es*Es-Ms*Ms)*(-4.0*Mdau*Mdau+4.0*Mmes*Mmes+4.0*Ms*Ms)*upbnd +8.0*Es*sqrt(0.25*pow(Mdau,4)-Es*Es*Mmes*Mmes-0.5*Mdau*Mdau*Mmes*Mmes+0.25*pow(Mmes,4)-0.5*Mdau*Mdau*Ms*Ms+0.5*Mmes*Mmes*Ms*Ms+0.25*pow(Ms,4)+Es*Es*Mmes*Mmes*upbnd*upbnd-Mmes*Mmes*Ms*Ms*upbnd*upbnd) )/  (-8.0*Es*Es + 8.0* Es*Es*upbnd*upbnd - 8.0*Ms*Ms*upbnd)  ;


	std::vector<double> result;
	result.push_back(low1);
	result.push_back(low2);
	result.push_back(up1);	
	result.push_back(up2);		
return result;
}


// Structures for integration regimes below
struct fluxparams {double Mmes; double Es; double Ms; double Len;};
struct Tfluxparams {double Mmes; double Es; double Ms;};



//cosint (Ie, i have integrated cos, (via dirac) and its just a function of Ppar and Len)
double cosint (double Ppar, void * p) {
	 struct fluxparams * params  = (struct fluxparams *)p;
	 double Mmes = (params->Mmes);
  	 double Ms = (params->Ms);
 	 double Es = (params->Es);
	 double length = (params->Len);
	 double cosint = 0.0;
	if(Mmes==MP){
		cosint = SWpionflux(Ppar)*FPI*FPI*VUD*VUD*Gf*Gf*pow(Mmes,4)*(pow(Ms/Mmes,2)+pow(Mdau/Mmes,2)-pow(pow(Ms/Mmes,2)-pow(Mdau/Mmes,2),2))*Mmes*exp(-length*Mmes*Piontotal/Ppar)    /(Ppar*Ppar*sqrt(Ppar*Ppar+Mmes*Mmes));
	} else if(Mmes==MKao){
		cosint = Kflux_FeynmanScaling(Ppar,0.0)*FKA*FKA*VUS*VUS*Gf*Gf*pow(Mmes,4)*(pow(Ms/Mmes,2)+pow(Mdau/Mmes,2)-pow(pow(Ms/Mmes,2)-pow(Mdau/Mmes,2),2))*Mmes*exp(-length*Mmes*Piontotal/Ppar)    /(Ppar*Ppar*sqrt(Ppar*Ppar+Mmes*Mmes));

	} else {
		std::cout<<"##################### ERROR #####################"<<std::endl;
		std::cout<<"Currently cosint onlt takes pion and kaon"<<std::endl;
		std::cout<<"##################### ERROR #####################"<<std::endl;
		cosint = 0;
	}
	return cosint;
}



//Pint, so now I am integrating over P, with known kinematic bounds, so now lits just a cfunction of LEN
double pint (double len,  void * p) {
	 struct Tfluxparams * params  = (struct Tfluxparams *)p;
	 double Mmes = (params->Mmes);
  	 double Ms = (params->Ms);
 	 double Es = (params->Es);
	if (Es<=Ms){ 
		std::cout<<"##################### ERROR #####################"<<std::endl;
		std::cout<<"In pint, Es<=Ms, Es: "<<Es<<" and Ms: "<<Ms<<std::endl;
		std::cout<<"maybe check ebining?"<<std::endl;
		std::cout<<"##################### ERROR #####################"<<std::endl;
	}

	std::vector<double> intBounds = removeNAN( kinconregions(CosIntBoundLow,CosIntBoundHigh,Es,Ms,Mmes) );
	double numbound = intBounds.size();
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        struct fluxparams params2 = {Mmes, Es, Ms, len};
	double result,resultpart1,resultpart2,error,errorpart1,errorpart2;
	gsl_function F;
        F.function = &cosint;
	F.params = &params2;

	if( numbound==2){
		//std::cout<<intBounds[0]<<"  "<<intBounds[1]<<std::endl;
		gsl_integration_qags (&F,intBounds[0],intBounds[1], 0, 1e-4, 1000, w, &result, &error); 
	} else if (numbound==4) {
		//std::cout<<intBounds[0]<<"  "<<intBounds[1]<<"  "<<intBounds[2]<<"  "<<intBounds[3]<<" es "<<Es<<" Ms "<<Ms<<std::endl;
		gsl_integration_qags (&F,intBounds[0],intBounds[1], 0, 1e-4, 1000, w, &resultpart1, &errorpart1); 
		gsl_integration_qags (&F,intBounds[2],intBounds[3], 0, 1e-4, 1000, w, &resultpart2, &errorpart2); 
		result=resultpart1+resultpart2;
	} else {
		std::cout<<"##################### ERROR #####################"<<std::endl;
		std::cout<<"integration bounds do not equal 2 or 4.  size of bounds: "<<numbound<<std::endl;
		std::cout<<"##################### ERROR #####################"<<std::endl;
	}

	gsl_integration_workspace_free (w);
	return result;
}

// Lint, so im integrating over L, now its just a nice simple function (no struct) of just M and Es
// Lint is the end resuts for a flux of steriles Ms and Es

double lint (double Ms, double Es, double Mmes) {

	if(Es<=Ms){	//Just incase I pass too many!
		return 0.0;
	}
	if (Ms>MP-Mdau && Mmes==MP){		//Again, pions should be safe this way, must class this up.
		return 0.0;
	}

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        struct Tfluxparams params2 = {Mmes, Es, Ms};
	double result,error;
	gsl_function F;
        F.function = &pint;
	F.params = &params2;
	
	gsl_integration_qags (&F, 0.0, DECAYTUNNEL ,0, 1e-4, 1000, w, &result, &error); 
	gsl_integration_workspace_free (w);

	return result;
}





