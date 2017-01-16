#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <globes/globes.h>
#ifndef  BASEFNC_H  
#include "basefunc.h"
#endif


#define INITIAL 1;
#define FINAL 1;
#define PANTI 1;


extern double theta12 ; //asin(sqrt(0.311));
extern 	double theta13 ; //asin(sqrt(0.0255));
extern	double theta23 ;
extern	double sdm ;
extern	double ldm_N;
extern	double ldm_I;
extern	double deltacp ;
extern double eBinStep;
extern double eBinMax;
extern double eBinMin;
extern std::vector<double > ebin;


//void changeGlobesFluxFile (std::string filename){
//	std::string s1 = "sed -ie 's/'tester.dat'/" , s2 = filename ,s3 ="/' temp.glb" , result;
//	result = s1+s2+s3;
//	system("cp template.glb temp.glb");
//	system(result.c_str());
//}

double diffdecay_transit (double Ms, double Es, double Ev, double U, double Gdec){ 
	// Checked with "Workign"  mathematica script.

	double result = 999999999;
	double G = gamma(Ms,Es);
	double B = beta(Ms,Es);
	
	double bplus = Ms/(2*G*(1+B));
	double bminus = Ms/(2*G*(1-B));

	double a = 4.0/3.0; // alpha is 4/3 for nu and e, 2 for nuebar. 

	if (Ev > 0 && Ev < bplus ){
	//	std::cout<<"case1: "<<bplus<<"  "<<bminus<<"  "<<Ev<<std::endl;
		result = pow(Gf*U*Ms*Ev,2)*G/(2*pow(PI,3))*(	1+Ev*a*G*(1+B*B/3)/Ms		);
	} else if (Ev >bplus && Ev < bminus ) {

		result = pow(Gf*U,2)*pow(Ms,4)/(8*pow(PI,3)*G*B)* (      3*a/4*pow(Ev*G*(1-B)/Ms, 3 )-pow(Ev*G*(1-B)/Ms, 2 )+  (3-a)/12		);
	} else {

		result = 0.0;
	}

	return result*Gdec/Gvvv(Ms,U);
}


double transitspectrum_withE (double Ms, double Ev, double U, double Gdec, const std::vector<std::vector<double > > &fluxptr) {
	double result = 0.0;
	double dectot = decaytotal(Ms,U,Gdec); 

	for (int ie = whichbinM(Ms)+1; ie<ebin.size(); ie++) {

	if(ebin[ie]<Ms) {std::cout<<"TRANSITSPECTRUM ERROR Es<Ms"<<std::endl;}

		double flux_now = fluxptr[(int)floor(Ms/0.001)-1][ie];
		double Es = ebin[ie]+0.025;
		double G = gamma(Ms,Es);
		double B = beta(Ms,Es);

		
		double expsinh_inter = 0.5*dectot*M2G(540)/(G*B);
		if (expsinh_inter > 10){
			result =result+ U*U*flux_now*2*diffdecay_transit(Ms,Es,Ev,U,Gdec)*0.5*eBinStep/dectot;
		} else {
			result=result+ U*U*flux_now*2*diffdecay_transit(Ms,Es,Ev,U,Gdec)*exp(-expsinh_inter)*sinh(expsinh_inter)*eBinStep/dectot;
		}
	//	std::cout<<" Addition "<<addition<<" diffcom "<<diffdecay_transit(Ms,Es,Ev,U,Gdec)<<" dectot "<<dectot<<" flux_now "<<flux_now<<" exp "<<exp(-0.5*dectot*M2G(540)/(G*B))<<" sinh "<<sinh(0.5*dectot*M2G(540)/(G*B) )<<" G "<<G<<" B "<<B<<" Ms "<<Ms<<" Es "<<Es<<" Ev "<<Ev<<" inside sihn "<<0.5*dectot*M2G(540)/(G*B)<<std::endl;

	}
	
	return result;
}

std::vector<double > transitspectrum (double Ms, double U, double Gdec, const std::vector<std::vector<double > > &fluxptr) {
	std::vector<double > result;
	for (int ee=0; ee< ebin.size(); ee++) {
		result.push_back(transitspectrum_withE(Ms,ebin[ee],U,Gdec, fluxptr)  );
	}
return result;
}


//###########################################################################
//##########		Globes related from flux 	   ##################
//###########################################################################
void flux2GlobesFile (std::vector<double > myflux,int whichneut) {
	std::ofstream teststream;  
	teststream.open("temp_nu_flux_file.dat");

	double nue=0.0, nuebar=0.0,numubar=0.0,nutau=0.0,nutaubar=0.0,numu = 0.0;

	// Globes flux file. 501 lines with columns      Energy / Phi_nue / Phi_numu / Phi_nutau / Phi_nuebar / Phi_numubar / Phi_nutaubar
	if ( whichneut == 0) {
		for (int i = 0; i<501; i++){
			if (i<myflux.size()){
				teststream<<i*eBinStep+eBinStep<<"  "<<myflux[i]<<"  "<<numu<<"  "<<nutau<<"  "<<nuebar<<"  "<<numubar<<"  "<<nutaubar<<std::endl;
			} else {
				teststream<<i*eBinStep+eBinStep<<"  "<<nue<<"  "<<numu<<"  "<<nutau<<"  "<<nuebar<<"  "<<numubar<<"  "<<nutaubar<<std::endl;
			}
		}
	} else {
		for (int i = 0; i<501; i++){
			if (i<myflux.size()){
				teststream<<i*eBinStep+eBinStep<<"  "<<nue<<"  "<<myflux[i]<<"  "<<nutau<<"  "<<nuebar<<"  "<<numubar<<"  "<<nutaubar<<std::endl;
			} else {
				teststream<<i*eBinStep+eBinStep<<"  "<<nue<<"  "<<numu<<"  "<<nutau<<"  "<<nuebar<<"  "<<numubar<<"  "<<nutaubar<<std::endl;
			}
		}

	}

	
	teststream.close();
}


std::vector<double >  globes_spectrum(std::vector<double > myflux){
	std::vector<double > result;

	//char * myine = "na.glb";

	flux2GlobesFile(myflux,0);
	glbClearExperimentList();

	glbDefineAEDLVariable("glbEmin",eBinMin+eBinStep);
	glbDefineAEDLVariable("glbEmax",eBinMax+eBinStep);
	glbDefineAEDLVariable("glbEnum",eBinMax/eBinStep);
 

 	glbInitExperiment((char *)"template.glb",&glb_experiment_list[0],&glb_num_of_exps);
  
	glbSetBaselineInExperiment(0,0.54);

  	glb_params true_values_N = glbAllocParams();

  	glbDefineParams(true_values_N,theta12,theta13,theta23,deltacp,sdm,ldm_N);
  	glbSetDensityParams(true_values_N,1.0,GLB_ALL);

  	glbSetOscillationParameters(true_values_N);
	glbSetRates();

//	int l = INITIAL; // initial neutrino flavour
//	int m = FINAL; // final neutrino flavour
//	int panti = PANTI; // (anti)neutrino flag

	int i;
	int n_bins = glbGetNumberOfBins(0);
	double *signal_mu = glbGetSignalRatePtr(0,1);
	double *signal_e = glbGetSignalRatePtr(0,4);
	double *background_mu = glbGetBGRatePtr(0,1);
	double *background_e = glbGetBGRatePtr(0,4);
	double *bin_c = glbGetBinCentersListPtr(0);
	double *bin_w = glbGetBinSizeListPtr(0);

	for(i=0;i<n_bins;i++)
	{
		result.push_back(signal_e[i]+background_e[i]);
		//printf("%.5lf %.6g %.6g %.6g %.6g\n", bin_c[i]-0.5*bin_w[i], signal_mu[i], background_mu[i], signal_e[i], background_e[i]);
		//std::cout<<bin_c[i]-0.5*bin_w[i]<<"  "<<signal_e[i]+background_e[i]<<std::endl;
	}

	glbFreeParams(true_values_N);

return result;
}

std::vector<double> decaytransit (double Ms, double U, double Gdec, const std::vector<std::vector<double > > &fluxptr) {
	std::vector<double > result;
	std::vector<double > globed =  globes_spectrum(   transitspectrum(Ms,U,Gdec, fluxptr));
	for (int ee=2; ee<=38; ee=ee+2) {
		//std::cout<<ebin[ee]<<"  "<<globed[ee] + globed[ee+1]<<std::endl;
		result.push_back(    globed[ee] + globed[ee+1]   );
	}
return result;
}


// ######################################################################################################
// ############################## NORMALIZATION #########################################################


std::vector<double >  globes_norming(std::vector<double > myflux){
	std::vector<double > result;

	//char * myine = "na.glb";

	flux2GlobesFile(myflux,1);
	glbClearExperimentList();

	glbDefineAEDLVariable("glbEmin",eBinMin+eBinStep);
	glbDefineAEDLVariable("glbEmax",eBinMax+eBinStep);
	glbDefineAEDLVariable("glbEnum",eBinMax/eBinStep);
 

 	glbInitExperiment((char *)"template.glb",&glb_experiment_list[0],&glb_num_of_exps);
  
	glbSetBaselineInExperiment(0,0.54);

  	glb_params true_values_N = glbAllocParams();

  	glbDefineParams(true_values_N,theta12,theta13,theta23,deltacp,sdm,ldm_N);
  	glbSetDensityParams(true_values_N,1.0,GLB_ALL);

  	glbSetOscillationParameters(true_values_N);
	glbSetRates();

	int l = 2; // initial neutrino flavour
	int m = 1; // final neutrino flavour
	int panti = PANTI; // (anti)neutrino flag

	int i;
	int n_bins = glbGetNumberOfBins(0);
	double *signal_mu = glbGetSignalRatePtr(0,1);
	double *signal_e = glbGetSignalRatePtr(0,0);
	double *background_mu = glbGetBGRatePtr(0,1);
	double *background_e = glbGetBGRatePtr(0,0);
	double *bin_c = glbGetBinCentersListPtr(0);
	double *bin_w = glbGetBinSizeListPtr(0);

	for(i=0;i<n_bins;i++)
	{
		result.push_back(signal_mu[i]+background_mu[i]);
		//printf("%.5lf %.6g %.6g %.6g %.6g\n", bin_c[i]-0.5*bin_w[i], signal_mu[i], background_mu[i], signal_e[i], background_e[i]);
		//std::cout<<bin_c[i]-0.5*bin_w[i]<<"  "<<signal_e[i]+background_e[i]<<std::endl;
	}

	glbFreeParams(true_values_N);

return result;
}



