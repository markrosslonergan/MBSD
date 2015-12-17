#include <cstdlib>
#include <cstdio>
#include <numeric>
#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include <ctime>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "LR.h"
#include "bounds.h"

#define SIG_ON 1.0
#define SIG_OFF 0.0

#define ENERGY_FLAG 0
#define ANGULAR_FLAG 1
#define QE_FLAG 2

#define POSNU 0
#define POSNUBAR 1
#define NEGNU 2
#define NEGNUBAR 3


const char *NAMES[3] = { "Energy", "Angular", "Quasi-elastic" };

double getEvents(CL_input input, double events[][NUM_EVENT_OBS])
{
	double mS = input.mS;
	double mZprime = input.mZprime;

	FILE *ptr_file;
    	
	char buf[3000];

	char * pch;

	int n = 1;
	int m = 0;
	char s[100];
	sprintf(s,"%.3lf_%.3lf.dat", mS, mZprime);
	//char filename[500] = "../decay/data/\0";
	char filename[500] = "../decay/data12oct/anti/\0";
//	char filename[500] = "/scratch/ross/git/MBSD/decay_new/data/\0";

	strcat(filename,s);
//	printf("Filename: %s\n",filename);
	ptr_file =fopen(filename,"r");

    	if (!ptr_file)
       	{			
		printf("ERROR LOADING MC EVENTS 1\n");
		exit(1);
	}

    	while (fgets(buf,3000, ptr_file)!=NULL)
	{
		pch = strtok(buf," \t");
		n=1;
 		while (pch != NULL)
		{
			if(n==1){	//printf("%.7g, ",strtof(pch,NULL));
					events[m][0] = strtof(pch,NULL);	//E_sum 
				}
			if(n==2){	//printf("%.7g, ",strtof(pch,NULL));
					events[m][1] = strtof(pch,NULL);	//Th_sum
				}
			if(n==3){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][2] = strtof(pch,NULL);	//AngSep
				}
			if(n==4){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][3] = strtof(pch,NULL);	//E_sterile
				}
			if(n==6){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][4] = strtof(pch,NULL);	//E_high
				}
			if(n==7){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][6] = strtof(pch,NULL);	//E_low
				}
			if(n==8){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][5] = strtof(pch,NULL);	//E_low
				}

			pch = strtok(NULL," \t");
			n++;	
		}
		m++;
	}
	fclose(ptr_file);

//	printf("%s\n",flux_temp);

//printf("Total lines: %d\n", m);

return 0;
}

double getEvents_NU(CL_input input, double events[][NUM_EVENT_OBS])
{
	double mS = input.mS;
	double mZprime = input.mZprime;

	FILE *ptr_file;
    	
	char buf[3000];

	char * pch;

	int n = 1;
	int m = 0;
	char s[100];
	sprintf(s,"%.3lf_%.3lf.dat", mS, mZprime);
	//char filename[500] = "../decay/data/\0";
		char  filename[500] = "../decay/data28oct/\0";
	
//	char filename[500] = "/scratch/ross/git/MBSD/decay_new/data/\0";

	strcat(filename,s);
//	printf("Filename: %s\n",filename);
	ptr_file =fopen(filename,"r");

    	if (!ptr_file)
       	{			
		printf("ERROR LOADING MC EVENTS 1 getEvents_NU\n");
		exit(1);
	}

    	while (fgets(buf,3000, ptr_file)!=NULL)
	{
		pch = strtok(buf," \t");
		n=1;
 		while (pch != NULL)
		{
			if(n==1){	//printf("%.7g, ",strtof(pch,NULL));
					events[m][0] = strtof(pch,NULL);	//E_sum 
				}
			if(n==2){	//printf("%.7g, ",strtof(pch,NULL));
					events[m][1] = strtof(pch,NULL);	//Th_sum
				}
			if(n==3){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][2] = strtof(pch,NULL);	//AngSep
				}
			if(n==4){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][3] = strtof(pch,NULL);	//E_sterile
				}
			if(n==6){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][4] = strtof(pch,NULL);	//E_high
				}
			if(n==7){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][6] = strtof(pch,NULL);	//E_low
				}
			if(n==8){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][5] = strtof(pch,NULL);	//E_low
				}

			pch = strtok(NULL," \t");
			n++;	
		}
		m++;
	}
	fclose(ptr_file);

//	printf("%s\n",flux_temp);

//printf("Total lines: %d\n", m);

return 0;
}
double getEvents_NUBAR(CL_input input, double events[][NUM_EVENT_OBS])
{
	double mS = input.mS;
	double mZprime = input.mZprime;

	FILE *ptr_file;
    	
	char buf[3000];

	char * pch;

	int n = 1;
	int m = 0;
	char s[100];
	sprintf(s,"%.3lf_%.3lf.dat", mS, mZprime);
	//char filename[500] = "../decay/data/\0";
		char  filename[500] = "../decay/data28oct/anti/\0";
	
//	char filename[500] = "/scratch/ross/git/MBSD/decay_new/data/\0";

	strcat(filename,s);
//	printf("Filename: %s\n",filename);
	ptr_file =fopen(filename,"r");

    	if (!ptr_file)
       	{			
		printf("ERROR LOADING MC EVENTS 1 getEvents_NUBAR\n");
		exit(1);
	}

    	while (fgets(buf,3000, ptr_file)!=NULL)
	{
		pch = strtok(buf," \t");
		n=1;
 		while (pch != NULL)
		{
			if(n==1){	//printf("%.7g, ",strtof(pch,NULL));
					events[m][0] = strtof(pch,NULL);	//E_sum 
				}
			if(n==2){	//printf("%.7g, ",strtof(pch,NULL));
					events[m][1] = strtof(pch,NULL);	//Th_sum
				}
			if(n==3){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][2] = strtof(pch,NULL);	//AngSep
				}
			if(n==4){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][3] = strtof(pch,NULL);	//E_sterile
				}
			if(n==6){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][4] = strtof(pch,NULL);	//E_high
				}
			if(n==7){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][6] = strtof(pch,NULL);	//E_low
				}
			if(n==8){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][5] = strtof(pch,NULL);	//E_low
				}

			pch = strtok(NULL," \t");
			n++;	
		}
		m++;
	}
	fclose(ptr_file);

//	printf("%s\n",flux_temp);

//printf("Total lines: %d\n", m);

return 0;
}
int wipeEventArray(double events[][NUM_EVENT_OBS])
{
	int i;
	for (i=0;i<=NUMEVENTS-1;i++)
	{ 
		events[i][0]=0.0;
		events[i][1]=0.0;
		events[i][2]=0.0;
		events[i][3]=0.0;
		events[i][4]=0.0;
		events[i][5]=0.0;
	}
return 0;
}

double applyObservableCuts(CL_input input, double events[][NUM_EVENT_OBS])
{
	double ECut = input.eCut;
	double EFloor = input.eFloor;
	double ERatio = input.eRatio;
	double thCut = input.thCut;

	int i,m;
	int good=0;

	for (i=0;i<=NUMEVENTS-1;i++)
	{ 
		//For each event; does it pass cuts? If so, make sure E_sum is deposited energy and Th_sum is emission angle.
		if(events[i][0]>ECut) //Total deposited energy is greater than threshold: E_sum > ECut.
		{
			if(events[i][2] < thCut) //AngSep < thCut
			{ 
				good++;
			}
			else if (events[i][5]<ERatio*events[i][4] && events[i][5]< EFloor) //Energy asymmetry cut
			{
				events[i][1] = events[i][6];
				good++;
			}
			else // wipe event
			{
				for(m=0;m<NUM_EVENT_OBS;m++){events[i][m]=0.0;} 
			}
		}
		else // wipe event
		{
			for(m=0;m<NUM_EVENT_OBS;m++){events[i][m]=0.0;} 
		}
	}

return ((double)good)/((double)NUMEVENTS);
}


BF_RESULT * mcmc_stats_fit_spectra_indiv(CL_input in, double cutEfficiency, double events[NUMEVENTS][NUM_EVENT_OBS])
{	
	int which_var = in.which_var;
	int pos_selector = in.pos_selector;

	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r,std::time(0));

	//chi_histogrammer_indiv(in,0.01,0.02,cutEfficiency,events,which_var);


	double zeta_b =0.0;
	double sigma_s = 1.0;
	double sigma_zeta = in.Sigma_Zeta;
	
	double N=0.0;
	double lambda=0.0;
	double sum =0.0;
	
	int BINS = which_BINS(which_var);
	std::vector<double > spec_obs (BINS,0);
	std::vector<double > spec_bkg (BINS,0);
	double OBSEVENTS = 0;

	switch (which_var)
	{	
		case 0:	// Energy Spectrum
			if (pos_selector == POSNU || pos_selector == POSNUBAR){
				spec_obs = {204, 280, 214, 99, 83, 59, 51, 33, 37, 23, 19, 21, 12, 16, 4, 9, 4, 7, 3};
				spec_bkg = {151.5, 218.8, 155.6, 108.7, 72.5, 57.6, 45, 38.5, 31.4,22.2, 20.4, 17.2, 14.1, 10.2, 9.1, 8.2, 5.6, 5.7, 2.9};
			} else if (pos_selector == NEGNU || pos_selector == NEGNUBAR){
				spec_obs ={93,130,85,68,45,40,14,18,11,14,12,12,12,2,4,7,3,2,4};
				spec_bkg ={ 74.2,107.5,73.5,49.3,36.7,27.8,25.1,20.4,18.6,13.9,13.5,9.8,8.9,7.8,5.3,5,3.9,3.8,1.9};
			} 
			break;
		case 1: //Angular Spectrum
			if (pos_selector == POSNU || pos_selector == POSNUBAR){
				spec_obs = {22,34,43,41,60,87,90,139,237,429};
				spec_bkg = {19.9,23.1,28.8,32.1,46.4,63.1,86.1,121,196.8,390};
			} else if (pos_selector == NEGNU || pos_selector == NEGNUBAR){
				spec_obs = {10,13,16,20,24,36,41,70,94,263};
				spec_bkg = {9.2,11.2,13.5,16,18.7,24.2,36,52.1,94.9,237.1};
			} 
	
			break;
		case 2: //Quasi_Elastic Reco Spectrum
			if (pos_selector == POSNU || pos_selector == POSNUBAR){
				spec_obs ={232,156,156,79,81,70,63,65,62,34,70};
		    	        spec_bkg = {181.1,108.4,120.4,64.2,90.3,67.7,70.4,57.5,52.3,39,70.2};
			} else if (pos_selector == NEGNU || pos_selector == NEGNUBAR){// BROKE not actual CCQE
				spec_obs ={232,156,156,79,81,70,63,65,62,34,70};
		    	        spec_bkg = {181.1,108.4,120.4,64.2,90.3,67.7,70.4,57.5,52.3,39,70.2};
			} 
	
			break;
		default:
			std::cout<<"ERROR: which_var has to be 0-2 corresponding to E,A and QE"<<std::endl;
	}

	OBSEVENTS = std::accumulate(spec_obs.begin(), spec_obs.end(), 0); 

	double Gram[BINS];
	for(int i=0;i<BINS;i++)
	{
		Gram[i]=0.0;
	}

	std::vector<double > best_spectrum;
	double best_zeta_b = 0.0;

	double tmp_tot=0.0;
	double logchiU,chiU,contEfficiency;
	double best_contEfficiency=1e50;
	double best_N_events=1e-50;
	double best_chiU = 1e50;
	double best = 1e5;
	double temp_mS = in.mS;
	double temp_mZprime = in.mZprime;

	double N_events=0;
	double N_sig_events=0;
	double N_bg_events=0;
	double best_N_sig_events = 0;
	double best_N_bg_events = 0;

	double logchiU_step = 0.025;
	double logchiU_start = -7.0;
//#########################################################################

        std::vector<double > VGram(Gram, Gram + BINS);
        std::vector<double > bf_zeta_b;  
        double bf_chi ;

// ########################   Calculate LogLiki (Dchi^2) for Background only with SYtematics ################

        std::vector<double > Zeros(BINS, 0.0);
        std::vector<double > BF_bkg_only_zeta_b = {0};
        double BF_bkg_only_chi = 10000;
        nuisMarginalize(& BF_bkg_only_zeta_b, & BF_bkg_only_chi, &Zeros, which_var, sigma_zeta,pos_selector);
        std::cout<<"# which: "<<which_var<<" ,sigma_zeta Bf: "<< BF_bkg_only_zeta_b[0]<<" BF Chi: "<<BF_bkg_only_chi<<std::endl;

//######################################################################### BEGIN MCMC 

	int tORt = 3 ;
	int mcCount=0; 	int mcRawCount = 0; int stuckCount=0;
	int mcNumRun = 1000; 
	int endStart = 0;
	int endEnd = 2000; 	
	int mcMaxStuck = -1;
	int mcRawMax = 20000; 
	double mcSumLast=1e10;
	double mcTemp=0.4;
	double mcProb=0.0;
	double mcGenRan = gsl_rng_uniform(r);
	double finalScale = getTotalNumEvents(in);
	bool endStep = false;
	bool validBound = true;

	std::vector<double > mcS (tORt,4); 	int mcNumVar = mcS.size();
	std::vector<double > mcMin (tORt,-5);
	std::vector<double > mcBestVar (tORt,99);
	double mm = 0.0;
	std::vector<double > mcMax (tORt,mm);
		mcMax[0]= std::min(mm,log(boundBASEu(temp_mS))/log(10.0));
		mcMax[1]= std::min(0.0,log(boundBASEzp(temp_mZprime))/log(10.0));
		//mcMax[2]= log(boundBASEu(temp_mS))/log(10.0);
	std::vector<double > mcVarLast = mcMax;
	for(int i=0; i< mcNumVar; i++) 
	{
		mcVarLast[i]=mcMin[0]+gsl_rng_uniform(r)*(mcMax[0]-mcMin[0]);
	}
	std::vector<double > mcTempVar = mcVarLast;
	bool boolnone = true;

	while(( mcCount < mcNumRun || (endStep && ( (mcRawCount - endStart) < endEnd)) ) && mcRawCount < mcRawMax && boolnone ) 
	{
       		mcRawCount++;
		sum=0.0;

		// Initialise T = 0.0 end run.
		if((mcCount == mcNumRun-1 || mcRawCount == mcRawMax-endEnd) && best < 0.99*BF_bkg_only_chi){

			endStep = true;
			std::fill (mcS.begin(),mcS.end(),0.025);
			mcVarLast=mcBestVar;
			mcSumLast=best;
			endStart = mcRawCount;
			
		} 	

		// Once a non-forced fake (schwartz boundary) point is found, reduce stepsize approtiately
		if(mcSumLast < 0.99*BF_bkg_only_chi && !endStep)
		{ 	

			std::fill (mcS.begin(),mcS.end(),0.125);
		}


		if( mcNumVar == 2){
			validBound = bound_is_legit_order1(pow(10,mcTempVar[0]),pow(10,mcTempVar[1]),temp_mS,temp_mZprime);
		} else  {
		        validBound = bound_is_legit_tau(pow(10,mcTempVar[0]),pow(10,mcTempVar[2]),pow(10,mcTempVar[1]),temp_mS,temp_mZprime); 
		}

		if( !validBound || mcTempVar[0]>mcMax[0] || mcTempVar[1]>mcMax[1]  ) {  
			sum = 1e6;
			stuckCount++;
		} else {
			in.mS = temp_mS;
			in.mZprime = temp_mZprime;


			double tUd=0;	double tUp=0;	double tChi=0;



			if( mcNumVar == 2){
				tUp=pow(10,mcTempVar[0]);
				tChi=pow(10,mcTempVar[1]);
				tUd=1.0;
		
			} else {
				tUp=pow(10,mcTempVar[0]);
				tChi=pow(10,mcTempVar[1]);
				tUd=pow(10,mcTempVar[2]);
			}
		        
					
			bf_zeta_b = {0};	//re initalise for safetly
			bf_chi=10000;
	
			contEfficiency = histogrammer_indiv2(in,tUp,tUd,tChi,cutEfficiency,events,Gram,which_var,finalScale);
		        VGram.assign(Gram, Gram+BINS);
		        nuisMarginalize(&bf_zeta_b, &bf_chi, &VGram,which_var,sigma_zeta,pos_selector);
		        zeta_b = bf_zeta_b[0];
			/*for(int i =0;i<BINS;i++){
				std::cout<<zeta_b<<" "<<VGram[i]<<" "<<Gram[i]<<std::endl;
			}*/

			N_events = 0;
			N_sig_events = 0;
			N_bg_events = 0;

			int bin = 0;
			double temp_sig=0.0;
			double temp_bg =0.0 ;

			
			if(sum!=0.0){ std::cout<<"# ERROR: Sum is not Zero: "<<sum<<std::endl;}

				for(bin=0;bin<BINS;bin++)
				{
					temp_sig = Gram[bin];
					temp_bg = (1.0+zeta_b)*spec_bkg[bin];
					lambda = temp_sig + temp_bg;
					N = spec_obs[bin]; 
				
					N_events += lambda;
					N_sig_events += temp_sig;
					N_bg_events += temp_bg;
					sum= sum + 2.0*(lambda-N) + 2.0*N*log(N/lambda);
				}
				sum = sum + pow((zeta_b/sigma_zeta),2.0); 
				
				if(abs(sum-bf_chi)>1e-3){
					std::cout<<"#Error: sum != bf_chi in mcmc. Sum : "<<sum<<" bf_chi: "<<bf_chi<<std::endl;
				}
		}

		if(sum < best)			
		{ 
			best = sum; 
			best_contEfficiency = contEfficiency; 
			best_N_events = N_events;
			best_N_sig_events = N_sig_events;
			best_N_bg_events = N_bg_events;
			best_spectrum = VGram;
			best_zeta_b = zeta_b;
			for(int i=0; i < mcNumVar;i++){ 
				mcBestVar[i] = mcTempVar[i];
			}
		}


		
		if(!endStep){

			mcGenRan = gsl_rng_uniform(r);
			//Probability of accepting a new point in the chain.
			mcProb = exp(-(sum-mcSumLast)/mcTemp); 

			if(-(sum-mcSumLast)/mcTemp < -10.0)
			{ 
				mcProb = 0;
			} else if( sum < mcSumLast)
			{ 
				mcProb = 1;
			}
			

//std::cout<<mcCount<<" T: "<<mcTemp<<" Step "<<mcS[0]<<" random "<<mcGenRan<<" prob: "<<mcProb<<" exp "<<-(sum-mcSumLast)/mcTemp<<"  sum "<<sum<<" sumLast "<<mcSumLast<<" var ";
//mcPrintVar(mcTempVar);
//std::cout<<" BEST: "<<best<<" raw#: "<<mcRawCount<<" NOR"<<std::endl;
		
			if(mcGenRan < mcProb) { 
				mcSumLast=sum;
				mcVarLast=mcTempVar;
				mcCount++;
				sum = 0.0;
			}
			  //otherwise choose a new point

		} else {  // If in Endstep

//std::cout<<mcCount<<" T: "<<mcTemp<<" Step "<<mcS[0]<<" random "<<mcGenRan<<" prob: "<<mcProb<<" exp "<<-(sum-mcSumLast)/mcTemp<<"  sum "<<sum<<" sumLast "<<mcSumLast<<" var ";
//mcPrintVar(mcTempVar);
//std::cout<<" BEST: "<<best<<" raw#: "<<mcRawCount<<" END"<<std::endl;

			  mcTemp= 1/(0.1*double((mcRawCount - endStart)+1)+1)+0.05;
			  if(sum < mcSumLast){ 
				mcSumLast=sum;
				mcVarLast=mcTempVar;
				mcCount++;
				sum = 0.0;
			} 

		} 	



		for(int i=0; i < mcNumVar;i++){ //Initialise MCMC variables
			 mcGenRan= gsl_rng_uniform(r);
			 int test = 0;
			 mcTempVar[i] = mcVarLast[i]+mcS[i]*(mcGenRan-0.5)*(mcMin[i]-mcMax[i]);
			 while((mcTempVar[i] > mcMax[i] || mcTempVar[i] < mcMin[i] ) && boolnone ){	
						 mcGenRan = gsl_rng_uniform(r);
						 test++; 
						 if(test > 50000){
							boolnone = false;
						 }
						 mcTempVar[i] = mcVarLast[i]+mcS[i]*(mcGenRan-0.5)*(mcMin[i]-mcMax[i]);
					//	 std::cout<<"Thats the one: No points at all in whole range!: "<<test<<std::endl;
			 } 	

		}

		sum=0.0;


}//End MCMC loop
	double decaywidth=Gsterile2vvv(mcBestVar[0],mcBestVar[2],mcBestVar[1],temp_mS,temp_mZprime);
	double valGoF1=GoF1(best_spectrum,spec_obs,spec_bkg,best_zeta_b);
	double valGoF2=GoF2(best_spectrum,spec_obs,spec_bkg,best_zeta_b);

	std::cout<<"#From Mark: "<<getTotalNumEvents(in)<<"; Variable:  "<<NAMES[which_var]<<";  "<<best_N_events<<" = ("<<best_N_sig_events<<" + "<<best_N_bg_events<<")"<<std::endl;
	std::cout<<"# Decay Related Checks: Gtotal (GeV): "<<decaywidth<<" 1/Gtotal (m): "<<G2M(pow(decaywidth,-1))<<" 1/Gtotal (s): "<<invG2S(pow(decaywidth,-1))<<std::endl;
       	std::cout<<"# GoF tests. GoF1: "<<valGoF1<<" GoF1/nDoF: "<<valGoF1/spec_obs.size()<<" pVal1: "<<pval(valGoF1,spec_obs.size())<<" GoF2: "<<pow(best_N_events - OBSEVENTS ,2)/OBSEVENTS<<std::endl;
	std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<BF_bkg_only_chi-best<<" ";
	mcPrintVar(mcBestVar);
	std::cout<<" "<<best_contEfficiency<<" "<<cutEfficiency<<" "<<best_zeta_b<<std::endl;
gsl_rng_free(r);

BF_RESULT* output = new BF_RESULT();
//BF_RESULT * output = (BF_RESULT *)malloc(sizeof(BF_RESULT));

output->stats_bf=best_zeta_b;

switch(which_var) 
{
	case 0:
		output->E_spectrum = best_spectrum;
		output->E_bf_Up = pow(10,mcBestVar[0]);
		output->E_bf_Chi = pow(10,mcBestVar[1]);
		if(mcNumVar == 2) {
			output->E_bf_Ud = 1.0;
		} else {
			output->E_bf_Ud = pow(10,mcBestVar[2]);
		}

		break;
	case 1:
		output->A_spectrum = best_spectrum;
		output->A_bf_Up = pow(10,mcBestVar[0]);
		output->A_bf_Chi = pow(10,mcBestVar[1]);
		if(mcNumVar == 2) {
			output->A_bf_Ud = 1.0;
		} else {
			output->A_bf_Ud = pow(10,mcBestVar[2]);
		}

		break;
	case 2:	
		output->QE_spectrum= best_spectrum;
		output->QE_bf_Up = pow(10,mcBestVar[0]);
		output->QE_bf_Chi = pow(10,mcBestVar[1]);
		if(mcNumVar == 2) {
			output->QE_bf_Ud = 1.0;
		} else {
			output->QE_bf_Ud = pow(10,mcBestVar[2]);
		}

		break;
	default:							
		std::cout<<"ERROR: which_var in mcmc must be 0(E), 1(A) or 2(QE)"<<std::endl;
}



return output;

}

BF_RESULT * mcmc_stats_dual(CL_input in, double cutEfficiency_NU, double cutEfficiency_NUBAR, double events_NU[NUMEVENTS][NUM_EVENT_OBS], double events_NUBAR[NUMEVENTS][NUM_EVENT_OBS])
{	
	int which_var = in.which_var;
	int pos_selector = in.pos_selector;

	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set(r,std::time(0));

	//chi_histogrammer_indiv(in,0.01,0.02,cutEfficiency,events,which_var);


	double zeta_b =0.0;
	double zeta_b_NUBAR = 0.0;
	double sigma_s = 1.0;
	double sigma_zeta = in.Sigma_Zeta;
	
	double N=0.0;
	double lambda=0.0;
	double sum =0.0;
	
	int BINS = 2*which_BINS(which_var);
	std::vector<double > spec_obs (BINS,0);
	std::vector<double > spec_bkg (BINS,0);
	double OBSEVENTS = 0;


	switch (which_var)
	{	
		case 0:	// Energy Spectrum
				spec_obs = {204, 280, 214, 99, 83, 59, 51, 33, 37, 23, 19, 21, 12, 16, 4, 9, 4, 7, 3,93,130,85,68,45,40,14,18,11,14,12,12,12,2,4,7,3,2,4};
				spec_bkg = {151.5, 218.8, 155.6, 108.7, 72.5, 57.6, 45, 38.5, 31.4,22.2, 20.4, 17.2, 14.1, 10.2, 9.1, 8.2, 5.6, 5.7, 2.9, 74.2,107.5,73.5,49.3,36.7,27.8,25.1,20.4,18.6,13.9,13.5,9.8,8.9,7.8,5.3,5,3.9,3.8,1.9};
			break;
		case 1: //Angular Spectrum
				spec_obs = {22,34,43,41,60,87,90,139,237,429,10,13,16,20,24,36,41,70,94,263};
				spec_bkg = {19.9,23.1,28.8,32.1,46.4,63.1,86.1,121,196.8,390,9.2,11.2,13.5,16,18.7,24.2,36,52.1,94.9,237.1};
			break;
		case 2: //Quasi_Elastic Reco Spectrum
				spec_obs ={232,156,156,79,81,70,63,65,62,34,70,232,156,156,79,81,70,63,65,62,34,70};
		    	        spec_bkg = {181.1,108.4,120.4,64.2,90.3,67.7,70.4,57.5,52.3,39,70.,181.1,108.4,120.4,64.2,90.3,67.7,70.4,57.5,52.3,39,70.2};
			break;
		default:
			std::cout<<"ERROR: which_var has to be 0-2 corresponding to E,A and QE"<<std::endl;
	}

	OBSEVENTS = std::accumulate(spec_obs.begin(), spec_obs.end(), 0); 

	double Gram[BINS];
	for(int i=0;i<BINS;i++)
	{
		Gram[i]=0.0;
	}

	std::vector<double > best_spectrum = spec_obs;
	double best_zeta_b = 0.0;	double best_zeta_b_NUBAR = 0.0;

	
	double tmp_tot=0.0;
	double logchiU,chiU,contEfficiency;
	double best_contEfficiency=1e50;
	double best_N_events=1e-50;
	double best_chiU = 1e50;
	double best = 1e5;
	double temp_mS = in.mS;
	double temp_mZprime = in.mZprime;

	double N_events=0;
	double N_sig_events=0;
	double N_bg_events=0;
	double best_N_sig_events = 0;
	double best_N_bg_events = 0;

//#########################################################################

        std::vector<double > VGram(Gram, Gram + BINS);
        std::vector<double > bf_zeta_b;  
        double bf_chi ;

	if(VGram.size() != BINS && BINS != spec_obs.size()){std::cout<<"ERROR, dimesnons mismatch! "<<BINS<<" "<<VGram.size()<<" "<<spec_obs.size()<<std::endl;}


// ########################   Calculate LogLiki (Dchi^2) for Background only with SYtematics ################

        std::vector<double > Zeros(BINS, 0.0);
        std::vector<double > BF_bkg_only_zeta_b = {0.0,0.0};
        double BF_bkg_only_chi = 10000;
	
        nuisMarginalize_dual(& BF_bkg_only_zeta_b, & BF_bkg_only_chi, &Zeros, which_var, sigma_zeta);
        std::cout<<"# which: "<<which_var<<" ,sigma_zeta Bf: "<< BF_bkg_only_zeta_b[0]<<" BF Chi: "<<BF_bkg_only_chi<<std::endl;

//######################################################################### BEGIN MCMC 

	int tORt = 3 ;
	int mcCount=0; 	int mcRawCount = 0; int stuckCount=0;
	int mcNumRun = 3000; //1000 
	int endStart = 0;
	int endEnd = 6000; 	//2000
	int mcMaxStuck = -1;
	int mcRawMax = 30000; //20000
	double mcSumLast=1e10;
	double mcTemp=0.4;
	double mcProb=0.0;
	double mcGenRan = gsl_rng_uniform(r);
	double finalScale_NU = getTotalNumEvents_NU(in);
	double finalScale_NUBAR = getTotalNumEvents_NUBAR(in);
	bool endStep = false;
	bool validBound = true;

	std::vector<double > mcS (tORt,4); 	int mcNumVar = mcS.size();
	std::vector<double > mcMin (tORt,-5);
	std::vector<double > mcBestVar (tORt,99);
	double mm = -0.6505; // a lowest of Ut4^2=0.05 
	mcMin[2]=-3;
	std::vector<double > mcMax (tORt,mm);
		mcMax[0]= std::min(mm,log(boundBASEu(temp_mS))/log(10.0));
		mcMax[1]= std::min(-1.0,log(boundBASEzp(temp_mZprime))/log(10.0));
		//mcMax[2]= log(boundBASEu(temp_mS))/log(10.0);
	std::vector<double > mcVarLast = mcMax;
	for(int i=0; i< mcNumVar; i++) 
	{
		mcVarLast[i]=mcMin[0]+gsl_rng_uniform(r)*(mcMax[0]-mcMin[0]);
	}
	std::vector<double > mcTempVar = mcVarLast;
	bool boolnone = true;

	while(( mcCount < mcNumRun || (endStep && ( (mcRawCount - endStart) < endEnd)) ) && mcRawCount < mcRawMax && boolnone ) 
	{
       		mcRawCount++;
		sum=0.0;
		
		//std::cout<<mcRawCount<<" ";
		//mcPrintVar(mcTempVar);	
		//std::cout<<"  "<<mcSumLast<<std::endl;

		// Initialise T = 0.0 end run.
		if((mcCount == mcNumRun-1 || mcRawCount == mcRawMax-endEnd) && best < 0.99*BF_bkg_only_chi){

			endStep = true;
			std::fill (mcS.begin(),mcS.end(),0.025);
			mcVarLast=mcBestVar;
			mcSumLast=best;
			endStart = mcRawCount;
			
		} 	

		// Once a non-forced fake (schwartz boundary) point is found, reduce stepsize approtiately
		if(mcSumLast < 0.99*BF_bkg_only_chi && !endStep)
		{ 	

			std::fill (mcS.begin(),mcS.end(),0.125);
		}


		if( mcNumVar == 2){
			validBound = bound_is_legit_order1(pow(10,mcTempVar[0]),pow(10,mcTempVar[1]),temp_mS,temp_mZprime);
		} else  {
		        validBound = bound_is_legit_tau(pow(10,mcTempVar[0]),pow(10,mcTempVar[2]),pow(10,mcTempVar[1]),temp_mS,temp_mZprime); 
		}

		if( !validBound || mcTempVar[0]>mcMax[0] || mcTempVar[1]>mcMax[1]  ) {  
			sum = 1e6;
			stuckCount++;
		} else {
			in.mS = temp_mS;
			in.mZprime = temp_mZprime;


			double tUd=0;	double tUp=0;	double tChi=0;



			if( mcNumVar == 2){
				tUp=pow(10,mcTempVar[0]);
				tChi=pow(10,mcTempVar[1]);
				tUd=1.0;
		
			} else {
				tUp=pow(10,mcTempVar[0]);
				tChi=pow(10,mcTempVar[1]);
				tUd=pow(10,mcTempVar[2]);
			}
		        
					
			bf_zeta_b = {0};	//re initalise for safetly
			bf_chi=10000;
	
			contEfficiency = histogrammer_indiv2_dual(in,tUp,tUd,tChi,cutEfficiency_NU,cutEfficiency_NUBAR,events_NU,events_NUBAR,Gram,which_var,finalScale_NU,finalScale_NUBAR);
		        
			VGram.assign(Gram, Gram+BINS);
		        nuisMarginalize_dual(&bf_zeta_b, &bf_chi, &VGram,which_var,sigma_zeta);
		        zeta_b = bf_zeta_b[0];
			zeta_b_NUBAR = bf_zeta_b[1];
/*			for(int i =0;i<BINS;i++){
				std::cout<<zeta_b<<" "<<VGram[i]<<" "<<Gram[i]<<std::endl;
			}
*/

			N_events = 0;
			N_sig_events = 0;
			N_bg_events = 0;

			int bin = 0;
			double temp_sig=0.0;
			double temp_bg =0.0 ;

			
			if(sum!=0.0){ std::cout<<"# ERROR: Sum is not Zero: "<<sum<<std::endl;}

				for(bin=0;bin<BINS/2;bin++)
				{
					temp_sig = Gram[bin];
					temp_bg = (1.0+zeta_b)*spec_bkg[bin];
					lambda = temp_sig + temp_bg;
					N = spec_obs[bin]; 
				
					N_events += lambda;
					N_sig_events += temp_sig;
					N_bg_events += temp_bg;
					sum= sum + 2.0*(lambda-N) + 2.0*N*log(N/lambda);
				}
				for(bin=BINS/2;bin<BINS;bin++)
				{
					temp_sig = Gram[bin];
					temp_bg = (1.0+zeta_b_NUBAR)*spec_bkg[bin];
					lambda = temp_sig + temp_bg;
					N = spec_obs[bin]; 
				
					N_events += lambda;
					N_sig_events += temp_sig;
					N_bg_events += temp_bg;
					sum= sum + 2.0*(lambda-N) + 2.0*N*log(N/lambda);
				}
				sum = sum + pow((zeta_b/sigma_zeta),2.0)+pow((zeta_b_NUBAR/sigma_zeta),2.0); 
				
				if(abs(sum-bf_chi)>1e-3){
					std::cout<<"# Error: sum != bf_chi in mcmc. Sum : "<<sum<<" bf_chi: "<<bf_chi<<std::endl;
				}
		}

		if(sum < best)			
		{ 
			best = sum; 
			best_contEfficiency = contEfficiency; 
			best_N_events = N_events;
			best_N_sig_events = N_sig_events;
			best_N_bg_events = N_bg_events;
			best_spectrum = VGram;
			best_zeta_b = zeta_b;
			best_zeta_b_NUBAR = zeta_b_NUBAR;
			for(int i=0; i < mcNumVar;i++){ 
				mcBestVar[i] = mcTempVar[i];
			}
		}


		
		if(!endStep){

			mcGenRan = gsl_rng_uniform(r);
			//Probability of accepting a new point in the chain.
			mcProb = exp(-(sum-mcSumLast)/mcTemp); 

			if(-(sum-mcSumLast)/mcTemp < -10.0)
			{ 
				mcProb = 0;
			} else if( sum < mcSumLast)
			{ 
				mcProb = 1;
			}
			

//std::cout<<mcCount<<" T: "<<mcTemp<<" Step "<<mcS[0]<<" random "<<mcGenRan<<" prob: "<<mcProb<<" exp "<<-(sum-mcSumLast)/mcTemp<<"  sum "<<sum<<" sumLast "<<mcSumLast<<" var ";
//mcPrintVar(mcTempVar);
//std::cout<<" BEST: "<<best<<" raw#: "<<mcRawCount<<" NOR"<<std::endl;
		
			if(mcGenRan < mcProb) { 
				mcSumLast=sum;
				mcVarLast=mcTempVar;
				mcCount++;
				sum = 0.0;
			}
			  //otherwise choose a new point

		} else {  // If in Endstep

//std::cout<<mcCount<<" T: "<<mcTemp<<" Step "<<mcS[0]<<" random "<<mcGenRan<<" prob: "<<mcProb<<" exp "<<-(sum-mcSumLast)/mcTemp<<"  sum "<<sum<<" sumLast "<<mcSumLast<<" var ";
//mcPrintVar(mcTempVar);
//std::cout<<" BEST: "<<best<<" raw#: "<<mcRawCount<<" END"<<std::endl;

			  mcTemp= 1/(0.1*double((mcRawCount - endStart)+1)+1)+0.05;
			  if(sum < mcSumLast){ 
				mcSumLast=sum;
				mcVarLast=mcTempVar;
				mcCount++;
				sum = 0.0;
			} 

		} 	



		for(int i=0; i < mcNumVar;i++){ //Initialise MCMC variables
			 mcGenRan= gsl_rng_uniform(r);
			 int test = 0;
			 mcTempVar[i] = mcVarLast[i]+mcS[i]*(mcGenRan-0.5)*(mcMin[i]-mcMax[i]);
			 while((mcTempVar[i] > mcMax[i] || mcTempVar[i] < mcMin[i] ) && boolnone ){	
						 mcGenRan = gsl_rng_uniform(r);
						 test++; 
						 if(test > 50000){
							boolnone = false;
						 }
						 mcTempVar[i] = mcVarLast[i]+mcS[i]*(mcGenRan-0.5)*(mcMin[i]-mcMax[i]);
					//	 std::cout<<"Thats the one: No points at all in whole range!: "<<test<<std::endl;
			 } 	

		}

		sum=0.0;


}//End MCMC loop
	double decaywidth=Gsterile2vvv(pow(10,mcBestVar[0]),pow(10,mcBestVar[2]),pow(10,mcBestVar[1]),temp_mS,temp_mZprime);
	double valGoF1=GoF1(best_spectrum,spec_obs,spec_bkg,best_zeta_b);
	//double valGoF2=GoF2(best_spectrum,spec_obs,spec_bkg,best_zeta_b);

	std::cout<<"# From Mark: nu"<<getTotalNumEvents_NU(in)<<" nubar: "<<getTotalNumEvents_NUBAR(in)<<"; Variable:  "<<NAMES[which_var]<<";  "<<best_N_events<<" = ( "<<best_N_sig_events<<" + "<<best_N_bg_events<<" )"<<std::endl;
	std::cout<<"# Decay Related Checks: Gtotal (GeV): "<<decaywidth<<" 1/Gtotal (m): "<<G2M(pow(decaywidth,-1))<<" 1/Gtotal (s): "<<invG2S(pow(decaywidth,-1))<<std::endl;
       	std::cout<<"# GoF tests. GoF1: "<<valGoF1<<" GoF1/nDoF: "<<valGoF1/spec_obs.size()<<" pVal1: "<<pval(valGoF1,spec_obs.size())<<" GoF2: "<<pow(best_N_events - OBSEVENTS ,2)/OBSEVENTS<<std::endl;
	std::cout<<"# bf spectrum: ";
	mcPrintVar(best_spectrum);
	std::cout<<std::endl;

	std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<BF_bkg_only_chi-best<<" ";
	mcPrintVar(mcBestVar);
	std::cout<<" "<<best_contEfficiency<<" "<<cutEfficiency_NU<<" "<<best_zeta_b<<" "<<best_zeta_b_NUBAR<<std::endl;


gsl_rng_free(r);

BF_RESULT* output = new BF_RESULT();
//BF_RESULT * output = (BF_RESULT *)malloc(sizeof(BF_RESULT));

output->stats_bf=best_zeta_b;

switch(which_var) 
{
	case 0:
		output->E_spectrum = best_spectrum;
		output->E_bf_Up = pow(10,mcBestVar[0]);
		output->E_bf_Chi = pow(10,mcBestVar[1]);
		if(mcNumVar == 2) {
			output->E_bf_Ud = 1.0;
		} else {
			output->E_bf_Ud = pow(10,mcBestVar[2]);
		}

		break;
	case 1:
		output->A_spectrum = best_spectrum;
		output->A_bf_Up = pow(10,mcBestVar[0]);
		output->A_bf_Chi = pow(10,mcBestVar[1]);
		if(mcNumVar == 2) {
			output->A_bf_Ud = 1.0;
		} else {
			output->A_bf_Ud = pow(10,mcBestVar[2]);
		}

		break;
	case 2:	
		output->QE_spectrum= best_spectrum;
		output->QE_bf_Up = pow(10,mcBestVar[0]);
		output->QE_bf_Chi = pow(10,mcBestVar[1]);
		if(mcNumVar == 2) {
			output->QE_bf_Ud = 1.0;
		} else {
			output->QE_bf_Ud = pow(10,mcBestVar[2]);
		}

		break;
	default:							
		std::cout<<"ERROR: which_var in mcmc must be 0(E), 1(A) or 2(QE)"<<std::endl;
}



return output;

}


BF_RESULT * plot_surface(CL_input in, double cutEfficiency, double events[NUMEVENTS][NUM_EVENT_OBS])
{	
	int which_var = in.which_var;

	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);


	//chi_histogrammer_indiv(in,0.01,0.02,cutEfficiency,events,which_var);


	double zeta_b =0.0;
	double sigma_s = 1.0;
	double sigma_zeta = in.Sigma_Zeta;
	int pos_selector=in.pos_selector;	
	double N=0.0;
	double lambda=0.0;
	double sum =0.0;
	
	int BINS = which_BINS(which_var);
	double spec_obs[BINS];
	double spec_bkg[BINS];

	switch (which_var)
	{	
	case 0:
		spec_obs[0] = 204.0;spec_obs[1] = 280.0;spec_obs[2] = 214.0;spec_obs[3] = 99.0;spec_obs[4] = 83.0;spec_obs[5] = 59.0;spec_obs[6] = 51.0;spec_obs[7] = 33.0;spec_obs[8] = 37.0;spec_obs[9] = 23.0; spec_obs[10] = 19.0;spec_obs[11] = 21.0;spec_obs[12] = 12.0; spec_obs[13] = 16.0; spec_obs[14] = 4.0;spec_obs[15] = 9.0; spec_obs[16] = 4.0; spec_obs[17] = 7.0; spec_obs[18] = 3.0;
		spec_bkg[0] = 151.5;	spec_bkg[1] = 218.8;	spec_bkg[2] = 155.6;spec_bkg[3] = 108.7;spec_bkg[4] = 72.5;spec_bkg[5] = 57.6;spec_bkg[6] = 45;spec_bkg[7] = 38.5;spec_bkg[8] = 31.4;spec_bkg[9] = 22.2;spec_bkg[10] = 20.4;spec_bkg[11] = 17.2;		spec_bkg[12] = 14.1;	spec_bkg[13] = 10.2;	spec_bkg[14] = 9.1;	spec_bkg[15] = 8.2;	spec_bkg[16] = 5.6;	spec_bkg[17] = 5.7;	spec_bkg[18] = 2.9;
		break;
	case 1: 
		spec_obs[0] = 22; 	spec_obs[1] = 34; 	spec_obs[2] = 43; 	spec_obs[3] = 41; 	spec_obs[4] = 60; 	spec_obs[5] = 87; 	spec_obs[6] = 90; 	spec_obs[7] = 139; 	spec_obs[8] = 237; 	spec_obs[9] = 429;  	
	spec_bkg[0] = 19.9; 	spec_bkg[1] = 23.1; 	spec_bkg[2] = 28.8; 	spec_bkg[3] = 32.1; 	spec_bkg[4] = 46.4; 	spec_bkg[5] = 63.1; 	spec_bkg[6] = 86.1; 	spec_bkg[7] = 121; 	spec_bkg[8] = 196.8; 	spec_bkg[9] = 390;
		break;

	case 2:
	      	spec_obs[0] =  232; spec_obs[1] = 156 ;spec_obs[2] = 156 ;spec_obs[3] = 79 ;spec_obs[4] = 81 ;spec_obs[5] = 70 ;spec_obs[6] = 63 ;spec_obs[7] = 65 ;spec_obs[8] = 62 ;spec_obs[9] = 34 ;spec_obs[10] = 70;    
	     spec_bkg[0] = 181.1 ;spec_bkg[1] = 108.4;spec_bkg[2] = 120.4;spec_bkg[3] = 64.2;spec_bkg[4] = 90.3; spec_bkg[5] = 67.7;spec_bkg[6] = 70.4;spec_bkg[7] = 57.5;spec_bkg[8] = 52.3; spec_bkg[9] = 39;spec_bkg[10] = 70.2;
		break;
	default:
		std::cout<<"ERROR: which_var has to be 0-2"<<std::endl;


	}

	double Gram[BINS];
	for(int i=0;i<BINS;i++)
	{
		Gram[i]=0.0;
	}
	double tmp_tot=0.0;

	double logchiU,chiU,contEfficiency;
	double best_contEfficiency=1e50;
	double best_N_events=1e-50;
	double best_chiU = 1e50;
	double best = 1e5;
	double temp_mS = in.mS;
	double temp_mZprime = in.mZprime;

	double N_events=0;
	double N_sig_events=0;
	double N_bg_events=0;
	double best_N_sig_events = 0;
	double best_N_bg_events = 0;

	double logchiU_step = 0.025;
	double logchiU_start = -7.0;
//#########################################################################

        std::vector<double > VGram(Gram, Gram + BINS);
        std::vector<double > bf_zeta_b;  
        double bf_chi ;

// ########################   Calculate LogLiki (Dchi^2) for Background only with SYtematics ################

        std::vector<double > Zeros(EBINS, 0.0);
        std::vector<double > BF_bkg_only_zeta_b = {0};
        double BF_bkg_only_chi = 10000;
        nuisMarginalize(& BF_bkg_only_zeta_b,& BF_bkg_only_chi, &Zeros, which_var, sigma_zeta,pos_selector);
        std::cout<<"# which: "<<which_var<<" ,sigma_zeta Bf: "<< BF_bkg_only_zeta_b[0]<<" BF Chi: "<<BF_bkg_only_chi<<std::endl;

//######################################################################### BEGIN MCMC

	double fudge1=1;
	int tORt = 3 ;

	std::vector<double > mcS (tORt,0.5); 	int mcNumVar = mcS.size();
	std::vector<double > mcMin (tORt,-5);
	std::vector<double > mcBestVar (tORt,99);

	std::vector<double > mcMax (tORt,0.0);
		mcMax[0]= fudge1*log(boundBASEu(temp_mS))/log(10.0);
		mcMax[1]= fudge1*log(boundBASEzp(temp_mZprime))/log(10.0);
	

	std::vector<double > mcVarLast = mcMax;
	for(int i=0; i< mcNumVar; i++) 
	{
		mcVarLast[i]=mcMin[0]+gsl_rng_uniform(r)*(mcMax[0]-mcMin[0]);
	}

	std::vector<double > mcTempVar = mcVarLast;




	int mcCount=0; 
	int mcRawCount = 0;
	int mcNumRun = 1000; //200
	double mcSumLast=1e10;
	double mcTemp=0.4;
	double mcProb=0.0;

	bool endStep = false;
	int endStart = 0;
	int endEnd = 2000; //How many function calls after main run, for endPhase 2000
	

	int stuckCount=0; int mcMaxStuck = -1;
	int mcRawMax = 25000; //20k

	double mcGenRan = gsl_rng_uniform(r);
	double finalScale = getTotalNumEvents(in);


	bool impbool=false;
	int impcount=0;
	double tUd=0;	double tUp=0;	double tChi=0;
	double impUd=0;	double impUp=0;	double impChi=0;
while(true){
	//for(tUp=-4;tUp<=mcMax[0];tUp=tUp+0.02){
	//for(tUd=-4;tUd<=mcMax[2];tUd=tUd+0.02){
	//for(tChi=-4;tChi<=mcMax[1];tChi=tChi+0.02){
		if(sum < 100 && !impbool){
			impbool = true;
			impcount=0;
			impUp=tUp;
			impUd=tUd;
			impChi=tChi;
		}		

		if(impbool){
			tUp=-4.0+gsl_rng_uniform(r)*(mcMax[0]+4);
			tChi=-4.0+gsl_rng_uniform(r)*(mcMax[1]+4);
			tUd=-4.0+gsl_rng_uniform(r)*(0.0+4);
			impcount++;
			if(impcount==200){ impbool=false;}
		} else {
			tUp=-4.0+gsl_rng_uniform(r)*(mcMax[0]+4);
			tChi=-4.0+gsl_rng_uniform(r)*(mcMax[1]+4);
			tUd=-4.0+gsl_rng_uniform(r)*(0.0+4);
		}


			in.mS = temp_mS;
			in.mZprime = temp_mZprime;

			sum=0.0;

	
			contEfficiency = histogrammer_indiv2(in,pow(10,tUp),pow(10,tUd),pow(10,tChi),cutEfficiency,events,Gram,which_var,finalScale);
		       	VGram.assign(Gram, Gram+BINS);
		        nuisMarginalize(&bf_zeta_b, &bf_chi, &VGram,which_var,sigma_zeta,pos_selector);
		        zeta_b = bf_zeta_b[0];

			N_events = 0;
			N_sig_events = 0;
			N_bg_events = 0;

			int bin = 0;
			double temp_sig,temp_bg;


		bool failedBound = false;
		if( mcNumVar == 2){
			failedBound = !bound_is_legit_order1(pow(10,tUp),pow(10,tChi),temp_mS,temp_mZprime);
		} else {
		        failedBound = !bound_is_legit_tau(pow(10,tUp),pow(10,tUd),pow(10,tChi),temp_mS,temp_mZprime); 
		}


		if(failedBound ) {  
			// Then it will never be able to decay before so be below the bound, so while its not lets		
			// make a fake chi (bad) for it so MCMC never visits it.
			sum = 1e6;
			stuckCount++;

		} else {
				for(bin=0;bin<BINS;bin++)
				{
					temp_sig = sigma_s*Gram[bin];
					temp_bg = (1.0+zeta_b)*spec_bkg[bin];
					lambda = temp_sig + temp_bg;
					N = spec_obs[bin]; //MB has seen O[] events.
				
					N_events += lambda;
					N_sig_events += temp_sig;
					N_bg_events += temp_bg;
					sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
				}
			
				sum += pow((zeta_b/sigma_zeta),2.0); //add on sys bkg
		}

		
		std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<tUp<<" "<<tUd<<" "<<tChi<<" "<<sum<<std::endl;
			
}// end infinite while
//}
//}
//} //end fors


	//std::cout<<"# From Mark: "<<getTotalNumEvents(in)<<"; Variable:  "<<NAMES[which_var]<<";  "<<best_N_events<<" = ("<<best_N_sig_events<<" + "<<best_N_bg_events<<")"<<std::endl;

      // 	std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<BF_bkg_only_chi-best<<" ";
	//mcPrintVar(mcBestVar);
	//std::cout<<" "<<best_contEfficiency<<" "<<cutEfficiency<<std::endl;



gsl_rng_free(r);

BF_RESULT* output = new BF_RESULT();
//BF_RESULT * output = (BF_RESULT *)malloc(sizeof(BF_RESULT));

switch(which_var) 
{
	case 0:
		output->E_bf_Up = pow(10,mcBestVar[0]);
		output->E_bf_Chi = pow(10,mcBestVar[1]);
		if(mcNumVar == 2) {
			output->E_bf_Ud = 1.0;
		} else {
			output->E_bf_Ud = pow(10,mcBestVar[2]);
		}

		break;
	case 1:
		output->A_bf_Up = pow(10,mcBestVar[0]);
		output->A_bf_Chi = pow(10,mcBestVar[1]);
		if(mcNumVar == 2) {
			output->A_bf_Ud = 1.0;
		} else {
			output->A_bf_Ud = pow(10,mcBestVar[2]);
		}

		break;
	case 2:
		output->QE_bf_Up = pow(10,mcBestVar[0]);
		output->QE_bf_Chi = pow(10,mcBestVar[1]);
		if(mcNumVar == 2) {
			output->QE_bf_Ud = 1.0;
		} else {
			output->QE_bf_Ud = pow(10,mcBestVar[2]);
		}

		break;
	default:							
		std::cout<<"ERROR: which_var in mcmc must be 0(E), 1(A) or 2(QE)"<<std::endl;
}



return output;

}

double spectra_indiv(CL_input in, double Up, double Ud, double Chi)
{	
	int which_var = in.which_var;
	double sigma_zeta = in.Sigma_Zeta;
	double finalScale = 0;
	int pos_selector = in.pos_selector;
	static double events[NUMEVENTS][NUM_EVENT_OBS];
	wipeEventArray(events); //wipe the array.

	if(pos_selector==POSNU || pos_selector == POSNUBAR){
		getEvents_NU(in,events); // populates event array. 
		finalScale = getTotalNumEvents_NU(in);
	} else if (pos_selector == NEGNU || pos_selector ==NEGNUBAR){
		getEvents_NUBAR(in,events);
		finalScale = getTotalNumEvents_NUBAR(in);
	}
	
	double cutEfficiency = applyObservableCuts(in,events); //applys cuts. (Post-cut events may have lots of blank slots).

	int BINS = which_BINS(which_var);

	double Gram[BINS];
	for(int i=0;i<BINS;i++)
	{
		Gram[i]=0.0;
	}

        std::vector<double > VGram(Gram, Gram + BINS);
	std::vector<double > zeta_b = {0,0};  
        double bf_chi = 10000; 

	//	std::cout<<" whichvar: "<<which_var<<" sigma_eta: "<<sigma_zeta<<" finalscale: "<<finalScale<<" cutEfficiency "<<cutEfficiency<<std::endl;
	std::vector<double > Zeros(EBINS, 0.0);
        std::vector<double > bkg_only_zeta_b = {0,0};
        double bkg_only_chi = 10000;
        nuisMarginalize(&bkg_only_zeta_b, &bkg_only_chi, &Zeros, which_var, sigma_zeta,pos_selector);

	histogrammer_indiv2(in,Up,Ud,Chi,cutEfficiency,events,Gram,which_var,finalScale);

        VGram.assign(Gram, Gram+BINS);
        nuisMarginalize(&zeta_b, &bf_chi, &VGram,which_var,sigma_zeta,pos_selector);
	bool validBound = false;
        validBound = bound_is_legit_tau(Up,Ud,Chi,in.mS, in.mZprime);
	if (!validBound || Up > boundBASEu(in.mS) || Chi > boundBASEzp(in.mZprime))
	{
		bf_chi=1000+bkg_only_chi; 
	}

	
	return bkg_only_chi-bf_chi;
}



int main(int argc, char * argv[])
{

	//std::cout<<whatsmaxUXorder1(0.09,0.3)<<"  "<<whatsmaxUXorder1(0.03,0.6)<<"  "<<whatsmaxUXorder1(0.05,0.45)<<std::endl;		

	static CL_input in;

	in.eCut = 0.0;
	in.thCut = 180.0;
	in.eFloor = 0.0; // I believe 0.0 for Floor and Ratio mean that the cut is never passed.
	in.eRatio = 0.0;
        in.Sigma_Zeta = 0.20;
	in.which_var = 0;
	double chiU = 1e-10;
	double Up = 1e-10;
	double Ud = 1e-10;
	double Chi = 1e-10;

	int c;
	int printFlag = 0;
	int modeFlag = 0;
	int statsFlag = 0;

	opterr = 0;

	while ((c = getopt (argc, argv, "m:Z:X:c:S:B:I:D:TPFERKOAQ")) != -1)
   	{ 
	switch (c) 
      	{
      		case 'm':
			in.mS = strtof(optarg,NULL);
        		break;
      		case 'Z':
			in.mZprime = strtof(optarg,NULL);
        		break;
      		case 'c':
        		in.eCut = 0.14;
        		in.thCut = strtof(optarg,NULL);
        		break;
       
      		case 'R':
        		in.eFloor = 0.02; //100 MeV
        		in.eRatio = 0.2; // Lowest is less than 10% of Highest now 20%
        		break;
		case 'E':
			modeFlag = 1;
			break;
		case 'A':
			modeFlag = 2;
			break;
    		case 'P':
        		printFlag = 1;
			break;
    		case 'F':
        		modeFlag = 3;
			break;
    		case 'O':
        		modeFlag = 4;
			break;
      		case 'X':
			//Up =  strtof(optarg,NULL);
			//Ud =  strtof(optarg,NULL);
			//Chi =  strtof(optarg,NULL);
			if (sscanf(optarg, "%lf:%lf:%lf", &Up,&Ud,&Chi) != 2)
			Up=pow(10,Up);
			Ud=pow(10,Ud);
			Chi=pow(10,Chi);
			//istd::cout<<"1 :"<<one<<" 2: "<<two<<" 3: "<<three<<std::endl;
			break;
    		case 'B':
			if(!strcmp(optarg,"E")){ modeFlag = 5; 	     in.which_var = ENERGY_FLAG; }
			else if(!strcmp(optarg,"A")){ modeFlag = 6;  in.which_var = ANGULAR_FLAG;  }
			else if(!strcmp(optarg,"Q")){ modeFlag = 9;  in.which_var= QE_FLAG; }
			else { printf("Very bad thing!\nAborting...\n\nYou're probably not using the -B flag correctly.\n\n"); exit(1); }
			break;
		case 'I':
			modeFlag=10;
			if(!strcmp(optarg,"E")){      in.which_var = ENERGY_FLAG; }
			else if(!strcmp(optarg,"A")){ in.which_var = ANGULAR_FLAG; }
			else if(!strcmp(optarg,"Q")){ in.which_var= QE_FLAG; }
			else { printf("Very bad thing!\nAborting...\n\nYou're probably not using the -I flag correctly.\n\n"); exit(1); }
			break;
		case 'S':
                        statsFlag = 1;
			in.Sigma_Zeta =  strtof(optarg,NULL);
			break;
		case 'T':
			modeFlag=7;
			break;
		case 'Q':
			modeFlag=8;
			break;
		case 'K':
			modeFlag=11;
			break;
		case 'D':
			in.pos_selector=strtof(optarg,NULL);
			break;
      		case '?':
			printf("Abandon hope all ye who enter this value: %c\n",optopt);
			printf("Allowed arguments:\n"
				"\t-m\tSets mS mass.\n"
				"\t-Z\tSets mZprime mass.\n"
				"\t-X\tSets Up, Ud, Chi in the form ``-X up:ud:chi``\n"
				"\t-c\tSets cuts (0.14,5.0).\n"
				"\t-P\tPrints input parameters.\n"
				"\t-E(-A)(-Q)\tPrints the (E)nergy, (A)ngular spectrum or (Q)e spectrum.\n"
				"\t-I? \tIndividual run of either the (IE)nergy, (IA)ngular or (IQ)e minimizer (NEED STATS ON).\n"
				"\t-B? \t'Brints' the (BE)nergy, (BA)ngular spectrum or (BQ)e spectrum.\n"
				"\t-F\tPrints mimima.\n"
				"\t-R\tTurns on energy asymmetry cut.\n"
				"\t-O\tPrints the cut efficiency.\n"
				"\t-S\tToggles Systematics of bkg (Default off).\n"
				"\t-K\tPLots the (K)ontainment efficiency. (Needs -X to define chiU.)\n"
				"\t-T\tTesting ground, god knows what you will find.\n"
				"\t-D\t D? anti Deutrino Mode flag\n");
                  	return 1;
      		default:
			printf("I don't know how you got here.\n");
        		abort ();
 	}}

	//If option -P is given, we print the input information.		
	if(printFlag==1)
	{ 	
		if(!modeFlag){
			if(!statsFlag){
				printf("# mS = %.5lf\n# mZprime = %.5lf\n# eCut = %.5lf\n# thCut = %.5lf\n# eFloor = %.5lf\n# eRatio = %.5lf\n",in.mS, in.mZprime, in.eCut, in.thCut, in.eFloor, in.eRatio);}
			else{
				printf("# mS = %.5lf\n# mZprime = %.5lf\n# eCut = %.5lf\n# thCut = %.5lf\n# eFloor = %.5lf\n# eRatio = %.5lf\n# sigma_zeta = %.5lf\n",in.mS, in.mZprime, in.eCut, in.thCut, in.eFloor, in.eRatio, in.Sigma_Zeta);}
		}
		else 
		{
		printf("# mS = %.5lf\n# mZprime = %.5lf\n# eCut = %.5lf\n# thCut = %.5lf\n# eFloor = %.5lf\n# eRatio = %.5lf\n# chiU = %.5g\n# Up = %.5g\n# Ud = %.5g\n# Chi = %.5g\n",in.mS, in.mZprime, in.eCut, in.thCut,in.eFloor,in.eRatio,chiU,Up,Ud,Chi);
		}
	}
	//Set up an empty event array and populate it from the decay files.
	static double events[NUMEVENTS][NUM_EVENT_OBS];
	wipeEventArray(events); //wipe the array.
	getEvents_NU(in,events); // populates event array. 

	static double events_NU[NUMEVENTS][NUM_EVENT_OBS];
	wipeEventArray(events_NU);
	getEvents_NU(in,events_NU);

	static double events_NUBAR[NUMEVENTS][NUM_EVENT_OBS];
	wipeEventArray(events_NUBAR);
	getEvents_NUBAR(in,events_NUBAR);


	double cutEfficiency = applyObservableCuts(in,events); //applys cuts. (Post-cut events may have lots of blank slots).


	double cutEfficiency_NU = applyObservableCuts(in,events_NU); //applys cuts. (Post-cut events may have lots of blank slots).
	double cutEfficiency_NUBAR = applyObservableCuts(in,events_NUBAR); //applys cuts. (Post-cut events may have lots of blank slots).

	if(modeFlag)
	{


		//std::cout<<"# Start Mode Flag"<<std::endl;
		int i;
		double eGram[EBINS];
		double cosGram[COSBINS];
		double qeGram[QEBINS];
		double contEfficiency;
		double finalScale = getTotalNumEvents_NU(in);
	
		std::vector<double> bfspectrum;
		double bfzeta=0;

		//A slightly hacky way to print the best-fit spectra.
		if(modeFlag == 5 || modeFlag == 6 || modeFlag == 9)// 5: print best fitting energy spec. 6: print best fitting angular spec. 8: print best fitting QE spectrum
		{
			BF_RESULT * bestfit = new BF_RESULT();	
			bestfit = mcmc_stats_fit_spectra_indiv(in,cutEfficiency,events);


			switch(modeFlag){
				case 5:
					Up = bestfit->E_bf_Up;
					Ud = bestfit->E_bf_Ud;
					Chi = bestfit->E_bf_Chi;
					modeFlag = 1; statsFlag=0;
					bfspectrum = bestfit->E_spectrum;
					bfzeta=bestfit->stats_bf;
					//contEfficiency   = histogrammer_indiv2(in,Up,Ud,Chi,cutEfficiency,events,eGram,ENERGY_FLAG,finalScale);
					//std::cout<<"# zeta_bkg "<<bfzeta<< " chi: "<<pow((bfzeta/0.03),2)<<std::endl;

				
					break;
				case 6:
					Up = bestfit->A_bf_Up;
					Ud = bestfit->A_bf_Ud;
					Chi = bestfit->A_bf_Chi;
					modeFlag=2; statsFlag=0;
					bfspectrum = bestfit->A_spectrum;
					//contEfficiency   = histogrammer_indiv2(in,Up,Ud,Chi,cutEfficiency,events,cosGram,ANGULAR_FLAG,finalScale);
					break;
				case 9:
					Up = bestfit->QE_bf_Up;
					Ud = bestfit->QE_bf_Ud;
					Chi = bestfit->QE_bf_Chi;
					modeFlag=8; statsFlag=0;
					bfspectrum = bestfit->QE_spectrum;
					//contEfficiency   = histogrammer_indiv2(in,Up,Ud,Chi,cutEfficiency,events,qeGram,QE_FLAG,finalScale);
					break;

				default:							
					std::cout<<"Error: Its definitely impossible to be in this switch case."<<std::endl;
			}

			


		}
		

		// The main modeFlag fractionating column.
		if(modeFlag == 1)
		{
			//std::cout<<"# zeta_bkg "<<bfzeta<< " chi: "<<pow((bfzeta/0.03),2)<<std::endl;
			//contEfficiency   = histogrammer_indiv2(in,Up,Ud,Chi,cutEfficiency,events,eGram,0,finalScale);
			for(i=0;i<EBINS;i++)
			{
				std::cout<<BintoCentralE(i)<<" "<<bfspectrum[i]<<std::endl;
				//std::cout<<BintoCentralE(i)<<" "<<eGram[i]<<std::endl;
			}
		}
		else if(modeFlag == 2)
		{ 
			//contEfficiency   = histogrammer_indiv2(in,Up,Ud,Chi,cutEfficiency,events,cosGram,1,finalScale);
			for(i=0;i<COSBINS;i++)
			{
				std::cout<<BintoCentralCos(i)<<" "<<bfspectrum[i]<<std::endl;
			}
		}
		else if(modeFlag == 8)
		{ 
			//contEfficiency   = histogrammer_indiv2(in,Up,Ud,Chi,cutEfficiency,events,qeGram,2,finalScale);
			for(i=0;i<QEBINS;i++)
			{
				std::cout<<BintoCentalQE(i)<<" "<<bfspectrum[i]<<std::endl;
			}
		}
		else if(modeFlag == 3)
		{ 
		std::cout<<" Removed plot_minimization_spectrum here"<<std::endl;
		//	plot_minimization_spectrum(in,cutEfficiency,events);
		}
		else if(modeFlag == 4)
		{ 
			std::cout<<in.mS<<" "<<in.mZprime<<" "<<cutEfficiency<<std::endl;
		}
		else if(modeFlag == 7)
		{
			//std::cout<<"Test:"<<QEfromEandCos(0.2,0.9)<<std::endl;
			//for(double ee=0.1; ee<3; ee+=0.1){
			//	std::cout<<ee<<" "<<QEfromEandCos(ee,0)<<"  "<<QEtoBin(QEfromEandCos(ee,0))<<std::endl;
			//}
			//
			//stats_fit_spectra_indiv(in,cutEfficiency,events,1);
			//for(double k =0.001; k<=0.2; k+=0.001){
			//	std::cout<<k<<" "<<boundUpeaky(k,0.1)<<std::endl;
			//}


		} else if (modeFlag == 10) {
			
		//	mcmc_stats_fit_spectra_indiv(in,cutEfficiency,events);
			mcmc_stats_dual(in,cutEfficiency_NU, cutEfficiency_NUBAR,events_NU,events_NUBAR);
			//plot_surface(in,cutEfficiency,events);

		} else if (modeFlag == 11) {
			
			std::cout<<in.mS<<" "<<in.mZprime<<" "<<chiU<<" "<<contEfficiency<<std::endl;

		} else { std::cout<<"BAD THING!"<<std::endl; }
	}
	else 
	{
		if(statsFlag)
		{
		
			if( modeFlag == 1){    in.which_var = ENERGY_FLAG; }
			if( modeFlag == 2){  in.which_var = ANGULAR_FLAG;  }
			if( modeFlag == 8){  in.which_var= QE_FLAG; }
			in.which_var=1;	

			double delchi = spectra_indiv(in,Up,Ud,Chi);
			std::cout<<"# ms: "<<in.mS<<" mz: "<<in.mZprime<<" Up: "<<Up<<" Ud: "<<Ud<<" Chi: "<<Chi<<" delchi: "<<delchi<<std::endl;
			std::cout<<in.mS<<" "<<in.mZprime<<" "<<Up<<" "<<Ud<<" "<<Chi<<" "<<delchi<<std::endl;
		}
		else
		{
			//ifit_spectra(in,cutEfficiency,events);
			//std::cout<<"# Begin End of chain,  fit_spectra"<<std::endl;
			std::cout<<"Removed  fit_spectra"<<std::endl;
		}	
	}

return 0;
}
