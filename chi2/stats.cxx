#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <unistd.h>

#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/VariableMetricMinimizer.h>

#include "LR.h"

#define SIG_ON 1.0
#define SIG_OFF 0.0

using namespace ROOT::Minuit2;

double getEvents(CL_input input, double events[][4])
{
	double mS = input.mS;
	double mZprime = input.mZprime;

	FILE *ptr_file;
    	
	char buf[3000];

	char * pch;

	int n = 1;
	int m = 0;
	char s[100];
	char filename[500] = "../decay/data/\0";
	sprintf(s,"%.4lf_%.4lf.dat", mS, mZprime);
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

int wipeEventArray(double events[][4])
{
	int i;
	for (i=0;i<=NUMEVENTS-1;i++)
	{ 
		events[i][0]=0.0;
		events[i][1]=0.0;
		events[i][2]=0.0;
		events[i][3]=0.0;
	}
return 0;
}

double applyObservableCuts(CL_input input, double events[][4])
{
	double ECut = input.eCut;
	double thCut = input.thCut;

	int i,m;
	int good=0;

	double goodEvents[NUMEVENTS][4];
	wipeEventArray(goodEvents);

	for (i=0;i<=NUMEVENTS-1;i++)
	{ 
		if(events[i][0]>ECut) //E_sum > ECut	
		{ 
			if(events[i][2]<thCut) // Ang.Sep. < thCut
			{ 	
				for(m=0;m<4;m++){goodEvents[good][m]=events[i][m];} 
				good++;
			}
		}
	}

	for(i=0;i<=NUMEVENTS-1;i++)
	{
		for(m=0;m<4;m++){events[i][m]=goodEvents[i][m];} 
	}

return ((double)good)/((double)NUMEVENTS);
}


double exclude_point(CL_input in, double sigma_s)
{
	double events[NUMEVENTS][4];
	wipeEventArray(events); //wipe the array.
	getEvents(in,events); // populates event array. 

	double cutEfficiency = applyObservableCuts(in,events); //applys cuts; post-cut events have lots of blank slots.

	double min_signal = 1e5;
	double min_background_only = 1e5;

	{
		MnUserParameters upar;

		upar.Add("mS",in.mS);
		upar.Add("mZprime",in.mZprime);
		upar.Add("logchiU",-4.0,0.1,-7.0,-1.0);
		upar.Add("zeta_b",0.0,0.05,-2.0,2.0);
		upar.Add("sigma_s",sigma_s);
		upar.Add("cutEff",cutEfficiency);

		E_log_likelihood fcn(events);

		MnMigrad migrad(fcn,upar);

		migrad.Fix("mS");
		migrad.Fix("mZprime");
		migrad.Fix("logchiU");
//		migrad.Fix("zeta_b");
		migrad.Fix("sigma_s");
		migrad.Fix("cutEff");
	
		FunctionMinimum min = migrad();

		min_signal = min.Fval();

		std::cout<<min.Fval()<<std::endl;

	}

return min_signal;
}
double plot_minimization_spectrum(CL_input in, double cutEfficiency, double events[NUMEVENTS][4])
{
	double zeta_b = 0.0;
	double sigma_s = 1.0;
	
	double N=0.0;
	double lambda=0.0;
	double E_sum=0.0;
	double A_sum=0.0;
	
	double eGram[EBINS];
	double cosGram[COSBINS];

	int i;
	for(i=0;i<EBINS;i++)
	{
		eGram[i]=0.0;
	}
	for(i=0;i<COSBINS;i++)
	{
		cosGram[i]=0.0;
	}

	double eO[EBINS];
	double eB[EBINS];

		//Ebinlow(GeV) 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9
		eO[0] = 204.0;
		eO[1] = 280.0; 
		eO[2] = 214.0; 
		eO[3] = 99.0; 
		eO[4] = 83.0; 
		eO[5] = 59.0; 
		eO[6] = 51.0; 
		eO[7] = 33.0; 
		eO[8] = 37.0; 
		eO[9] = 23.0; 	
		eO[10] = 19.0;
 		eO[11] = 21.0; 
		eO[12] = 12.0; 
		eO[13] = 16.0; 		
		eO[14] = 4.0; 
		eO[15] = 9.0; 
		eO[16] = 4.0; 
		eO[17] = 7.0; 
		eO[18] = 3.0;

		eB[0] = 151.5;
		eB[1] = 218.8;
		eB[2] = 155.6;
		eB[3] = 108.7;
		eB[4] = 72.5;
		eB[5] = 57.6;
		eB[6] = 45;
		eB[7] = 38.5;
		eB[8] = 31.4;
		eB[9] = 22.2;
		eB[10] = 20.4;
		eB[11] = 17.2;
		eB[12] = 14.1;
		eB[13] = 10.2;
		eB[14] = 9.1;
		eB[15] = 8.2;
		eB[16] = 5.6;
		eB[17] = 5.7;
		eB[18] = 2.9;

	double cosO[COSBINS];
	double cosB[COSBINS];

	cosO[0] = 22;
	cosO[1] = 34;
	cosO[2] = 43;
	cosO[3] = 41;
	cosO[4] = 60;
	cosO[5] = 87;
	cosO[6] = 90;
	cosO[7] = 139;
	cosO[8] = 237;
	cosO[9] = 429;

	cosB[0] = 19.9;
	cosB[1] = 23.1;
	cosB[2] = 28.8;
	cosB[3] = 32.1;
	cosB[4] = 46.4;
	cosB[5] = 63.1;
	cosB[6] = 86.1;
	cosB[7] = 121;
	cosB[8] = 196.8;
	cosB[9] = 390;

	double temp_tot_E =0.0;
	double temp_tot_A =0.0;
//	int j;
//	for(j=0;j<COSBINS;j++)
//	{	
//		temp_tot_A += cosO[j] - cosB[j];
//	} // temp_tot_A = 174.7
//	temp_tot_A = 174.7;
//	for(j=0;j<EBINS;j++)
//	{	
//		temp_tot_E += eO[j] - eB[j];
//	} // temp_to_E = 182.8
//	temp_tot_E = 182.8;


//	std::cout<<"# Ob. Excess Tot: (Energy: "<<temp_tot_E<<") (Angular: "<<temp_tot_A<<")"<<std::endl;


	double logchiU,chiU,contEfficiency;
	double best_E_contEfficiency=1e50;
	double best_A_contEfficiency=1e50;
	double best_E_N_events=1e-50;
	double best_A_N_events=1e-50;
	double E_best_chiU = 1e50;
	double E_best = 1e4;
	double A_best_chiU = 1e50;
	double A_best = 1e4;
	double temp_mS = in.mS;
	double temp_mZprime = in.mZprime;

	double E_N_events=0;
	double A_N_events=0;
	double E_N_sig_events=0;
	double A_N_sig_events=0;
	double E_N_bg_events=0;
	double A_N_bg_events=0;
	double E_best_N_sig_events = 0;
	double A_best_N_sig_events = 0;
	double E_best_N_bg_events = 0;
	double A_best_N_bg_events = 0;

	double logchiU_step = 0.01;
	double logchiU_start = -7.0;
//int loopflag = 1;
//while(loopflag)
//{
	logchiU=logchiU_start;
	while(E_sum < 1e3 && A_sum < 1e3)
	{

		in.mS = temp_mS;
		in.mZprime = temp_mZprime;

		chiU=pow(10.0,logchiU);
		E_sum=0.0;
		A_sum=0.0;

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram);

		E_N_events = 0;
		E_N_sig_events = 0;
		E_N_bg_events = 0;
		A_N_sig_events = 0;
		A_N_bg_events = 0;
		A_N_events = 0;

		int bin = 0;
		double temp_sig,temp_bg;
	//	if (sigma_s > 1e-5){
			for(bin=0;bin<EBINS;bin++)
			{
				temp_sig = sigma_s*eGram[bin];
				temp_bg = (1.0+zeta_b)*eB[bin];
				lambda = temp_sig + temp_bg;
				N = eO[bin]; //MB has seen O[] events.
				
				E_N_events += lambda;
				E_N_sig_events += temp_sig;
				E_N_bg_events += temp_bg;
			//	E_sum+= (lambda-N)*(lambda-N)/lambda;
				E_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
			}

			for(bin=0;bin<COSBINS;bin++)
			{
				temp_sig = sigma_s*cosGram[bin];
				temp_bg = (1.0+zeta_b)*cosB[bin];
				lambda = temp_sig + temp_bg;
				N = cosO[bin]; //MB has seen O[] events.
				
				A_N_events += lambda;
				A_N_sig_events += temp_sig;
				A_N_bg_events += temp_bg;
			//	A_sum+= (lambda-N)*(lambda-N)/lambda;
				A_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
			}
	//	}

		//a guesstimate of the systematic error on the background.
		double sigma_zeta = 0.05;
		// add prior on zeta.
		E_sum += pow((zeta_b/sigma_zeta),2.0);
		A_sum += pow((zeta_b/sigma_zeta),2.0);
		std::cout<<chiU<<" "<<E_sum<<" "<<A_sum<<std::endl;

		logchiU += logchiU_step;

	}

return 0.0;
}


double fit_E_spectrum(CL_input in, double cutEfficiency, double events[NUMEVENTS][4])
{
	double zeta_b = 0.0;
	double sigma_s = 1.0;
	
	double N=0.0;
	double lambda=0.0;
	double E_sum=0.0;
	double A_sum=0.0;
	
	double eGram[EBINS];
	double cosGram[COSBINS];

	int i;
	for(i=0;i<EBINS;i++)
	{
		eGram[i]=0.0;
	}
	for(i=0;i<COSBINS;i++)
	{
		cosGram[i]=0.0;
	}

	double eO[EBINS];
	double eB[EBINS];

		//Ebinlow(GeV) 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9
		eO[0] = 204.0;
		eO[1] = 280.0; 
		eO[2] = 214.0; 
		eO[3] = 99.0; 
		eO[4] = 83.0; 
		eO[5] = 59.0; 
		eO[6] = 51.0; 
		eO[7] = 33.0; 
		eO[8] = 37.0; 
		eO[9] = 23.0; 	
		eO[10] = 19.0;
 		eO[11] = 21.0; 
		eO[12] = 12.0; 
		eO[13] = 16.0; 		
		eO[14] = 4.0; 
		eO[15] = 9.0; 
		eO[16] = 4.0; 
		eO[17] = 7.0; 
		eO[18] = 3.0;

		eB[0] = 151.5;
		eB[1] = 218.8;
		eB[2] = 155.6;
		eB[3] = 108.7;
		eB[4] = 72.5;
		eB[5] = 57.6;
		eB[6] = 45;
		eB[7] = 38.5;
		eB[8] = 31.4;
		eB[9] = 22.2;
		eB[10] = 20.4;
		eB[11] = 17.2;
		eB[12] = 14.1;
		eB[13] = 10.2;
		eB[14] = 9.1;
		eB[15] = 8.2;
		eB[16] = 5.6;
		eB[17] = 5.7;
		eB[18] = 2.9;

	double cosO[COSBINS];
	double cosB[COSBINS];

	cosO[0] = 22;
	cosO[1] = 34;
	cosO[2] = 43;
	cosO[3] = 41;
	cosO[4] = 60;
	cosO[5] = 87;
	cosO[6] = 90;
	cosO[7] = 139;
	cosO[8] = 237;
	cosO[9] = 429;

	cosB[0] = 19.9;
	cosB[1] = 23.1;
	cosB[2] = 28.8;
	cosB[3] = 32.1;
	cosB[4] = 46.4;
	cosB[5] = 63.1;
	cosB[6] = 86.1;
	cosB[7] = 121;
	cosB[8] = 196.8;
	cosB[9] = 390;

	double temp_tot_E =0.0;
	double temp_tot_A =0.0;
//	int j;
//	for(j=0;j<COSBINS;j++)
//	{	
//		temp_tot_A += cosO[j] - cosB[j];
//	} // temp_tot_A = 174.7
//	temp_tot_A = 174.7;
//	for(j=0;j<EBINS;j++)
//	{	
//		temp_tot_E += eO[j] - eB[j];
//	} // temp_to_E = 182.8
//	temp_tot_E = 182.8;


//	std::cout<<"# Ob. Excess Tot: (Energy: "<<temp_tot_E<<") (Angular: "<<temp_tot_A<<")"<<std::endl;


	double logchiU,chiU,contEfficiency;
	double best_E_contEfficiency=1e50;
	double best_A_contEfficiency=1e50;
	double best_E_N_events=1e-50;
	double best_A_N_events=1e-50;
	double E_best_chiU = 1e50;
	double E_best = 1e4;
	double A_best_chiU = 1e50;
	double A_best = 1e4;
	double temp_mS = in.mS;
	double temp_mZprime = in.mZprime;

	double E_N_events=0;
	double A_N_events=0;
	double E_N_sig_events=0;
	double A_N_sig_events=0;
	double E_N_bg_events=0;
	double A_N_bg_events=0;
	double E_best_N_sig_events = 0;
	double A_best_N_sig_events = 0;
	double E_best_N_bg_events = 0;
	double A_best_N_bg_events = 0;

	double logchiU_step = 0.1;
	double logchiU_start = -7.0;
//int loopflag = 1;
//while(loopflag)
//{
	logchiU=logchiU_start;
	while(E_sum < 80 && A_sum < 80)
	{

		in.mS = temp_mS;
		in.mZprime = temp_mZprime;

		chiU=pow(10.0,logchiU);
		E_sum=0.0;
		A_sum=0.0;

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram);

		E_N_events = 0;
		E_N_sig_events = 0;
		E_N_bg_events = 0;
		A_N_sig_events = 0;
		A_N_bg_events = 0;
		A_N_events = 0;

		int bin = 0;
		double temp_sig,temp_bg;
	//	if (sigma_s > 1e-5){
			for(bin=0;bin<EBINS;bin++)
			{
				temp_sig = sigma_s*eGram[bin];
				temp_bg = (1.0+zeta_b)*eB[bin];
				lambda = temp_sig + temp_bg;
				N = eO[bin]; //MB has seen O[] events.
				
				E_N_events += lambda;
				E_N_sig_events += temp_sig;
				E_N_bg_events += temp_bg;
			//	E_sum+= (lambda-N)*(lambda-N)/lambda;
				E_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
			}

			for(bin=0;bin<COSBINS;bin++)
			{
				temp_sig = sigma_s*cosGram[bin];
				temp_bg = (1.0+zeta_b)*cosB[bin];
				lambda = temp_sig + temp_bg;
				N = cosO[bin]; //MB has seen O[] events.
				
				A_N_events += lambda;
				A_N_sig_events += temp_sig;
				A_N_bg_events += temp_bg;
			//	A_sum+= (lambda-N)*(lambda-N)/lambda;
				A_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
			}
	//	}

		//a guesstimate of the systematic error on the background.
		double sigma_zeta = 0.05;
		// add prior on zeta.
		E_sum += pow((zeta_b/sigma_zeta),2.0);
		A_sum += pow((zeta_b/sigma_zeta),2.0);
//		std::cout<<chiU<<" "<<E_sum<<" "<<A_sum<<std::endl;

		logchiU += logchiU_step;

	}

logchiU -= logchiU_step;
double breakpoint_logchiU = logchiU;
double E_sum_previous = E_sum + 1e-8;
double A_sum_previous = A_sum + 1e-8;
double temp_E_sum = 0.0;
double temp_A_sum = 0.0;
logchiU_step = 0.001;

//printf("I did a thing.\n");
	while( E_sum < E_sum_previous)
	{

		temp_E_sum = E_sum;
	
		logchiU -= logchiU_step;

//printf("I did a thing. Again.\n");
		in.mS = temp_mS;
		in.mZprime = temp_mZprime;

		chiU=pow(10.0,logchiU);
		E_sum=0.0;

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram);

		E_N_events = 0;
		E_N_sig_events = 0;
		E_N_bg_events = 0;

		int bin = 0;
		double temp_sig,temp_bg;
	//	if (sigma_s > 1e-5){
			for(bin=0;bin<EBINS;bin++)
			{
				temp_sig = sigma_s*eGram[bin];
				temp_bg = (1.0+zeta_b)*eB[bin];
				lambda = temp_sig + temp_bg;
				N = eO[bin]; //MB has seen O[] events.
				
				E_N_events += lambda;
				E_N_sig_events += temp_sig;
				E_N_bg_events += temp_bg;
			//	E_sum+= (lambda-N)*(lambda-N)/lambda;
				E_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
			}
	//	}

		//a guesstimate of the systematic error on the background.
		double sigma_zeta = 0.05;
		// add prior on zeta.
		E_sum += pow((zeta_b/sigma_zeta),2.0);
//		std::cout<<chiU<<" "<<E_sum<<" "<<A_sum<<std::endl;

		if(E_sum<E_best)
		{ 
			E_best=E_sum; 
			E_best_chiU = chiU; 
			best_E_contEfficiency = contEfficiency; 
			best_E_N_events = E_N_events;
			E_best_N_sig_events = E_N_sig_events;
			E_best_N_bg_events = E_N_bg_events;
		}
		
		E_sum_previous = temp_E_sum;

	}

	logchiU = breakpoint_logchiU;
	while( A_sum < A_sum_previous)
	{

		temp_A_sum = A_sum;
	
		logchiU -= logchiU_step;

//printf("I did a thing. Again.\n");
		in.mS = temp_mS;
		in.mZprime = temp_mZprime;

		chiU=pow(10.0,logchiU);
		A_sum=0.0;

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram);

		A_N_sig_events = 0;
		A_N_bg_events = 0;
		A_N_events = 0;

		int bin = 0;
		double temp_sig,temp_bg;
	//	if (sigma_s > 1e-5){
			for(bin=0;bin<COSBINS;bin++)
			{
				temp_sig = sigma_s*cosGram[bin];
				temp_bg = (1.0+zeta_b)*cosB[bin];
				lambda = temp_sig + temp_bg;
				N = cosO[bin]; //MB has seen O[] events.
				
				A_N_events += lambda;
				A_N_sig_events += temp_sig;
				A_N_bg_events += temp_bg;
			//	A_sum+= (lambda-N)*(lambda-N)/lambda;
				A_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
			}
	//	}

		//a guesstimate of the systematic error on the background.
		double sigma_zeta = 0.05;
		// add prior on zeta.
		A_sum += pow((zeta_b/sigma_zeta),2.0);
//		std::cout<<chiU<<" "<<E_sum<<" "<<A_sum<<std::endl;

		if(A_sum<A_best)
		{ 
			A_best=A_sum; 
			A_best_chiU = chiU; 
			best_A_contEfficiency = contEfficiency; 
			best_A_N_events = A_N_events;
			A_best_N_sig_events = A_N_sig_events;
			A_best_N_bg_events = A_N_bg_events;

		}
		
		A_sum_previous = temp_A_sum;

	}


	if(fabs(E_best-65.1554) < 1e-3) { E_best=65.1554; }
	if(fabs(A_best-38.9753) < 1e-3) { A_best=38.9753; }

//	double E_check = E_best_chiU*E_best_chiU*best_E_contEfficiency*getTotalNumEvents(in)*cutEfficiency;
//	double A_check = A_best_chiU*A_best_chiU*best_A_contEfficiency*getTotalNumEvents(in)*cutEfficiency;

//	std::cout<<E_check<<" "<<A_check<<std::endl;
	std::cout<<"# From Mark: "<<getTotalNumEvents(in)<<" Energy: "<<best_E_N_events<<"= ("<<E_best_N_sig_events<<" + "<<E_best_N_bg_events<<") Angle: "<<best_A_N_events<<"= ("<<A_best_N_sig_events<<" + "<<A_best_N_bg_events<<") "<<std::endl;

//	std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<65.1554-E_best<<" "<<E_best_chiU<<" "<<38.9753-A_best<<" "<<A_best_chiU<<std::endl;
	std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<65.1554-E_best<<" "<<E_best_chiU<<" "<<best_E_contEfficiency<<" "<<39.9753-A_best<<" "<<A_best_chiU<<" "<<best_A_contEfficiency<<" "<<cutEfficiency<<std::endl;
return 0.0;
}

int main(int argc, char * argv[])
{
	static CL_input in;
	in.eCut = 0.0;
	in.thCut = 180.0;
	double chiU = 1e-10;

	int c;
	int printFlag = 0;
	int modeFlag = 0;

	opterr = 0;

	while ((c = getopt (argc, argv, "m:Z:X:c:PFEA")) != -1)
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
      		case 'X':
        		chiU = strtof(optarg,NULL);
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
      		case '?':
			printf("Abandon hope all ye who enter this value: %c\n",optopt);
			printf("Allowed arguments:\n\t-m\tsets mS mass.\n\t-Z\tsets mZprime mass.\n\t-X\tsets chiU\n\t-c\tsets cuts (0.14,5.0).\n\t-P\tprints input parameters.\n\t-E(-A)\tprints the (E)nergy or (A)ngular spectrum.\n");
                  	return 1;
      		default:
			printf("I don't know how you got here.\n");
        		abort ();
 	}}

	//If option -P is given, we print the input information.		
	if(printFlag==1)
	{ 	
		if(!modeFlag){	
		printf("# mS = %.5lf\n# mZprime = %.5lf\n# eCut = %.5lf\n# thCut = %.5lf\n",in.mS, in.mZprime, in.eCut, in.thCut);
		}
		else 
		{
		printf("# mS = %.5lf\n# mZprime = %.5lf\n# eCut = %.5lf\n# thCut = %.5lf\n# chiU = %.5g\n",in.mS, in.mZprime, in.eCut, in.thCut,chiU);
		}
	}


	//Set up an empty event array and populate it from the decay files.
	static double events[NUMEVENTS][4];
	wipeEventArray(events); //wipe the array.
	getEvents(in,events); // populates event array. 

	double cutEfficiency = applyObservableCuts(in,events); //applys cuts. (Post-cut events may have lots of blank slots).

	if(modeFlag)
	{
		int i;
		double eGram[EBINS];
		double cosGram[COSBINS];

		double contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram);

		if(modeFlag == 1)
		{
			for(i=0;i<EBINS;i++)
			{
				std::cout<<BintoCentralE(i)<<" "<<eGram[i]<<std::endl;
			}
		}
		else if(modeFlag == 2)
		{ 
			for(i=0;i<COSBINS;i++)
			{
				std::cout<<BintoCentralCos(i)<<" "<<cosGram[i]<<std::endl;
			}
		}
		else if(modeFlag == 3)
		{ 
			plot_minimization_spectrum(in,cutEfficiency,events);
		}

		else { std::cout<<"BAD THING!"<<std::endl; }
	}
	else 
	{ 
		fit_E_spectrum(in,cutEfficiency,events);
	}

return 0;
}

