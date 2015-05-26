#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include <unistd.h>

#include "LR.h"

#define SIG_ON 1.0
#define SIG_OFF 0.0

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


double plot_minimization_spectrum(CL_input in, double cutEfficiency, double events[NUMEVENTS][NUM_EVENT_OBS])
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

double fit_spectra(CL_input in, double cutEfficiency, double events[NUMEVENTS][NUM_EVENT_OBS]){
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
	//while(E_sum < 80 && A_sum < 80 && logchiU < log(sqrt(boundU(temp_mS)*boundUtau(temp_mS))*sqrt(boundChi(temp_mZprime)))/log(10.0)) 
	while(E_sum < 80 && A_sum < 80 && logchiU < log(sqrt(boundU(temp_mS))*sqrt(boundChi(temp_mZprime)))/log(10.0)) 
	//while(E_sum < 80 && A_sum < 80) 
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

double stats_fit_spectra(CL_input in, double cutEfficiency, double events[NUMEVENTS][NUM_EVENT_OBS])
{
	double Ezeta_b = 0.0;
	double Azeta_b = 0.0;
	double sigma_s = 1.0;
	double sigma_zeta = in.Sigma_Zeta;

	
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
//#########################################################################

        std::vector<double > eVGram(eGram, eGram + EBINS);
        std::vector<double > cosVGram(cosGram, cosGram + COSBINS);
        
        std::vector<double > bf_zeta_b;  
        double bf_chi ;
        std::vector<double > abf_zeta_b; 
        double abf_chi ;
//      marginalize(&bf_zeta_b,&bf_chi, eVGram);
//      std::cout<<"Bf: "<<bf_zeta_b[0]<<" BFChi: "<<bf_chi<<std::endl;

        
// ########################   Calculate LogLiki (Dchi^2) for Background only with SYtematics ################

        std::vector<double > eZeros(EBINS, 0.0);
        std::vector<double > EBF_bkg_only_zeta_b = {0};
        double EBF_bkg_only_chi = 10000;
        nuisMarginalize(& EBF_bkg_only_zeta_b,&EBF_bkg_only_chi, &eZeros,0, sigma_zeta);
        std::cout<<"# Energy sigma_zeta Bf: "<< EBF_bkg_only_zeta_b[0]<<" BF Chi: "<<EBF_bkg_only_chi<<std::endl;

        std::vector<double > aZeros(COSBINS, 0.0);
        std::vector<double > ABF_bkg_only_zeta_b = {0};
        double ABF_bkg_only_chi = 10000;
        nuisMarginalize(& ABF_bkg_only_zeta_b,&ABF_bkg_only_chi, &aZeros,1,sigma_zeta);
        std::cout<<"# Angle sigma_zeta Bf: "<<ABF_bkg_only_zeta_b[0]<<" BFChi: "<<ABF_bkg_only_chi<<std::endl;


//#########################################################################

	logchiU=logchiU_start;
	//while(E_sum < 80 && A_sum < 80 && logchiU < log(sqrt(boundU(temp_mS)*boundUtau(temp_mS))*sqrt(boundChi(temp_mZprime)))/log(10.0)) 
	while(E_sum < 80 && A_sum < 80 && logchiU < log(sqrt(boundU(temp_mS))*sqrt(boundChi(temp_mZprime)))/log(10.0)) 
	//while(E_sum < 80 && A_sum < 80) 
	{

		in.mS = temp_mS;
		in.mZprime = temp_mZprime;

		chiU=pow(10.0,logchiU);
		E_sum=0.0;
		A_sum=0.0;

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram);
                        eVGram.assign(eGram, eGram+EBINS);
                        nuisMarginalize(&bf_zeta_b, &bf_chi, &eVGram,0,sigma_zeta);
                        Ezeta_b=bf_zeta_b[0];
                
                        cosVGram.assign(cosGram, cosGram+COSBINS);
                        nuisMarginalize(&abf_zeta_b, &abf_chi, &cosVGram,1,sigma_zeta);
                        Azeta_b=abf_zeta_b[0];

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
				temp_bg = (1.0+Ezeta_b)*eB[bin];
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
				temp_bg = (1.0+Azeta_b)*cosB[bin];
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
	//	double sigma_zeta = 0.05;
		// add prior on zeta.
		E_sum += pow((Ezeta_b/sigma_zeta),2.0);
		A_sum += pow((Azeta_b/sigma_zeta),2.0);
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
                        eVGram.assign(eGram, eGram+EBINS);
                        nuisMarginalize(&bf_zeta_b, &bf_chi, &eVGram,0,sigma_zeta);
                        Ezeta_b=bf_zeta_b[0];

		E_N_events = 0;
		E_N_sig_events = 0;
		E_N_bg_events = 0;

		int bin = 0;
		double temp_sig,temp_bg;
	//	if (sigma_s > 1e-5){
			for(bin=0;bin<EBINS;bin++)
			{
				temp_sig = sigma_s*eGram[bin];
				temp_bg = (1.0+Ezeta_b)*eB[bin];
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
		//double sigma_zeta = 0.05;
		// add prior on zeta.
		E_sum += pow((Ezeta_b/sigma_zeta),2.0);
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
                        cosVGram.assign(cosGram, cosGram+COSBINS);
                        nuisMarginalize(&abf_zeta_b, &abf_chi, &cosVGram,1,sigma_zeta);
                        Azeta_b = abf_zeta_b[0];

		A_N_sig_events = 0;
		A_N_bg_events = 0;
		A_N_events = 0;

		int bin = 0;
		double temp_sig,temp_bg;
	//	if (sigma_s > 1e-5){
			for(bin=0;bin<COSBINS;bin++)
			{
				temp_sig = sigma_s*cosGram[bin];
				temp_bg = (1.0+Azeta_b)*cosB[bin];
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
	//	double sigma_zeta = 0.05;
		// add prior on zeta.
		A_sum += pow((Azeta_b/sigma_zeta),2.0);
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


//	if(fabs(E_best-65.1554) < 1e-3) { E_best=65.1554; }
//	if(fabs(A_best-38.9753) < 1e-3) { A_best=38.9753; }
        if(fabs(E_best-EBF_bkg_only_chi) < 1e-3) { E_best=EBF_bkg_only_chi; }
        if(fabs(A_best-ABF_bkg_only_chi) < 1e-3) { A_best=ABF_bkg_only_chi; }


//	double E_check = E_best_chiU*E_best_chiU*best_E_contEfficiency*getTotalNumEvents(in)*cutEfficiency;
//	double A_check = A_best_chiU*A_best_chiU*best_A_contEfficiency*getTotalNumEvents(in)*cutEfficiency;

//	std::cout<<E_check<<" "<<A_check<<std::endl;
	std::cout<<"# From Mark: "<<getTotalNumEvents(in)<<" Energy: "<<best_E_N_events<<"= ("<<E_best_N_sig_events<<" + "<<E_best_N_bg_events<<") Angle: "<<best_A_N_events<<"= ("<<A_best_N_sig_events<<" + "<<A_best_N_bg_events<<") "<<std::endl;

//	std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<65.1554-E_best<<" "<<E_best_chiU<<" "<<38.9753-A_best<<" "<<A_best_chiU<<std::endl;
//std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<65.1554-E_best<<" "<<E_best_chiU<<" "<<best_E_contEfficiency<<" "<<39.9753-A_best<<" "<<A_best_chiU<<" "<<best_A_contEfficiency<<" "<<cutEfficiency<<std::endl;
        std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<EBF_bkg_only_chi-E_best<<" "<<E_best_chiU<<" "<<best_E_contEfficiency<<" "<<ABF_bkg_only_chi-A_best<<" "<<A_best_chiU<<" "<<best_A_contEfficiency<<" "<<cutEfficiency<<std::endl;

return 0.0;
}

int main(int argc, char * argv[])
{
	static CL_input in;
	in.eCut = 0.0;
	in.thCut = 180.0;
	in.eFloor = 0.0; // I believe 0.0 for Floor and Ratio mean that the cut is never passed.
	in.eRatio = 0.0;
        in.Sigma_Zeta = 0.20;
	double chiU = 1e-10;

	int c;
	int printFlag = 0;
	int modeFlag = 0;
	int statsFlag = 0;

	opterr = 0;

	while ((c = getopt (argc, argv, "m:Z:X:c:S:PFEROA")) != -1)
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
        		in.eFloor = 0.1; //100 MeV
        		in.eRatio = 0.1; // Lowest is less than 10% of Highest
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
    		case 'O':
        		modeFlag = 4;
			break;
		case 'S':
                        statsFlag = 1;
			in.Sigma_Zeta =  strtof(optarg,NULL);
			break;
      		case '?':
			printf("Abandon hope all ye who enter this value: %c\n",optopt);
			printf("Allowed arguments:\n\t-m\tsets mS mass.\n\t-Z\tsets mZprime mass.\n\t-X\tsets chiU\n\t-c\tsets cuts (0.14,5.0).\n\t-P\tprints input parameters.\n\t-E(-A)\tprints the (E)nergy or (A)ngular spectrum.\n\t-F\tprints mimima\n\t-R\tturns on energy asymmetry cut\n\t-O\tprints the cut efficiency.\n\t-S\t Toggles Systematics of bkg (Default off)\n");
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
		printf("# mS = %.5lf\n# mZprime = %.5lf\n# eCut = %.5lf\n# thCut = %.5lf\n# eFloor = %.5lf\n# eRatio = %.5lf\n# chiU = %.5g\n",in.mS, in.mZprime, in.eCut, in.thCut,in.eFloor,in.eRatio,chiU);
		}
	}

	//Set up an empty event array and populate it from the decay files.
	static double events[NUMEVENTS][NUM_EVENT_OBS];
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
		else if(modeFlag == 4)
		{ 
			std::cout<<in.mS<<" "<<in.mZprime<<" "<<cutEfficiency<<std::endl;
		}
		else { std::cout<<"BAD THING!"<<std::endl; }
	}
	else 
	{
		if(statsFlag)
		{
			stats_fit_spectra(in,cutEfficiency,events);
		}
		else
		{
			fit_spectra(in,cutEfficiency,events);
		}	
	}

return 0;
}
