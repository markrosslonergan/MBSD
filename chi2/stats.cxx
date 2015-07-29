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
	sprintf(s,"%.3lf_%.3lf.dat", mS, mZprime);
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
	double QE_sum=0.0;

	double eGram[EBINS];
	double cosGram[COSBINS];
	double qeGram[QEBINS];

	int i;
	for(i=0;i<EBINS;i++)
	{
		eGram[i]=0.0;
	}
	for(i=0;i<COSBINS;i++)
	{
		cosGram[i]=0.0;
	}
	for(i=0;i<QEBINS;i++)
	{
		qeGram[i]=0.0;
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

	double qeO[QEBINS];
	double qeB[QEBINS];
      	qeO[0] =  232;qeO[1] = 156 ;qeO[2] = 156 ;qeO[3] = 79 ;qeO[4] = 81 ;qeO[5] = 70 ;qeO[6] = 63 ;qeO[7] = 65 ;qeO[8] = 62 ;qeO[9] = 34 ;qeO[10] = 70;
        qeB[0] = 181.1 ;qeB[1] = 108.4;qeB[2] = 120.4;qeB[3] = 64.2;qeB[4] = 90.3; qeB[5] = 67.7;qeB[6] = 70.4;qeB[7] = 57.5;qeB[8] = 52.3; qeB[9] = 39;qeB[10] = 70.2;



	double temp_tot_E =0.0;
	double temp_tot_A =0.0;
	double temp_tot_QE=0.0;
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
	double best_QE_contEfficiency=1e50;
	double best_E_N_events=1e-50;
	double best_A_N_events=1e-50;
	double best_QE_N_events=1e-50;
	double E_best_chiU = 1e50;
	double E_best = 1e4;
	double A_best_chiU = 1e50;
	double A_best = 1e4;
	double QE_best_chiU = 1e50;
	double QE_best = 1e6;
	double temp_mS = in.mS;
	double temp_mZprime = in.mZprime;

	double E_N_events=0;
	double A_N_events=0;
	double QE_N_events=0;
	double E_N_sig_events=0;
	double A_N_sig_events=0;
	double QE_N_sig_events=0;
	double E_N_bg_events=0;
	double A_N_bg_events=0;
	double QE_N_bg_events=0;
	double E_best_N_sig_events = 0;
	double A_best_N_sig_events = 0;
	double QE_best_N_sig_events = 0;
	double E_best_N_bg_events = 0;
	double A_best_N_bg_events = 0;
	double QE_best_N_bg_events = 0;

	double logchiU_step = 0.01;
	double logchiU_start = -7.0;
//int loopflag = 1;
//while(loopflag)
//{
	logchiU=logchiU_start; 
	while(E_sum < 1e3 && A_sum < 1e3 && QE_sum < 1e3)
	{

		in.mS = temp_mS;
		in.mZprime = temp_mZprime;

		chiU=pow(10.0,logchiU);
		E_sum=0.0;
		A_sum=0.0;
		QE_sum=0.0;

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram,qeGram);

		E_N_events = 0;
		E_N_sig_events = 0;
		E_N_bg_events = 0;
		A_N_sig_events = 0;
		A_N_bg_events = 0;
		A_N_events = 0;
		QE_N_sig_events = 0;
		QE_N_bg_events = 0;
		QE_N_events = 0;


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

			for(bin=0;bin<QEBINS;bin++)
			{
				temp_sig = sigma_s*qeGram[bin];
				temp_bg = (1.0+zeta_b)*qeB[bin];
				lambda = temp_sig + temp_bg;
				N = qeO[bin]; //MB has seen O[] events.
				
				QE_N_events += lambda;
				QE_N_sig_events += temp_sig;
				QE_N_bg_events += temp_bg;
			//	A_sum+= (lambda-N)*(lambda-N)/lambda;
				QE_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
			}
	//	}

		//a guesstimate of the systematic error on the background.
		double sigma_zeta = 0.05;
		// add prior on zeta.
		E_sum += pow((zeta_b/sigma_zeta),2.0);
		A_sum += pow((zeta_b/sigma_zeta),2.0);
		QE_sum += pow((zeta_b/sigma_zeta),2.0);
		std::cout<<chiU<<" "<<E_sum<<" "<<A_sum<<" "<<QE_sum<<std::endl;

		logchiU += logchiU_step;

	}

return 0.0;
}

BF_RESULT * fit_spectra(CL_input in, double cutEfficiency, double events[NUMEVENTS][NUM_EVENT_OBS]){
	double zeta_b = 0.0;
	double sigma_s = 1.0;
	

	double N=0.0;
	double lambda=0.0;
	double E_sum=0.0;
	double A_sum=0.0;
	double QE_sum=0.0;
	
	double eGram[EBINS];
	double cosGram[COSBINS];
	double qeGram[QEBINS];

	int i;
	for(i=0;i<EBINS;i++)
	{
		eGram[i]=0.0;
	}
	for(i=0;i<COSBINS;i++)
	{
		cosGram[i]=0.0;
	}
	for(i=0;i<QEBINS;i++)
	{
		qeGram[i]=0.0;
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
	double qeO[QEBINS];
	double qeB[QEBINS];
      	qeO[0] =  232; qeO[1] = 156 ;qeO[2] = 156 ;qeO[3] = 79 ;qeO[4] = 81 ;qeO[5] = 70 ;qeO[6] = 63 ;qeO[7] = 65 ;qeO[8] = 62 ;qeO[9] = 34 ;qeO[10] = 70;
        qeB[0] = 181.1 ;qeB[1] = 108.4;qeB[2] = 120.4;qeB[3] = 64.2;qeB[4] = 90.3; qeB[5] = 67.7;qeB[6] = 70.4;qeB[7] = 57.5;qeB[8] = 52.3; qeB[9] = 39;qeB[10] = 70.2;



	double temp_tot_E =0.0;
	double temp_tot_A =0.0;
	double temp_tot_QE =0.0;
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
	double best_QE_contEfficiency=1e50;
	double best_E_N_events=1e-50;
	double best_A_N_events=1e-50;
	double best_QE_N_events=1e-50;
	double E_best_chiU = 1e50;
	double E_best = 1e4;
	double A_best_chiU = 1e50;
	double A_best = 1e4;
	double QE_best_chiU = 1e50;
	double QE_best = 1e6;
	double temp_mS = in.mS;
	double temp_mZprime = in.mZprime;

	double E_N_events=0;
	double A_N_events=0;
	double QE_N_events=0;
	double E_N_sig_events=0;
	double A_N_sig_events=0;
	double QE_N_sig_events=0;
	double E_N_bg_events=0;
	double A_N_bg_events=0;
	double QE_N_bg_events=0;
	double E_best_N_sig_events = 0;
	double A_best_N_sig_events = 0;
	double QE_best_N_sig_events = 0;
	double E_best_N_bg_events = 0;
	double A_best_N_bg_events = 0;
	double QE_best_N_bg_events = 0;

	double logchiU_step = 0.1;
	double logchiU_start = -7.0;


// Whats the background only fit. for QE
//			double armalot = 0.0;
//			for(int bin=0;bin<QEBINS;bin++)
//			{
//				armalot += 2.0*(qeB[bin]-qeO[bin]) + 2.0*qeO[bin]*log(qeO[bin]/qeB[bin]);
//				std::cout<<"arm: "<<armalot<<std::endl;
//			} 
//			std::cout<<"# background fit is: "<<armalot<<std::endl;
	




//int loopflag = 1;
//while(loopflag)
//{
	logchiU=logchiU_start;
	while(E_sum < MCHI && A_sum < MCHI  && QE_sum < MCHI && logchiU < log(sqrt(boundU(temp_mS)*boundUtau(temp_mS))*sqrt(boundChi(temp_mZprime)))/log(10.0)) 
	//while(E_sum < MCHI && A_sum < MCHI  && QE_sum < MCHI  && logchiU < log(sqrt(boundU(temp_mS))*sqrt(boundChi(temp_mZprime)))/log(10.0)) 
	//while(E_sum < MCHI && A_sum < MCHI && QE_sum < MCHI) 
	{

		in.mS = temp_mS;
		in.mZprime = temp_mZprime;

		chiU=pow(10.0,logchiU);
		E_sum=0.0;
		A_sum=0.0;
		QE_sum=0.0;

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram,qeGram);

		E_N_events = 0;
		E_N_sig_events = 0;
		E_N_bg_events = 0;
		A_N_sig_events = 0;
		A_N_bg_events = 0;
		A_N_events = 0;
		QE_N_sig_events = 0;
		QE_N_bg_events = 0;
		QE_N_events = 0;

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

			for(bin=0;bin<QEBINS;bin++)
			{
				temp_sig = sigma_s*qeGram[bin];
				temp_bg = (1.0+zeta_b)*qeB[bin];
				lambda = temp_sig + temp_bg;
				N = qeO[bin]; //MB has seen O[] events.
				
				QE_N_events += lambda;
				QE_N_sig_events += temp_sig;
				QE_N_bg_events += temp_bg;
			//	QE_sum+= (lambda-N)*(lambda-N)/lambda;
				QE_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
			}
	//	}

		//a guesstimate of the systematic error on the background.
		double sigma_zeta = 0.05;
		// add prior on zeta.
		E_sum += pow((zeta_b/sigma_zeta),2.0);
		A_sum += pow((zeta_b/sigma_zeta),2.0);
		QE_sum += pow((zeta_b/sigma_zeta),2.0);
//		std::cout<<chiU<<" "<<E_sum<<" "<<A_sum<<std::endl;

		logchiU += logchiU_step;

	}

logchiU -= logchiU_step;
double breakpoint_logchiU = logchiU;
double E_sum_previous = E_sum + 1e-8;
double A_sum_previous = A_sum + 1e-8;
double QE_sum_previous = QE_sum + 1e-8;
double temp_E_sum = 0.0;
double temp_A_sum = 0.0;
double temp_QE_sum = 0.0;
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

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram,qeGram);

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

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram,qeGram);

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

	logchiU = breakpoint_logchiU; 
	
	while( QE_sum < QE_sum_previous)
	{

		temp_QE_sum = QE_sum;
	
		logchiU -= logchiU_step;

//printf("I did a thing. Again.\n");
		in.mS = temp_mS;
		in.mZprime = temp_mZprime;

		chiU=pow(10.0,logchiU);
		QE_sum=0.0;

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram,qeGram);

		QE_N_events = 0;
		QE_N_sig_events = 0;
		QE_N_bg_events = 0;

		int bin = 0;
		double temp_sig,temp_bg;
	//	if (sigma_s > 1e-5){
			for(bin=0;bin<QEBINS;bin++)
			{
				temp_sig = sigma_s*qeGram[bin];
				temp_bg = (1.0+zeta_b)*qeB[bin];
				lambda = temp_sig + temp_bg;
				N = qeO[bin]; //MB has seen O[] events.
				
				QE_N_events += lambda;
				QE_N_sig_events += temp_sig;
				QE_N_bg_events += temp_bg;
			//	QE_sum+= (lambda-N)*(lambda-N)/lambda;
				QE_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
			}
	//	}

		//a guesstimate of the systematic error on the background.
		double sigma_zeta = 0.05;
		// add prior on zeta.
		QE_sum += pow((zeta_b/sigma_zeta),2.0);
//		std::cout<<chiU<<" "<<E_sum<<" "<<A_sum<<std::endl;

		if(QE_sum < QE_best)
		{ 
			QE_best=QE_sum; 
			QE_best_chiU = chiU; 
			best_QE_contEfficiency = contEfficiency; 
			best_QE_N_events = QE_N_events;
			QE_best_N_sig_events = QE_N_sig_events;
			QE_best_N_bg_events = QE_N_bg_events;
		}
		
		QE_sum_previous = temp_QE_sum;

	} 


	if(fabs(E_best-65.1554) < 1e-3) { E_best=65.1554; }
	if(fabs(A_best-38.9753) < 1e-3) { A_best=38.9753; }
	if(fabs(QE_best-49.4823) < 1e-3) { QE_best=49.4823;}

//	double E_check = E_best_chiU*E_best_chiU*best_E_contEfficiency*getTotalNumEvents(in)*cutEfficiency;
//	double A_check = A_best_chiU*A_best_chiU*best_A_contEfficiency*getTotalNumEvents(in)*cutEfficiency;

//	std::cout<<E_check<<" "<<A_check<<std::endl;
	std::cout<<"# From Mark: "<<getTotalNumEvents(in)<<" Energy: "<<best_E_N_events<<"= ("<<E_best_N_sig_events<<" + "<<E_best_N_bg_events<<") Angle: "<<best_A_N_events<<"= ("<<A_best_N_sig_events<<" + "<<A_best_N_bg_events<<") QE: "<<best_QE_N_events<<"= ("<<QE_best_N_sig_events<<" + "<<QE_best_N_bg_events<<") "<<std::endl;

//	std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<65.1554-E_best<<" "<<E_best_chiU<<" "<<38.9753-A_best<<" "<<A_best_chiU<<std::endl;
	std::cout<<"#"<<"|ms|"<<" | "<<"mz"<<"| "<<"DelChiE"<<" | "<<"U at BF E"<<" | "<<"ContEff at BF E"<<" | "<<"DelChiA"<<" | "<<"U at BF A"<<" | "<<"ContEff at BF A"<<" | "<<"CutEff"<<" | "<<"DelChiQE"<<" | "<<"U at BF QE"<<" | "<<"ContEff at BF QE"<<std::endl;

	std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<65.1554-E_best<<" "<<E_best_chiU<<" "<<best_E_contEfficiency<<" "<<39.9753-A_best<<" "<<A_best_chiU<<" "<<best_A_contEfficiency<<" "<<cutEfficiency<<" "<<49.4823-QE_best<<" "<<QE_best_chiU<<" "<<best_QE_contEfficiency<<std::endl;

BF_RESULT * output = (BF_RESULT *)malloc(sizeof(BF_RESULT));
output->E_bf = E_best_chiU;
output->A_bf = A_best_chiU;
output->QE_bf = QE_best_chiU;

return output;
}

BF_RESULT * stats_fit_spectra(CL_input in, double cutEfficiency, double events[NUMEVENTS][NUM_EVENT_OBS])
{
	double Ezeta_b = 0.0;
	double Azeta_b = 0.0;
	double QEzeta_b = 0.0;
	double sigma_s = 1.0;
	double sigma_zeta = in.Sigma_Zeta;

	
	double N=0.0;
	double lambda=0.0;
	double E_sum=0.0;
	double A_sum=0.0;
	double QE_sum=0.0;
	
	double eGram[EBINS];
	double cosGram[COSBINS];
	double qeGram[QEBINS];

	int i;
	for(i=0;i<EBINS;i++)
	{
		eGram[i]=0.0;
	}
	for(i=0;i<COSBINS;i++)
	{
		cosGram[i]=0.0;
	}
	for(i=0;i<QEBINS;i++)
	{
		qeGram[i]=0.0;
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
	double qeO[QEBINS];
	double qeB[QEBINS];
      	qeO[0] =  232; qeO[1] = 156 ;qeO[2] = 156 ;qeO[3] = 79 ;qeO[4] = 81 ;qeO[5] = 70 ;qeO[6] = 63 ;qeO[7] = 65 ;qeO[8] = 62 ;qeO[9] = 34 ;qeO[10] = 70;
        qeB[0] = 181.1 ;qeB[1] = 108.4;qeB[2] = 120.4;qeB[3] = 64.2;qeB[4] = 90.3; qeB[5] = 67.7;qeB[6] = 70.4;qeB[7] = 57.5;qeB[8] = 52.3; qeB[9] = 39;qeB[10] = 70.2;


	double temp_tot_E =0.0;
	double temp_tot_A =0.0;
	double temp_tot_QE =0.0;
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
	double best_QE_contEfficiency=1e50;
	double best_E_N_events=1e-50;
	double best_A_N_events=1e-50;
	double best_QE_N_events=1e-50;
	double E_best_chiU = 1e50;
	double E_best = 1e4;
	double A_best_chiU = 1e50;
	double A_best = 1e4;
	double QE_best_chiU = 1e50;
	double QE_best = 1e6;
	double temp_mS = in.mS;
	double temp_mZprime = in.mZprime;

	double E_N_events=0;
	double A_N_events=0;
	double QE_N_events=0;
	double E_N_sig_events=0;
	double A_N_sig_events=0;
	double QE_N_sig_events=0;
	double E_N_bg_events=0;
	double A_N_bg_events=0;
	double QE_N_bg_events=0;
	double E_best_N_sig_events = 0;
	double A_best_N_sig_events = 0;
	double QE_best_N_sig_events = 0;
	double E_best_N_bg_events = 0;
	double A_best_N_bg_events = 0;
	double QE_best_N_bg_events = 0;

	double logchiU_step = 0.1;
	double logchiU_start = -7.0;
//int loopflag = 1;
//while(loopflag)
//{
//#########################################################################

        std::vector<double > eVGram(eGram, eGram + EBINS);
        std::vector<double > cosVGram(cosGram, cosGram + COSBINS);
        std::vector<double > qeVGram(qeGram, qeGram + QEBINS);
        
        std::vector<double > bf_zeta_b;  
        double bf_chi ;
        std::vector<double > abf_zeta_b; 
        double abf_chi ;
        std::vector<double > qebf_zeta_b; 
        double qebf_chi ;
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
     
        std::vector<double > qeZeros(QEBINS, 0.0);
        std::vector<double > QEBF_bkg_only_zeta_b = {0};
        double QEBF_bkg_only_chi = 10000;
        nuisMarginalize(& QEBF_bkg_only_zeta_b, &QEBF_bkg_only_chi, &qeZeros,2,sigma_zeta);
        std::cout<<"# QE sigma_zeta Bf: "<<QEBF_bkg_only_zeta_b[0]<<" BFChi: "<<QEBF_bkg_only_chi<<std::endl;


//#########################################################################

	logchiU=logchiU_start;
	while(E_sum < MCHI && A_sum < MCHI && logchiU < log(sqrt(boundU(temp_mS)*boundUtau(temp_mS))*sqrt(boundChi(temp_mZprime)))/log(10.0)) 
	//while(E_sum < MCHI && A_sum < MCHI && logchiU < log(sqrt(boundU(temp_mS))*sqrt(boundChi(temp_mZprime)))/log(10.0)) 
	//while(E_sum < MCHI && A_sum < MCHI ) 
	{

		in.mS = temp_mS;
		in.mZprime = temp_mZprime;

		chiU=pow(10.0,logchiU);
		E_sum=0.0;
		A_sum=0.0;
		QE_sum=0.0;

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram,qeGram);
                        eVGram.assign(eGram, eGram+EBINS);
                        nuisMarginalize(&bf_zeta_b, &bf_chi, &eVGram,0,sigma_zeta);
                        Ezeta_b=bf_zeta_b[0];
                
                        cosVGram.assign(cosGram, cosGram+COSBINS);
                        nuisMarginalize(&abf_zeta_b, &abf_chi, &cosVGram,1,sigma_zeta);
                        Azeta_b=abf_zeta_b[0];

                        qeVGram.assign(qeGram, qeGram+QEBINS);
                        nuisMarginalize(&qebf_zeta_b, &qebf_chi, &qeVGram,2,sigma_zeta);
                        QEzeta_b=qebf_zeta_b[0];

		E_N_events = 0;
		E_N_sig_events = 0;
		E_N_bg_events = 0;
		A_N_sig_events = 0;
		A_N_bg_events = 0;
		A_N_events = 0;
		QE_N_sig_events = 0;
		QE_N_bg_events = 0;
		QE_N_events = 0;

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

			for(bin=0;bin<QEBINS;bin++)
			{
				temp_sig = sigma_s*qeGram[bin];
				temp_bg = (1.0+QEzeta_b)*qeB[bin];
				lambda = temp_sig + temp_bg;
				N = qeO[bin]; //MB has seen O[] events.
				
				QE_N_events += lambda;
				QE_N_sig_events += temp_sig;
				QE_N_bg_events += temp_bg;
			//	QE_sum+= (lambda-N)*(lambda-N)/lambda;
				QE_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
			}
	//	}

		//a guesstimate of the systematic error on the background.
	//	double sigma_zeta = 0.05;
		// add prior on zeta.
		E_sum += pow((Ezeta_b/sigma_zeta),2.0);
		A_sum += pow((Azeta_b/sigma_zeta),2.0);
		QE_sum += pow((QEzeta_b/sigma_zeta),2.0);
//		std::cout<<chiU<<" "<<E_sum<<" "<<A_sum<<std::endl;

		logchiU += logchiU_step;

	}

logchiU -= logchiU_step;
double breakpoint_logchiU = logchiU;
double E_sum_previous = E_sum + 1e-8;
double A_sum_previous = A_sum + 1e-8;
double QE_sum_previous = QE_sum + 1e-8;
double temp_E_sum = 0.0;
double temp_A_sum = 0.0;
double temp_QE_sum = 0.0;
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

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram,qeGram);
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

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram,qeGram);
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

	logchiU = breakpoint_logchiU;
	while( QE_sum < QE_sum_previous)
	{

		temp_QE_sum = QE_sum;
	
		logchiU -= logchiU_step;

//printf("I did a thing. Again.\n");
		in.mS = temp_mS;
		in.mZprime = temp_mZprime;

		chiU=pow(10.0,logchiU);
		QE_sum=0.0;

		contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram,qeGram);
                        qeVGram.assign(qeGram, qeGram+QEBINS);
                        nuisMarginalize(&qebf_zeta_b, &qebf_chi, &qeVGram,2,sigma_zeta);
                        QEzeta_b=qebf_zeta_b[0];

		QE_N_events = 0;
		QE_N_sig_events = 0;
		QE_N_bg_events = 0;

		int bin = 0;
		double temp_sig,temp_bg;
	//	if (sigma_s > 1e-5){
			for(bin=0;bin<QEBINS;bin++)
			{
				temp_sig = sigma_s*qeGram[bin];
				temp_bg = (1.0+QEzeta_b)*qeB[bin];
				lambda = temp_sig + temp_bg;
				N = qeO[bin]; //MB has seen O[] events.
				
				QE_N_events += lambda;
				QE_N_sig_events += temp_sig;
				QE_N_bg_events += temp_bg;
			//	QE_sum+= (lambda-N)*(lambda-N)/lambda;
				QE_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
			}
	//	}

		//a guesstimate of the systematic error on the background.
		//double sigma_zeta = 0.05;
		// add prior on zeta.
		QE_sum += pow((QEzeta_b/sigma_zeta),2.0);
//		std::cout<<chiU<<" "<<E_sum<<" "<<A_sum<<std::endl;

		if(QE_sum<QE_best)
		{ 
			QE_best=QE_sum; 
			QE_best_chiU = chiU; 
			best_QE_contEfficiency = contEfficiency; 
			best_QE_N_events = QE_N_events;
			QE_best_N_sig_events = QE_N_sig_events;
			QE_best_N_bg_events = QE_N_bg_events;
		}
		
		QE_sum_previous = temp_QE_sum;

	}


//	if(fabs(E_best-65.1554) < 1e-3) { E_best=65.1554; }
//	if(fabs(A_best-38.9753) < 1e-3) { A_best=38.9753; }
        if(fabs(E_best-EBF_bkg_only_chi) < 1e-3) { E_best=EBF_bkg_only_chi; }
        if(fabs(A_best-ABF_bkg_only_chi) < 1e-3) { A_best=ABF_bkg_only_chi; }
        if(fabs(QE_best-QEBF_bkg_only_chi) < 1e-3) { QE_best=QEBF_bkg_only_chi; }


//	double E_check = E_best_chiU*E_best_chiU*best_E_contEfficiency*getTotalNumEvents(in)*cutEfficiency;
//	double A_check = A_best_chiU*A_best_chiU*best_A_contEfficiency*getTotalNumEvents(in)*cutEfficiency;

//	std::cout<<E_check<<" "<<A_check<<std::endl;
	std::cout<<"# From Mark: "<<getTotalNumEvents(in)<<" Energy: "<<best_E_N_events<<"= ("<<E_best_N_sig_events<<" + "<<E_best_N_bg_events<<") Angle: "<<best_A_N_events<<"= ("<<A_best_N_sig_events<<" + "<<A_best_N_bg_events<<") QE: "<<best_QE_N_events<<"= ("<<QE_best_N_sig_events<<" + "<<QE_best_N_bg_events<<") "<<std::endl;

//	std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<65.1554-E_best<<" "<<E_best_chiU<<" "<<38.9753-A_best<<" "<<A_best_chiU<<std::endl;
//std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<65.1554-E_best<<" "<<E_best_chiU<<" "<<best_E_contEfficiency<<" "<<39.9753-A_best<<" "<<A_best_chiU<<" "<<best_A_contEfficiency<<" "<<cutEfficiency<<std::endl;
        std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<EBF_bkg_only_chi-E_best<<" "<<E_best_chiU<<" "<<best_E_contEfficiency<<" "<<ABF_bkg_only_chi-A_best<<" "<<A_best_chiU<<" "<<best_A_contEfficiency<<" "<<cutEfficiency<<"  "<<QEBF_bkg_only_chi-QE_best<<" "<<QE_best_chiU<<" "<<best_QE_contEfficiency<<std::endl;

BF_RESULT * output = (BF_RESULT *)malloc(sizeof(BF_RESULT));
output->E_bf = E_best_chiU;
output->A_bf = A_best_chiU;
output->QE_bf = QE_best_chiU;

return output;
}

BF_RESULT * stats_fit_spectra_indiv(CL_input in, double cutEfficiency, double events[NUMEVENTS][NUM_EVENT_OBS])
{	
	int which_var = in.which_var;

	double zeta_b =0.0;
	double sigma_s = 1.0;
	double sigma_zeta = in.Sigma_Zeta;
	
	double N=0.0;
	double lambda=0.0;
	double sum =0.0;
	
	int BINS = which_BINS(which_var);
	double spec_obs[BINS];
	double spec_bkg[BINS];

	if (which_var == 0)
	{	
		spec_obs[0] = 204.0;spec_obs[1] = 280.0;spec_obs[2] = 214.0;spec_obs[3] = 99.0;spec_obs[4] = 83.0;spec_obs[5] = 59.0;spec_obs[6] = 51.0;spec_obs[7] = 33.0;spec_obs[8] = 37.0;spec_obs[9] = 23.0; spec_obs[10] = 19.0;spec_obs[11] = 21.0;spec_obs[12] = 12.0; spec_obs[13] = 16.0; spec_obs[14] = 4.0;spec_obs[15] = 9.0; spec_obs[16] = 4.0; spec_obs[17] = 7.0; spec_obs[18] = 3.0;
		spec_bkg[0] = 151.5;	spec_bkg[1] = 218.8;	spec_bkg[2] = 155.6;spec_bkg[3] = 108.7;spec_bkg[4] = 72.5;spec_bkg[5] = 57.6;spec_bkg[6] = 45;spec_bkg[7] = 38.5;spec_bkg[8] = 31.4;spec_bkg[9] = 22.2;spec_bkg[10] = 20.4;spec_bkg[11] = 17.2;		spec_bkg[12] = 14.1;	spec_bkg[13] = 10.2;	spec_bkg[14] = 9.1;	spec_bkg[15] = 8.2;	spec_bkg[16] = 5.6;	spec_bkg[17] = 5.7;	spec_bkg[18] = 2.9;


	} else if (which_var ==1)
	{
	spec_obs[0] = 22; 	spec_obs[1] = 34; 	spec_obs[2] = 43; 	spec_obs[3] = 41; 	spec_obs[4] = 60; 	spec_obs[5] = 87; 	spec_obs[6] = 90; 	spec_obs[7] = 139; 	spec_obs[8] = 237; 	spec_obs[9] = 429;  	
spec_bkg[0] = 19.9; 	spec_bkg[1] = 23.1; 	spec_bkg[2] = 28.8; 	spec_bkg[3] = 32.1; 	spec_bkg[4] = 46.4; 	spec_bkg[5] = 63.1; 	spec_bkg[6] = 86.1; 	spec_bkg[7] = 121; 	spec_bkg[8] = 196.8; 	spec_bkg[9] = 390;


	} else if (which_var ==2)
	{
      	spec_obs[0] =  232; spec_obs[1] = 156 ;spec_obs[2] = 156 ;spec_obs[3] = 79 ;spec_obs[4] = 81 ;spec_obs[5] = 70 ;spec_obs[6] = 63 ;spec_obs[7] = 65 ;spec_obs[8] = 62 ;spec_obs[9] = 34 ;spec_obs[10] = 70;    
     spec_bkg[0] = 181.1 ;spec_bkg[1] = 108.4;spec_bkg[2] = 120.4;spec_bkg[3] = 64.2;spec_bkg[4] = 90.3; spec_bkg[5] = 67.7;spec_bkg[6] = 70.4;spec_bkg[7] = 57.5;spec_bkg[8] = 52.3; spec_bkg[9] = 39;spec_bkg[10] = 70.2;


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
	double best = 1e4;
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
        nuisMarginalize(& BF_bkg_only_zeta_b,& BF_bkg_only_chi, &Zeros, which_var, sigma_zeta);
        std::cout<<"# which: "<<which_var<<" ,sigma_zeta Bf: "<< BF_bkg_only_zeta_b[0]<<" BF Chi: "<<BF_bkg_only_chi<<std::endl;

//#########################################################################

	logchiU=logchiU_start;
	while(sum <  MCHI && logchiU < log(whatsmaxUXorder1(temp_mS,temp_mZprime))/log(10.0))
	//while(sum < MCHI && logchiU < log(sqrt(boundU(temp_mS)*boundUtau(temp_mS))*sqrt(boundChi(temp_mZprime)))/log(10.0)) 
	//while(sum < MCHI && logchiU < log(sqrt(boundUpeaky(temp_mS,temp_mZprime)*boundUtau(temp_mS))*sqrt(boundChi(temp_mZprime)))/log(10.0))
	//while(sum < MCHI && logchiU < log(boundUpeaky(temp_mS,temp_mZprime)*sqrt(boundChi(temp_mZprime)))/log(10.0))
	//while(E_sum < MCHI && A_sum < MCHI && logchiU < log(sqrt(boundU(temp_mS))*sqrt(boundChi(temp_mZprime)))/log(10.0)) 
	//while(E_sum < MCHI && A_sum < MCHI ) 
	{

		//if( logchiU > log(sqrt(boundU(temp_mS)*boundUtau(temp_mS))*sqrt(boundChi(temp_mZprime)))/log(10.0) )
		//{
		//	logchiU = log(sqrt(boundU(temp_mS)*boundUtau(temp_mS))*sqrt(boundChi(temp_mZprime)))/log(10.0);
		//}

		in.mS = temp_mS;
		in.mZprime = temp_mZprime;

		chiU=pow(10.0,logchiU);
		sum=0.0;

		contEfficiency = histogrammer_indiv(in,chiU,cutEfficiency,events,Gram,which_var);
                        VGram.assign(Gram, Gram+BINS);
                        nuisMarginalize(&bf_zeta_b, &bf_chi, &VGram,which_var,sigma_zeta);
                        zeta_b=bf_zeta_b[0];

		N_events = 0;
		N_sig_events = 0;
		N_bg_events = 0;

		int bin = 0;
		double temp_sig,temp_bg;
	//	if (sigma_s > 1e-5){
			for(bin=0;bin<BINS;bin++)
			{
				temp_sig = sigma_s*Gram[bin];
				temp_bg = (1.0+zeta_b)*spec_bkg[bin];
				lambda = temp_sig + temp_bg;
				N = spec_obs[bin]; //MB has seen O[] events.
				
				N_events += lambda;
				N_sig_events += temp_sig;
				N_bg_events += temp_bg;
			//	E_sum+= (lambda-N)*(lambda-N)/lambda;
				sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
			}

		// add prior on zeta.
		sum += pow((zeta_b/sigma_zeta),2.0);
//		std::cout<<chiU<<" "<<E_sum<<" "<<A_sum<<std::endl;

		logchiU += logchiU_step;

	}

logchiU -= logchiU_step;
double breakpoint_logchiU = logchiU;
double sum_previous = sum + 1e-8;
double temp_sum = 0.0;
logchiU_step = 0.0025;

//printf("I did a thing.\n");
	while( sum < sum_previous)
	{

		temp_sum = sum;
	
		logchiU -= logchiU_step;

//printf("I did a thing. Again.\n");
		in.mS = temp_mS;
		in.mZprime = temp_mZprime;

		chiU=pow(10.0,logchiU);
		sum=0.0;

		contEfficiency = histogrammer_indiv(in,chiU,cutEfficiency,events,Gram,which_var);
                        VGram.assign(Gram, Gram+BINS);
                        nuisMarginalize(&bf_zeta_b, &bf_chi, &VGram,which_var,sigma_zeta);
                        zeta_b=bf_zeta_b[0];

		N_events = 0;
		N_sig_events = 0;
		N_bg_events = 0;

		int bin = 0;
		double temp_sig,temp_bg;
	//	if (sigma_s > 1e-5){
			for(bin=0;bin<BINS;bin++)
			{
				temp_sig = sigma_s*Gram[bin];
				temp_bg = (1.0+zeta_b)*spec_bkg[bin];
				lambda = temp_sig + temp_bg;
				N = spec_obs[bin]; //MB has seen O[] events.
				
				N_events += lambda;
				N_sig_events += temp_sig;
				N_bg_events += temp_bg;
			//	E_sum+= (lambda-N)*(lambda-N)/lambda;
				sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
			}
	//	}

		//a guesstimate of the systematic error on the background.
		//double sigma_zeta = 0.05;
		// add prior on zeta.
		sum += pow((zeta_b/sigma_zeta),2.0);
//		std::cout<<chiU<<" "<<E_sum<<" "<<A_sum<<std::endl;

		if(sum<best)
		{ 
			best=sum; 
			best_chiU = chiU; 
			best_contEfficiency = contEfficiency; 
			best_N_events = N_events;
			best_N_sig_events = N_sig_events;
			best_N_bg_events = N_bg_events;
		}
		
		sum_previous = temp_sum;

	}


//	if(fabs(E_best-65.1554) < 1e-3) { E_best=65.1554; }
//	if(fabs(A_best-38.9753) < 1e-3) { A_best=38.9753; }
        if(fabs(best-BF_bkg_only_chi) < 1e-3) { best=BF_bkg_only_chi; }


//	std::cout<<E_check<<" "<<A_check<<std::endl;
	std::cout<<"# From Mark: "<<getTotalNumEvents(in)<<" which variable:  "<<which_var<<" : "<<best_N_events<<" = ("<<best_N_sig_events<<" + "<<best_N_bg_events<<")"<<std::endl;
		//	std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<65.1554-E_best<<" "<<E_best_chiU<<" "<<38.9753-A_best<<" "<<A_best_chiU<<std::endl;
//std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<65.1554-E_best<<" "<<E_best_chiU<<" "<<best_E_contEfficiency<<" "<<39.9753-A_best<<" "<<A_best_chiU<<" "<<best_A_contEfficiency<<" "<<cutEfficiency<<std::endl;
        std::cout<<temp_mS<<" "<<temp_mZprime<<" "<<BF_bkg_only_chi-best<<" "<<best_chiU<<" "<<best_contEfficiency<<" "<<cutEfficiency<<std::endl;

BF_RESULT * output = (BF_RESULT *)malloc(sizeof(BF_RESULT));
output->E_bf = best_chiU;

return output;
}



int main(int argc, char * argv[])
{

	//std::cout<<whatsmaxUX(0.08,0.2)<<std::endl;		

	static CL_input in;
	in.eCut = 0.0;
	in.thCut = 180.0;
	in.eFloor = 0.0; // I believe 0.0 for Floor and Ratio mean that the cut is never passed.
	in.eRatio = 0.0;
        in.Sigma_Zeta = 0.20;
	in.which_var = 0;
	double chiU = 1e-10;

	int c;
	int printFlag = 0;
	int modeFlag = 0;
	int statsFlag = 0;

	opterr = 0;

	while ((c = getopt (argc, argv, "m:Z:X:c:S:B:I:TPFEROAQ")) != -1)
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
    		case 'B':
			if(!strcmp(optarg,"E")){ modeFlag = 5; }
			else if(!strcmp(optarg,"A")){ modeFlag = 6; }
			else if(!strcmp(optarg,"Q")){ modeFlag = 9; }
			else { printf("Very bad thing!\nAborting...\n\nYou're probably not using the -B flag correctly.\n\n"); exit(1); }
			break;
		case 'I':
			modeFlag=10;
			if(!strcmp(optarg,"E")){ in.which_var = 0; }
			else if(!strcmp(optarg,"A")){ in.which_var = 1; }
			else if(!strcmp(optarg,"Q")){ in.which_var= 2; }
			else { printf("Very bad thing!\nAborting...\n\nYou're probably not using the -B flag correctly.\n\n"); exit(1); }
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
      		case '?':
			printf("Abandon hope all ye who enter this value: %c\n",optopt);
			printf("Allowed arguments:\n\t-m\tsets mS mass.\n\t-Z\tsets mZprime mass.\n\t-X\tsets chiU\n\t-c\tsets cuts (0.14,5.0).\n\t-P\tprints input parameters.\n\t-E(-A)(-Q)\tprints the (E)nergy, (A)ngular spectrum or (Q)e spectrum.\n\t-I \t Individual run of either the (IB)nergy, (IA)ngular or (IQ)e minimizer (NEED STATS ON).\n\t-B \tBrints the (BE)nergy, (BA)ngular spectrum or (BQ)e spectrum.\n\t-F\tprints mimima\n\t-R\tturns on energy asymmetry cut\n\t-O\tprints the cut efficiency.\n\t-S\t Toggles Systematics of bkg (Default off)\n\t-T\t Testing ground, god knows what you will find.\n");
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

		//std::cout<<"# Start Mode Flag"<<std::endl;
		int i;
		double eGram[EBINS];
		double cosGram[COSBINS];
		double qeGram[QEBINS];
		double contEfficiency;

		//A slightly hacky way to print the best-fit spectra.
		if(modeFlag == 5 || modeFlag == 6 || modeFlag == 9)// 5: print best fitting energy spec. 6: print best fitting angular spec. 8: print best fitting QE spectrum
		{
			BF_RESULT * bestfit = (BF_RESULT *)malloc(sizeof(BF_RESULT));
			bestfit = stats_fit_spectra(in,cutEfficiency,events);

			if (modeFlag == 5){ chiU = bestfit->E_bf; modeFlag = 1; statsFlag=0;}
			else if (modeFlag ==6){ chiU = bestfit->A_bf; modeFlag = 2; statsFlag=0;}
			else {chiU = bestfit->QE_bf; modeFlag = 8; statsFlag=0; }


			contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram,qeGram);
		}
		else
		{
			contEfficiency = histogrammer(in,chiU,cutEfficiency,events,eGram,cosGram,qeGram);
		}


		// The main modeFlag fractionating column.
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
		else if(modeFlag == 8)
		{ 
			for(i=0;i<QEBINS;i++)
			{
				std::cout<<BintoCentalQE(i)<<" "<<qeGram[i]<<std::endl;
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
			
			stats_fit_spectra_indiv(in,cutEfficiency,events);


		} else { std::cout<<"BAD THING!"<<std::endl; }
	}
	else 
	{
		if(statsFlag)
		{
			stats_fit_spectra(in,cutEfficiency,events);
		}
		else
		{
			//std::cout<<"# Begin End of chain,  fit_spectra"<<std::endl;
			fit_spectra(in,cutEfficiency,events);
		}	
	}

return 0;
}
