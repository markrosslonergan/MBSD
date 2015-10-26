#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <nlopt.hpp>
#include <iomanip>
#include <sys/time.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_gamma.h>
//#include <gsl/gsl_sf_gamma_inc.h>
//#include <gsl/pow.h>

#include "LR.h"
//#include "Minuit2/FCNBase.h"

#define POSNU 0
#define POSNUBAR 1
#define NEGNU 2
#define NEGNUBAR 3


double GoF1(std::vector<double > Spectrum, std::vector<double > Obs, std::vector<double > bkg, double bkg_zeta){
	double ans = 0;
	double N=0;
	int length = Obs.size();
	for(int i = 0; i<length; i++)
	{	
		N=Spectrum[i]+(bkg_zeta+1.0)*bkg[i];
		ans += 2*(N-Obs[i]+Obs[i]*log(Obs[i]/N));
	}
	return ans;
}
double GoF2(std::vector<double > Spectrum, std::vector<double > Obs, std::vector<double > bkg, double bkg_zeta){
	double ans = 0;
	double N=0;
	int length = Obs.size();
	for(int i = 0; i<length; i++)
	{	
		N=Spectrum[i]+(1.0+bkg_zeta)*bkg[i];
		ans+=pow(N-Obs[i],2)/Obs[i];
	}
	return ans;
}

double pval(double chi2, double ndof){

	return 1.0-(gsl_sf_gamma_inc(ndof/2.0, 0.0)-gsl_sf_gamma_inc(ndof/2.0,chi2/2.0))/gsl_sf_gamma(ndof/2);

}



int mcPrintVar(std::vector<double > list){
for(int i = 0; i< list.size(); i++) 
{
	std::cout<<list[i]<<"  ";
}
return 1;
}


double intpow( double base, int exponent )
{
int i;
double out = base;
for( i=1 ; i < exponent ; i++ )
{
out *= base;
}
return out;
}


double getTotalNumEvents(CL_input input)
{
	double mS = input.mS;
	double mZprime = input.mZprime;

	FILE *ptr_file;
    	
	char buf[3000];

	char * pch;
	
	double totNum = -1e5;

	int n = 1;
	int m = 0;
	char s[100];
//	char filename[500] = "/scratch/ross/git/MBSD/MC/HIST_\0";
//	char filename[500] = "/scratch/ross/dataMC24sept/HIST_\0";
	char filename[500] = "/scratch/ross/dataMC12oct/anti/HIST_\0";
	sprintf(s,"%.3lf_%.3lf.dat", mS, mZprime);
	strcat(filename,s);
//	printf("Total Filename: %s\n",filename);
	ptr_file =fopen(filename,"r");

    	if (!ptr_file)
       	{			
		printf("ERROR LOADING MC EVENTS 2\n");
		return 1;
	}
    	while (fgets(buf,3000, ptr_file)!=NULL)
	{
		if(m==1)
		{
			pch = strtok(buf,"\t");
			n=1;
 			while (pch != NULL)
			{
//				printf("%.7g, ", strtof(pch,NULL));
				if(n==3){	//printf("%.7g, ",strtof(pch,NULL));
						totNum = strtof(pch,NULL);	 
					 	fclose(ptr_file);
						return totNum;
					}
				else { 
					pch = strtok(NULL,"\t"); 
				}
				n++;	
			}
		}
		m++;
	}

fclose(ptr_file);
return totNum;
}
double getTotalNumEvents_NU(CL_input input)
{
	double mS = input.mS;
	double mZprime = input.mZprime;

	FILE *ptr_file;
    	
	char buf[3000];

	char * pch;
	
	double totNum = -1e5;

	int n = 1;
	int m = 0;
	char s[100];

	char filename[500]= "/scratch/ross/dataMC12oct/HIST_\0";

	sprintf(s,"%.3lf_%.3lf.dat", mS, mZprime);
	strcat(filename,s);
//	printf("Total Filename: %s\n",filename);
	ptr_file =fopen(filename,"r");

    	if (!ptr_file)
       	{			
		printf("ERROR LOADING MC EVENTS 2\n");
		return 1;
	}
    	while (fgets(buf,3000, ptr_file)!=NULL)
	{
		if(m==1)
		{
			pch = strtok(buf,"\t");
			n=1;
 			while (pch != NULL)
			{
//				printf("%.7g, ", strtof(pch,NULL));
				if(n==3){	//printf("%.7g, ",strtof(pch,NULL));
						totNum = strtof(pch,NULL);	 
					 	fclose(ptr_file);
						return totNum;
					}
				else { 
					pch = strtok(NULL,"\t"); 
				}
				n++;	
			}
		}
		m++;
	}

fclose(ptr_file);
return totNum;
}
double getTotalNumEvents_NUBAR(CL_input input)
{
	double mS = input.mS;
	double mZprime = input.mZprime;
	double FUDGE = 3;
	FILE *ptr_file;
    	
	char buf[3000];

	char * pch;
	
	double totNum = -1e5;

	int n = 1;
	int m = 0;
	char s[100];
	char  filename[500]="/scratch/ross/dataMC12oct/anti/HIST_\0";
	

	sprintf(s,"%.3lf_%.3lf.dat", mS, mZprime);
	strcat(filename,s);
//	printf("Total Filename: %s\n",filename);
	ptr_file =fopen(filename,"r");

    	if (!ptr_file)
       	{			
		printf("ERROR LOADING MC EVENTS 2\n");
		return 1;
	}
    	while (fgets(buf,3000, ptr_file)!=NULL)
	{
		if(m==1)
		{
			pch = strtok(buf,"\t");
			n=1;
 			while (pch != NULL)
			{
//				printf("%.7g, ", strtof(pch,NULL));
				if(n==3){	//printf("%.7g, ",strtof(pch,NULL));
						totNum = strtof(pch,NULL);	 
					 	fclose(ptr_file);
						return totNum;
					}
				else { 
					pch = strtok(NULL,"\t"); 
				}
				n++;	
			}
		}
		m++;
	}

fclose(ptr_file);
return FUDGE*totNum;
}
int EtoBin(double E)
{
	return floor(E/0.1);
}
double BintoCentralE(int bin)
{
	return 0.1*(((double)bin)+0.5);
}

int CostoBin(double C)
{
	return floor((1+C)/0.2);
}
double BintoCentralCos(int b)
{
	return (((double)b)+0.5)*0.2-1.0;
}

int QEtoBin(double QE)
{
	int ans=0;
	std::vector<double > QEbinLow = {0.2,0.3,0.375,0.475,0.55,0.675,0.8,0.95,1.1,1.3,1.5};
	std::vector<double > QEbinDel = {0.1,0.075,0.1,0.075,0.125,0.125,0.15,0.15,0.2,0.2,1.5};
	for(int i = 0;i<QEBINS;i++){
		if(QEbinLow[i]<= QE &&  QE <= (QEbinLow[i]+QEbinDel[i]))
		{
			ans = i;
		}
	
	}
	if(QE < QEbinLow[0] || QE > (QEbinLow[QEBINS-1]+QEbinDel[QEBINS-1]) )
		{
		ans=999;
		}
	return ans;
}

double BintoCentalQE(int b)
{
	std::vector<double > QEbinLow = {0.2,0.3,0.375,0.475,0.55,0.675,0.8,0.95,1.1,1.3,1.5};
	std::vector<double > QEbinDel = {0.1,0.075,0.1,0.075,0.125,0.125,0.15,0.15,0.2,0.2,1.5};
	return QEbinLow[b]+QEbinDel[b]/2.0;
}

int which_BINS(int which_var){
	int BINS =0;

	if (which_var == 0)
	{	
		BINS=EBINS;
	} else if (which_var ==1)
	{ 
		BINS=COSBINS;
	} else if (which_var ==2)
	{
		BINS=QEBINS;
	}

return BINS;
}

double QEfromEandCos(double Evis, double costh){
	double Mn = 0.9395666;
	double Mp = 0.938272;
	double Dms= Mn*Mn-Mp*Mp;
	double EB = 0.00786; //Binging energy of carbon per nuclear
	double Me = 0.000510998;

	double c1= -0.1472; double c2=0.2788; double c3=0.0898;double c4=-0.0196;
	double a1= 0.9942; double a2=0.0113; //MC corrections

	double EvisCor= a1*(Evis-Me)+a2+Me;
	double beta = sqrt(1-Me*Me/(EvisCor*EvisCor));
	double ans1 =  0.5*((2.0*Mn+EB)*EvisCor-(Dms+2*Mn*EB+EB*EB+Me*Me))/((Mn+EB)-EvisCor+sqrt(EvisCor*EvisCor-Me*Me)*costh);
	double Q2=2*ans1*EvisCor*(1-beta*costh)-Me*Me;
	
	//return  0.5*((2.0*Mn+EB)*Evis-(Dms+2*Mn*EB+EB*EB+Me*Me))/((Mn+EB)-Evis+sqrt(EvisCor*Evis-Me*Me)*costh);
	return ans1-(c1+c2*Q2+c3*Q2*Q2+c4*intpow(Q2,3));
}


double decayProb(CL_input input, double chiU, double Es)
{
	double mS = input.mS;
	double mZprime = input.mZprime;

	double L = 1.5; // 1 metre.
	double g1 = 0.2;
	double e = sqrt(4*M_PI*1.0/137);
	double tw = 0.523;

	

	double alpha = mS*mS/(mZprime*mZprime-mS*mS);

	double prefac =  (1/(4*M_PI*M_PI))*2*M_PI*M_PI*mS*pow(  (chiU*cos(tw)*e)/(4*4*M_PI*M_PI), 2)*alpha*alpha/2; 
//        double  prefac   (chiU*chiU*(1.0/(137.05*137.05))*alpha*alpha*mS/(4.0*intpow(4.0*M_PI,3.0)*g1*g1));
	//double func = (3.0/(2.0*alpha*alpha*alpha))*((2+alpha-3*alpha*alpha)/(1+alpha))+(4*alpha*alpha-3)*log(1.0+alpha)/intpow(alpha,4.0);
	double func;
	if(alpha < 0.001)
	{
	
		double func_perturb_0 = (1.0/(1.0+alpha))*(7.0/4.0 + 41.0/60.0*alpha); 

		double func_perturb_rest = -0.18333*intpow(alpha,2)+0.22857*intpow(alpha,3)-0.23274*intpow(alpha,4)+0.22421*intpow(alpha,5)-0.21190*intpow(alpha,6)+0.19899*intpow(alpha,7)-0.18662*intpow(alpha,8)+0.17517*intpow(alpha,9)-0.16475*intpow(alpha,10)+0.15531*intpow(alpha,11);
		func = func_perturb_0+func_perturb_rest;
	}
	else 
	{		
		func = (3.0/(2.0*intpow(alpha,3)))*(2+alpha-3*alpha*alpha)/(1.0+alpha) + (4.0*alpha*alpha-3.0)*(log(1+alpha)/log(exp(1.0)))/intpow(alpha,4);
	}


//printf("### mS: %.5lf\tmZ: %.5lf\talpha: %.5lf\tprefac: %.5lf\tfunc: %.5lf\n",mS,mZprime,alpha,prefac,func); 

return 1.0-exp(-(0.51e16)*prefac*func*L*mS/Es);
} 


double histogrammer(CL_input in, double chiU, double cutEff, const double events[][NUM_EVENT_OBS], double eGram[], double cosGram[], double qeGram[])
{
	double finalScale = getTotalNumEvents(in);


	int i=0;
	for(i=0;i<COSBINS;i++)
	{
		cosGram[i]=0.0;
	}

	for(i=0;i<EBINS;i++)
	{
		eGram[i]=0.0;
	}
	
	for(i=0;i<QEBINS;i++)
	{
		qeGram[i]=0.0;
	}

	double total = 0.0;
	double probTotal = 0.0;
	double prob = 0.0;
	int qeWhichBin = 0;
	double QEprob=0.0;
	double QEtotal = 0.0;
	double QEprobTotal = 0.0;

	for(i=0;i<NUMEVENTS;i++)
	{
		if(events[i][0] > 1e-4)
		{			
			//These lines DO NOT CONSIDER containment effects.
//			eGram[EtoBin(events[i][0])] +=1.0;
//			cosGram[CostoBin(events[i][2])] +=1.0;
			
			//These lines account for containment effects.
			prob = decayProb(in,chiU,events[i][3]);
			QEprob=prob;
			eGram[EtoBin(events[i][0])] += prob;
			cosGram[CostoBin(cos((M_PI/180.0)*events[i][1]))] += prob;

			qeWhichBin = QEtoBin(QEfromEandCos(events[i][0],cos((M_PI/180.0)*events[i][1])));

			if(qeWhichBin == 999) //Error Bin 
			{
				QEprob=0;
				QEtotal+=0.0;
				QEprobTotal+=QEprob;
			} else {
				QEprob=prob;
				qeGram[qeWhichBin] += QEprob;	
				QEtotal+=1.0;
				QEprobTotal+=QEprob;
			}
			
			total+=1.0;
			probTotal+=prob;
			
//		        printf("event %d: prob = %.5g\tprobTot = %.5g\n",i,prob,probTotal);

		}
	}

	double cosTot = 0.0;

	for(i=0;i<COSBINS;i++)
	{
		cosGram[i]=chiU*chiU*cutEff*finalScale*cosGram[i]/total;
		cosTot+=cosGram[i];
	}

	for(i=0;i<EBINS;i++)
	{
		eGram[i]=chiU*chiU*cutEff*finalScale*eGram[i]/total;
	}

	for(i=0;i<QEBINS;i++)
	{
		qeGram[i]=chiU*chiU*cutEff*finalScale*qeGram[i]/total;
	}

//printf("chiU: %.5g\tTot: %.5g\tprobTot/total: %.5g\tcosTot: %.5g\n",chiU, total, probTotal/total, cosTot);

return probTotal/total;
}

double histogrammer_indiv(CL_input in, double chiU,  double cutEff, const double events[][NUM_EVENT_OBS], double Gram[], int which_var)
{
	double finalScale = getTotalNumEvents(in);

	int BINS=which_BINS(which_var);
	
	int i=0;
	for(i=0;i<BINS;i++)
	{
		Gram[i]=0.0;
	}
	
	double total = 0.0;
	double probTotal = 0.0;
	double prob = 0.0;
	int qeWhichBin = 0;
	

	for(i=0;i<1000;i++)
	{
		if(events[i][0] > 1e-4)
		{			
			//These lines DO NOT CONSIDER containment effects.
//			eGram[EtoBin(events[i][0])] +=1.0;
//			cosGram[CostoBin(events[i][2])] +=1.0;
			
			//These lines account for containment effects.
			
			
			prob = decayProb(in,chiU,events[i][3]);
			if( which_var ==0 ){

				
				
				      	Gram[EtoBin(events[i][0])] += prob;
					total+=1.0;
					probTotal+=prob;
	
			} else if (which_var==1)
		       	{
					Gram[CostoBin(cos((M_PI/180.0)*events[i][1]))] += prob;
					total+=1.0;
					probTotal+=prob;
			}  else if (which_var==2) 
			{
				qeWhichBin = QEtoBin(QEfromEandCos(events[i][0],cos((M_PI/180.0)*events[i][1])));
	
				if(qeWhichBin == 999) //Error Bin 
				{
					prob=0;
					total+=0.0;
				} else {
					Gram[qeWhichBin] += prob;	
					total+=1.0;
					probTotal+=prob;
				}
			}
			
		
//		        printf("event %d: prob = %.5g\tprobTot = %.5g\n",i,prob,probTotal);

		}
	}

	

	for(i=0;i<BINS;i++)
	{
		Gram[i]=chiU*chiU*cutEff*finalScale*Gram[i]/total;
	}
//printf("chiU: %.5g\tTot: %.5g\tprobTot/total: %.5g\tcosTot: %.5g\n",chiU, total, probTotal/total, cosTot);

return probTotal/total;
}


double histogrammer_indiv2(CL_input in, double Up, double Ud, double chi, double cutEff, const double events[][NUM_EVENT_OBS], double Gram[], int which_var, double finalScale)
{
	int BINS=which_BINS(which_var);
	
	int i=0;
	for(i=0;i<BINS;i++)
	{
		Gram[i]=0.0;
	}
	
	double total = 0.0;
	double probTotal = 0.0;
	double prob = 0.0;
	int qeWhichBin = 0;

	for(i=0;i<NUMEVENTS;i++)//NUMEVENTS
	{
		if(events[i][0] > 1e-4)
		{			
			//These lines DO NOT CONSIDER containment effects.
//			eGram[EtoBin(events[i][0])] +=1.0;
//			cosGram[CostoBin(events[i][2])] +=1.0;
			
			//These lines account for containment effects.
			
			
			prob = decayProb(in,Ud*chi,events[i][3]);
		
			switch(which_var) {
				case 0:
					Gram[EtoBin(events[i][0])] += prob;
					total+=1.0;
					probTotal+=prob;
					break;
				case 1:
					Gram[CostoBin(cos((M_PI/180.0)*events[i][1]))] += prob;
					total+=1.0;
					probTotal+=prob;
					break;
				case 2:
					qeWhichBin = QEtoBin(QEfromEandCos(events[i][0],cos((M_PI/180.0)*events[i][1])));
	
					if(qeWhichBin == 999) //Error Bin 
					{
						prob=0;
						total+=0.0;
					} else {
						Gram[qeWhichBin] += prob;	
						total+=1.0;
						probTotal+=prob;
					}
					break;

				default:							
					std::cout<<"ERROR: which_var in histogrammer must be 0(E), 1(A) or 2(QE)"<<std::endl;
			}



//		        printf("event %d: prob = %.5g\tprobTot = %.5g\n",i,prob,probTotal);

		}
	}

	

	for(i=0;i<BINS;i++)
	{
		Gram[i]=Up*Up*chi*chi*cutEff*finalScale*Gram[i]/total;
	}
//printf("chiU: %.5g\tTot: %.5g\tprobTot/total: %.5g\tcosTot: %.5g\n",chiU, total, probTotal/total, cosTot);

return probTotal/total;
}

double histogrammer_indiv2_dual(CL_input in, double Up, double Ud, double chi, double cutEff_NU,double cutEff_NUBAR, const double events_NU[][NUM_EVENT_OBS], const double events_NUBAR[][NUM_EVENT_OBS], double Gram[], int which_var, double finalScale_NU, double finalScale_NUBAR)
{
	int BINS=2*which_BINS(which_var);
	
	int i=0;
	for(i=0;i<BINS;i++)
	{
		Gram[i]=0.0;
	}

	int NUstarter = 0;
	int NUBARstarter = BINS/2;


	double total_NU = 0.0;
	double probTotal_NU = 0.0;
	double prob_NU = 0.0;
	double total_NUBAR = 0.0;
	double probTotal_NUBAR = 0.0;
	double prob_NUBAR = 0.0;
	int qeWhichBin = 0;

	for(int i=0;i<NUMEVENTS;i++)//NUMEVENTS first NU run
	{

		if(events_NU[i][0] > 1e-4)
		{			
			//These lines DO NOT CONSIDER containment effects.
//			eGram[EtoBin(events[i][0])] +=1.0;
//			cosGram[CostoBin(events[i][2])] +=1.0;
			
			//These lines account for containment effects.
			
			
		
			prob_NU = decayProb(in,Ud*chi,events_NU[i][3]);

			switch(which_var) {
				case 0:
					//std::cout<<"ebin: "<<EtoBin(events_NU[i][0])<<" prob: "<<prob_NU<<std::endl;
					Gram[EtoBin(events_NU[i][0])+NUstarter-1] += prob_NU;
					total_NU+=1.0;
					probTotal_NU+=prob_NU;
					break;
				case 1:
					Gram[CostoBin(cos((M_PI/180.0)*events_NU[i][1]))+NUstarter] += prob_NU;
					total_NU+=1.0;
					probTotal_NU+=prob_NU;
					break;
				case 2:
					qeWhichBin = QEtoBin(QEfromEandCos(events_NU[i][0],cos((M_PI/180.0)*events_NU[i][1])));
	
					if(qeWhichBin == 999) //Error Bin 
					{
						prob_NU=0;
						total_NU+=0.0;
					} else {
						Gram[qeWhichBin] += prob_NU;	
						total_NU+=1.0;
						probTotal_NU+=prob_NU;
					}
					break;

				default:							
					std::cout<<"ERROR: which_var in histogrammer must be 0(E), 1(A) or 2(QE)"<<std::endl;
			}



//		        printf("event %d: prob = %.5g\tprobTot = %.5g\n",i,prob,probTotal);

		}

		if(events_NUBAR[i][0]> 1e-4)
		{			
			//These lines DO NOT CONSIDER containment effects.
//			eGram[EtoBin(events[i][0])] +=1.0;
//			cosGram[CostoBin(events[i][2])] +=1.0;
			
			//These lines account for containment effects.
			
			
		
			prob_NUBAR = decayProb(in,Ud*chi,events_NUBAR[i][3]);

			switch(which_var) {
				case 0:
					Gram[EtoBin(events_NUBAR[i][0])+NUBARstarter-1] += prob_NUBAR;
					total_NUBAR+=1.0;
					probTotal_NUBAR+=prob_NUBAR;
					break;
				case 1:
					Gram[CostoBin(cos((M_PI/180.0)*events_NUBAR[i][1]))+NUBARstarter] += prob_NUBAR;
					total_NUBAR+=1.0;
					probTotal_NUBAR+=prob_NUBAR;
					break;
				case 2:
					qeWhichBin = QEtoBin(QEfromEandCos(events_NU[i][0],cos((M_PI/180.0)*events_NU[i][1])));
	
					if(qeWhichBin == 999) //Error Bin 
					{
						prob_NU=0;
						total_NU+=0.0;
					} else {
						Gram[qeWhichBin] += prob_NU;	
						total_NU+=1.0;
						probTotal_NU+=prob_NU;
					}
					break;

				default:							
					std::cout<<"ERROR: which_var in histogrammer must be 0(E), 1(A) or 2(QE)"<<std::endl;
			}



//		        printf("event %d: prob = %.5g\tprobTot = %.5g\n",i,prob,probTotal);

		}
	}


	for(int i=0;i<BINS/2;i++)
	{
		Gram[i]=Up*Up*chi*chi*cutEff_NU*finalScale_NU*Gram[i]/total_NU;
		Gram[i+NUBARstarter]=Up*Up*chi*chi*cutEff_NUBAR*finalScale_NUBAR*Gram[i+NUBARstarter]/total_NUBAR;
	}	

//printf("chiU: %.5g\tTot: %.5g\tprobTot/total: %.5g\tcosTot: %.5g\n",chiU, total, probTotal/total, cosTot);

return probTotal_NU/total_NU;
}


// ########################################### Onto Minimizer #############################



double nuisFuncE(const std::vector<double> &x, std::vector<double> &grad, void *my_data){
        
        nuisStruct *d = reinterpret_cast<nuisStruct*>(my_data);
        std::vector<double > eGram = d->egram;
        double sigma_zeta = d->Sigma_Zeta;

	int pos_selector = d->pos_selector;

        double zeta_b = x[0];
	std::vector<double > eO = {204, 280,214, 99, 83, 59, 51, 33, 37, 23, 19, 21, 12, 16, 4, 9, 4, 7, 3};
        std::vector<double > eB = {151.5, 218.8, 155.6, 108.7, 72.5, 57.6, 45, 38.5, 31.4,22.2, 20.4, 17.2, 14.1, 10.2, 9.1, 8.2, 5.6, 5.7, 2.9};


		if (pos_selector == POSNU || pos_selector == POSNUBAR){
				eO = {204, 280, 214, 99, 83, 59, 51, 33, 37, 23, 19, 21, 12, 16, 4, 9, 4, 7, 3};
				eB = {151.5, 218.8, 155.6, 108.7, 72.5, 57.6, 45, 38.5, 31.4,22.2, 20.4, 17.2, 14.1, 10.2, 9.1, 8.2, 5.6, 5.7, 2.9};
		} else if (pos_selector == NEGNU || pos_selector == NEGNUBAR){
				eO ={93,130,85,68,45,40,14,18,11,14,12,12,12,2,4,7,3,2,4};
				eB ={ 74.2,107.5,73.5,49.3,36.7,27.8,25.1,20.4,18.6,13.9,13.5,9.8,8.9,7.8,5.3,5,3.9,3.8,1.9};
		} 
		
                double temp_sig=0,temp_bg=0,lambda=0,N=0, E_N_events=0,E_N_sig_events=0,E_N_bg_events= 0,E_sum = 0;

        int bin = 0;
        for(bin=0;bin<EBINS;bin++)
                        {
                                temp_sig = eGram[bin];
                                temp_bg = (1.0+zeta_b)*eB[bin];
                                lambda = temp_sig + temp_bg;
                                N = eO[bin]; 
                                E_N_events += lambda;
                                E_N_sig_events += temp_sig;
                                E_N_bg_events += temp_bg;

                                E_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);

                        }
        E_sum+= pow((zeta_b/sigma_zeta),2.0);

return E_sum;
}


double nuisFuncA(const std::vector<double> &x, std::vector<double> &grad, void *my_data){
        
        nuisStruct *d = reinterpret_cast<nuisStruct*>(my_data);
        std::vector<double > cosGram = d->egram;
        std::vector<double > postGram(COSBINS,0);
int pos_selector = d->pos_selector;

        double sigma_zeta = d->Sigma_Zeta;
        double RunSigma = 1;
        //if(sigma_zeta == 0){RunSigma =0;};


        double zeta_b = x[0];
        
        std::vector<double > aO = {22,34,43,41,60,87,90,139,237,429};
        std::vector<double > aB = {19.9,23.1,28.8,32.1,46.4,63.1,86.1,121,196.8,390};
		if (pos_selector == POSNU || pos_selector == POSNUBAR){
				aO = {22,34,43,41,60,87,90,139,237,429};
				aB = {19.9,23.1,28.8,32.1,46.4,63.1,86.1,121,196.8,390};
			} else if (pos_selector == NEGNU || pos_selector == NEGNUBAR){
				aO = {10,13,16,20,24,36,41,70,94,263};
				aB = {9.2,11.2,13.5,16,18.7,24.2,36,52.1,94.9,237.1};
			} 
	

        double temp_sig=0,temp_bg=0,lambda=0,N=0, A_N_events=0,A_N_sig_events=0,A_N_bg_events= 0,A_sum = 0;
        double sigma_s = 1.0;
        int bin = 0;
        for(bin=0;bin<COSBINS;bin++)
                        {
                                temp_sig = sigma_s*cosGram[bin];
                                temp_bg = (1.0+zeta_b)*aB[bin];
                                lambda = temp_sig + temp_bg;
                                N = aO[bin]; //MB has seen O[] events.
                                
                                A_N_events += lambda;
                                A_N_sig_events += temp_sig;
                                A_N_bg_events += temp_bg;
                        //      E_sum+= (lambda-N)*(lambda-N)/lambda;
                                A_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
                                postGram[bin]=lambda;
                //              std::cout<<E_sum<<std::endl;
                        }
        A_sum+= pow((zeta_b/sigma_zeta),2.0);
        //std::cout<<std::setprecision(12)<<zeta_b<<"  "<<E_sum<<std::endl;
        //d->egram = postGram;

return A_sum;
}

double nuisFuncQE(const std::vector<double> &x, std::vector<double> &grad, void *my_data){
        
        nuisStruct *d = reinterpret_cast<nuisStruct*>(my_data);
        std::vector<double > qeGram = d->egram;
        double sigma_zeta = d->Sigma_Zeta;

	int pos_selector = d->pos_selector;

        double zeta_b = x[0];

        std::vector<double > qeO = {232,156,156,79,81,70,63,65,62,34,70};
        std::vector<double > qeB = {181.1,108.4,120.4,64.2,90.3,67.7,70.4,57.5,52.3,39,70.2};

        double temp_sig=0,temp_bg=0,lambda=0,N=0, QE_N_events=0,QE_N_sig_events=0,QE_N_bg_events= 0,QE_sum = 0;
        double sigma_s = 1.0;
        int bin = 0;
        for(bin=0;bin<QEBINS;bin++)
                        {
                                temp_sig = sigma_s*qeGram[bin];
                                temp_bg = (1.0+zeta_b)*qeB[bin];
                                lambda = temp_sig + temp_bg;
                                N = qeO[bin]; //MB has seen O[] events.
                                
                                QE_N_events += lambda;
                                QE_N_sig_events += temp_sig;
                                QE_N_bg_events += temp_bg;
                        //      E_sum+= (lambda-N)*(lambda-N)/lambda;
                                QE_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
                        }
        QE_sum+= pow((zeta_b/sigma_zeta),2.0);
        //std::cout<<std::setprecision(12)<<zeta_b<<"  "<<E_sum<<std::endl;

return QE_sum;
}



double nuisMarginalize(std::vector<double > * bf_zeta_b, double * chi, std::vector<double > * eVGram, int whi, double SIGMAZETA, int pos_selector){
        
        nuisStruct ddata = {(*eVGram),SIGMAZETA,pos_selector};


        nlopt::opt full(nlopt::LN_COBYLA,1);
        
        //full.set_maxeval(200000);
        full.set_xtol_abs(1e-5);
        full.set_lower_bounds(-1);
        full.set_upper_bounds(1);
        std::vector<double > xtem = {0.01};
        double ctem;
        if(whi==0){ 
                full.set_min_objective(nuisFuncE, &ddata);
        } else if(whi==1){
                full.set_min_objective(nuisFuncA,&ddata);
        } else if(whi==2){
		full.set_min_objective(nuisFuncQE,&ddata);
	} else {
		std::cout<<"ERROR in LR.c ~ nuisMarginalise: spectrum selector which is > 2"<<std::endl;
	}


        
        //std::cout<<"Starting GLOBAL_FULL"<<std::endl;
        nlopt::result result = full.optimize(xtem, ctem);
        //std::cout<<"Minimized GLOBAL_FULL"<<std::endl;
        //std::cout<<"GLOBAl_FULL: found minimum at "<<xtem[0]<<" and thats "<<ctem<<std::endl;

        (*bf_zeta_b)=xtem;
        (*chi)=ctem;
        (*eVGram) = ddata.egram;

return 0;
}



double nuisFuncE_dual(const std::vector<double> &x, std::vector<double> &grad, void *my_data){
        
        nuisStruct *d = reinterpret_cast<nuisStruct*>(my_data);
        std::vector<double > eGram = d->egram;
        double sigma_zeta = d->Sigma_Zeta;


        double zeta_b = x[0];
	std::vector<double > eO = {204, 280, 214, 99, 83, 59, 51, 33, 37, 23, 19, 21, 12, 16, 4, 9, 4, 7, 3,93,130,85,68,45,40,14,18,11,14,12,12,12,2,4,7,3,2,4};
        std::vector<double > eB = {151.5, 218.8, 155.6, 108.7, 72.5, 57.6, 45, 38.5, 31.4,22.2, 20.4, 17.2, 14.1, 10.2, 9.1, 8.2, 5.6, 5.7, 2.9, 74.2,107.5,73.5,49.3,36.7,27.8,25.1,20.4,18.6,13.9,13.5,9.8,8.9,7.8,5.3,5,3.9,3.8,1.9};

	
                double temp_sig=0,temp_bg=0,lambda=0,N=0, E_N_events=0,E_N_sig_events=0,E_N_bg_events= 0,E_sum = 0;

        int bin = 0;
        for(bin=0;bin<2*EBINS;bin++)
                        {
                                temp_sig = eGram[bin];
                                temp_bg = (1.0+zeta_b)*eB[bin];
                                lambda = temp_sig + temp_bg;
                                N = eO[bin]; 
                                E_N_events += lambda;
                                E_N_sig_events += temp_sig;
                                E_N_bg_events += temp_bg;

                                E_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);

                        }
        E_sum+= pow((zeta_b/sigma_zeta),2.0);

return E_sum;
}


double nuisFuncA_dual(const std::vector<double> &x, std::vector<double> &grad, void *my_data){
        
        nuisStruct *d = reinterpret_cast<nuisStruct*>(my_data);
        std::vector<double > cosGram = d->egram;

        double sigma_zeta = d->Sigma_Zeta;
        double RunSigma = 1;
        //if(sigma_zeta == 0){RunSigma =0;};


        double zeta_b = x[0];
        
	std::vector<double > aO = {22,34,43,41,60,87,90,139,237,429,10,13,16,20,24,36,41,70,94,263};
	std::vector<double > aB  = {19.9,23.1,28.8,32.1,46.4,63.1,86.1,121,196.8,390,9.2,11.2,13.5,16,18.7,24.2,36,52.1,94.9,237.1};


        double temp_sig=0,temp_bg=0,lambda=0,N=0, A_N_events=0,A_N_sig_events=0,A_N_bg_events= 0,A_sum = 0;
        double sigma_s = 1.0;
        int bin = 0;
        for(bin=0;bin<2*COSBINS;bin++)
                        {
                                temp_sig = sigma_s*cosGram[bin];
                                temp_bg = (1.0+zeta_b)*aB[bin];
                                lambda = temp_sig + temp_bg;
                                N = aO[bin]; //MB has seen O[] events.
                                
                                A_N_events += lambda;
                                A_N_sig_events += temp_sig;
                                A_N_bg_events += temp_bg;
                        //      E_sum+= (lambda-N)*(lambda-N)/lambda;
                                A_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
                //              std::cout<<E_sum<<std::endl;
                        }
        A_sum+= pow((zeta_b/sigma_zeta),2.0);
        //std::cout<<std::setprecision(12)<<zeta_b<<"  "<<E_sum<<std::endl;
        //d->egram = postGram;

return A_sum;
}

double nuisFuncQE_dual(const std::vector<double> &x, std::vector<double> &grad, void *my_data){
        
        nuisStruct *d = reinterpret_cast<nuisStruct*>(my_data);
        std::vector<double > qeGram = d->egram;
        double sigma_zeta = d->Sigma_Zeta;


        double zeta_b = x[0];


        std::vector<double > qeO = {232,156,156,79,81,70,63,65,62,34,70};
        std::vector<double > qeB = {181.1,108.4,120.4,64.2,90.3,67.7,70.4,57.5,52.3,39,70.2};

        double temp_sig=0,temp_bg=0,lambda=0,N=0, QE_N_events=0,QE_N_sig_events=0,QE_N_bg_events= 0,QE_sum = 0;
        double sigma_s = 1.0;
        int bin = 0;
        for(bin=0;bin<QEBINS;bin++)
                        {
                                temp_sig = sigma_s*qeGram[bin];
                                temp_bg = (1.0+zeta_b)*qeB[bin];
                                lambda = temp_sig + temp_bg;
                                N = qeO[bin]; //MB has seen O[] events.
                                
                                QE_N_events += lambda;
                                QE_N_sig_events += temp_sig;
                                QE_N_bg_events += temp_bg;
                        //      E_sum+= (lambda-N)*(lambda-N)/lambda;
                                QE_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
                        }
        QE_sum+= pow((zeta_b/sigma_zeta),2.0);
        //std::cout<<std::setprecision(12)<<zeta_b<<"  "<<E_sum<<std::endl;

return QE_sum;
}
double nuisMarginalize_dual(std::vector<double > * bf_zeta_b, double * chi, std::vector<double > * eVGram, int whi, double SIGMAZETA){
        
        nuisStruct ddata = {(*eVGram),SIGMAZETA};


        nlopt::opt full(nlopt::LN_COBYLA,1);
        
        //full.set_maxeval(200000);
        full.set_xtol_abs(1e-5);
        full.set_lower_bounds(-1);
        full.set_upper_bounds(1);
        std::vector<double > xtem = {0.01};
        double ctem;
        if(whi==0){ 
                full.set_min_objective(nuisFuncE_dual, &ddata);
        } else if(whi==1){
                full.set_min_objective(nuisFuncA_dual,&ddata);
        } else if(whi==2){
		full.set_min_objective(nuisFuncQE_dual,&ddata);
	} else {
		std::cout<<"ERROR in LR.c ~ nuisMarginalise: spectrum selector which is > 2"<<std::endl;
	}


        
        //std::cout<<"Starting GLOBAL_FULL"<<std::endl;
        nlopt::result result = full.optimize(xtem, ctem);
        //std::cout<<"Minimized GLOBAL_FULL"<<std::endl;
        //std::cout<<"GLOBAl_FULL: found minimum at "<<xtem[0]<<" and thats "<<ctem<<std::endl;

        (*bf_zeta_b)=xtem;
        (*chi)=ctem;
        (*eVGram) = ddata.egram;

return 0;
}

