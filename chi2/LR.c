#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "LR.h"
#include "Minuit2/FCNBase.h"

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
	char filename[500] = "../MC/HIST_\0";
	sprintf(s,"%.4lf_%.4lf.dat", mS, mZprime);
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

double ROOT::Minuit2::E_log_likelihood::operator()(const std::vector<double>& p) const
	{
		CL_input in;
		in.mS = p[0];
		in.mZprime = p[1];
		double logchiU = p[2];
		double chiU = pow(10.0,logchiU);
		double zeta_b = p[3];
		double sigma_s = p[4];
		double cutEff = p[5];

//		std::cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<" "<<p[3]<<" "<<p[4]<<" "<<p[5]<<std::endl;

		double N=0.0;
		double lambda=0.0;
		double sum=0.0;

		double eGram[EBINS];
		double cosGram[COSBINS];

		double O[EBINS];
		double B[EBINS];

		//Ebinlow(GeV) 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9
		O[0] =  204.0;
		O[1] = 280.0; 
		O[2] = 214.0; 
		O[3] = 99.0; 
		O[4] = 83.0; 
		O[5] = 59.0; 
		O[6] = 51.0; 
		O[7] = 33.0; 
		O[8] = 37.0; 
		O[9] = 23.0; 	
		O[10] = 19.0;
 		O[11] = 21.0; 
		O[12] = 12.0; 
		O[13] = 16.0; 		
		O[14] = 4.0; 
		O[15] = 9.0; 
		O[16] = 4.0; 
		O[17] = 7.0; 
		O[18] = 3.0;

		B[0] = 151.5;
		B[1] = 218.8;
		B[2] = 155.6;
		B[3] = 108.7;
		B[4] = 72.5;
		B[5] = 57.6;
		B[6] = 45;
		B[7] = 38.5;
		B[8] = 31.4;
		B[9] = 22.2;
		B[10] = 20.4;
		B[11] = 17.2;
		B[12] = 14.1;
		B[13] = 10.2;
		B[14] = 9.1;
		B[15] = 8.2;
		B[16] = 5.6;
		B[17] = 5.7;
		B[18] = 2.9;
		

		histogrammer(in,chiU,cutEff,events,eGram,cosGram);

		int bin = 0;
		if (sigma_s > 1e-5){
			for(bin=0;bin<EBINS;bin++)
			{
				lambda = sigma_s*eGram[bin] + (1.0+zeta_b)*B[bin];
				N = O[bin]; //MB has seen O[] events.
				
			//	sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
				sum+= (lambda-N)*(lambda-N)/lambda;
			}
		}

		//a guesstimate of the systematic error on the background.
		double sigma_zeta = 0.05;
		// add prior on zeta.
		sum += pow((zeta_b/sigma_zeta),2.0);
		return sum;
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

double decayProb(CL_input input, double chiU, double Es)
{
	double mS = input.mS;
	double mZprime = input.mZprime;

	double L = 1.0; // 1 metre.
	double g1 = 0.2;

	double alpha = mS*mS/(mZprime*mZprime-mS*mS);

	double prefac = (chiU*chiU*(1.0/(137.05*137.05))*alpha*alpha*mS/(4.0*pow(4.0*M_PI,3.0)*g1*g1));
	//double func = (3.0/(2.0*alpha*alpha*alpha))*((2+alpha-3*alpha*alpha)/(1+alpha))+(4*alpha*alpha-3)*log(1.0+alpha)/pow(alpha,4.0);
	double func;
	if(alpha < 0.01)
	{
	
		double func_perturb_0 = (1.0/(1.0+alpha))*(7.0/4.0 + 41.0/60.0*alpha); 

		double func_perturb_rest = -0.18333*pow(alpha,2.0)+0.22857*pow(alpha,3.0)-0.23274*pow(alpha,4.0)+0.22421*pow(alpha,5.0)-0.21190*pow(alpha,6.0)+0.19899*pow(alpha,7.0)-0.18662*pow(alpha,8.0)+0.17517*pow(alpha,9.0)-0.16475*pow(alpha,10.0)+0.15531*pow(alpha,11.0);
		func = func_perturb_0+func_perturb_rest;
	}
	else 
	{		
		func = (3.0/(2.0*pow(alpha,3.0)))*(2+alpha-3*pow(alpha,2.0))/(1.0+alpha) + (4.0*alpha*alpha-3.0)*(log(1+alpha)/log(exp(1.0)))/pow(alpha,4.0);
	}


//printf("### mS: %.5lf\tmZ: %.5lf\talpha: %.5lf\tprefac: %.5lf\tfunc: %.5lf\n",mS,mZprime,alpha,prefac,func); 

return 1.0-exp(-(0.51e16)*prefac*func*L*mS/Es);
} 

double histogrammer(CL_input in, double chiU, double cutEff, const double events[][4], double eGram[], double cosGram[])
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

	double total = 0.0;
	double probTotal = 0.0;
	double prob = 0.0;
	for(i=0;i<NUMEVENTS;i++)
	{
		if(events[i][0] > 1e-4)
		{			
			//These lines DO NOT CONSIDER containment effects.
//			eGram[EtoBin(events[i][0])] +=1.0;
//			cosGram[CostoBin(events[i][2])] +=1.0;
			
			//These lines account for containment effects.
			prob = decayProb(in,chiU,events[i][3]);
			eGram[EtoBin(events[i][0])] += prob;
			cosGram[CostoBin(events[i][2])] += prob;
//			
			total+=1.0;
			probTotal+=prob;
//		printf("event %d: prob = %.5g\tprobTot = %.5g\n",i,prob,probTotal);

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

//printf("chiU: %.5g\tTot: %.5g\tprobTot/total: %.5g\tcosTot: %.5g\n",chiU, total, probTotal/total, cosTot);

return probTotal/total;
}


