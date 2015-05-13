#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "LR.h"
//#include "Minuit2/FCNBase.h"

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

double histogrammer(CL_input in, double chiU, double cutEff, const double events[][NUM_EVENT_OBS], double eGram[], double cosGram[])
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
			cosGram[CostoBin(cos((M_PI/180.0)*events[i][1]))] += prob;
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

double boundU( double ms){
        double peakB=1;
        double fit1 = -1.7969808082063344*pow(10,-8); 
        double fit2 = 0.6078622816783561;
        

        if(ms < 0.003){
                peakB = 0.3*exp(-3000*ms);
        } else if (ms>=0.003 && ms<=0.03){
                peakB = 5*pow(10,-5)*exp(-100*ms)+fit1;
        } else if (ms>0.03){
               peakB = pow(10,-5)*exp(-30*ms)*fit2;
        }              


return sqrt(peakB);
}

double boundChi(double mz){
        double babarB=1;
        
        if(mz<0.25){
                babarB=-0.001*0.9987312619554833*pow(mz - 0.1,0.284) + 0.001;
        } else if(mz>=0.25 && mz<0.52){
                babarB=0.000392*exp(mz*mz);
        } else if(mz>=0.52){
                babarB=0.0012*pow(mz,1.3) + pow(8.671375207584429,-7);
        }

return babarB;

}


double boundChiU(double ms, double mz)
{
return boundChi(mz)*boundU(ms);
}


