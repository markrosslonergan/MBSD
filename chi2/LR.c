#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <nlopt.hpp>
#include <iomanip>


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

double nuisFuncE(const std::vector<double> &x, std::vector<double> &grad, void *my_data){
        
        nuisStruct *d = reinterpret_cast<nuisStruct*>(my_data);
        std::vector<double > eGram = d->egram;
        double sigma_zeta = d->Sigma_Zeta;


        double zeta_b = x[0];
                                 //{204, 280, 214, 99, 83, 59, 51, 33, 37, 23, 19, 21, 12, 16, 4, 9, 4, 7, 3}
        std::vector<double > eO = {204.0, 280.0,214.0, 99.0, 83.0, 59.0, 51.0, 33.0, 37.0, 23.0, 19.0, 21.0, 12.0, 16.0, 4.0, 9.0, 4.0, 7.0, 3.0};
        std::vector<double > eB = {151.5, 218.8, 155.6, 108.7, 72.5, 57.6, 45, 38.5, 31.4,22.2, 20.4, 17.2, 14.1, 10.2, 9.1, 8.2, 5.6, 5.7, 2.9};

        double temp_sig=0,temp_bg=0,lambda=0,N=0, E_N_events=0,E_N_sig_events=0,E_N_bg_events= 0,E_sum = 0;
        double sigma_s = 1.0;
        int bin = 0;
        for(bin=0;bin<EBINS;bin++)
                        {
                                temp_sig = sigma_s*eGram[bin];
                                temp_bg = (1.0+zeta_b)*eB[bin];
                                lambda = temp_sig + temp_bg;
                                N = eO[bin]; //MB has seen O[] events.
                                
                                E_N_events += lambda;
                                E_N_sig_events += temp_sig;
                                E_N_bg_events += temp_bg;
                        //      E_sum+= (lambda-N)*(lambda-N)/lambda;
                                E_sum+= 2.0*(lambda-N) + 2.0*N*log(N/lambda);
                //              std::cout<<E_sum<<std::endl;
                        }
        E_sum+= pow((zeta_b/sigma_zeta),2.0);
        //std::cout<<std::setprecision(12)<<zeta_b<<"  "<<E_sum<<std::endl;

return E_sum;
}


double nuisFuncA(const std::vector<double> &x, std::vector<double> &grad, void *my_data){
        
        nuisStruct *d = reinterpret_cast<nuisStruct*>(my_data);
        std::vector<double > cosGram = d->egram;
        std::vector<double > postGram(COSBINS,0);

        double sigma_zeta = d->Sigma_Zeta;
        double RunSigma = 1;
        //if(sigma_zeta == 0){RunSigma =0;};


        double zeta_b = x[0];
        
        std::vector<double > aO = {22,34,43,41,60,87,90,139,237,429};
        std::vector<double > aB = {19.9,23.1,28.8,32.1,46.4,63.1,86.1,121,196.8,390};

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
        A_sum+= RunSigma*pow((zeta_b/sigma_zeta),2.0);
        //std::cout<<std::setprecision(12)<<zeta_b<<"  "<<E_sum<<std::endl;
        //d->egram = postGram;

return A_sum;
}


double boundChiU(double ms, double mz)
{
return boundChi(mz)*boundU(ms);
}

double nuisMarginalize(std::vector<double > * bf_zeta_b, double * chi, std::vector<double > * eVGram, int whi, double SIGMAZETA){
        
        nuisStruct ddata = {(*eVGram),SIGMAZETA};

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

