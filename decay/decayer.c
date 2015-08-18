#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define NUMEVENTS 200000

struct PDF_CHOICE { double Enu; double cosThnu; double Phinu; };
typedef struct OBSERVABLES { double E_sum; double Th_sum; double AngSep; double E_sterile; double Th_sterile; double E_high; double Th_high; double E_low; double Th_low; double FS_AngSep; } OBSERVABLES;

double pdf_function_test(double x,double y,double mS, double mZprime, void * pointer)
{
	
	double ret = (3.0/2.0)*(1.0 - x*x);
	if(ret<0){ ret = 0.0; }

return ret; 
}

double pdf_function(double x, double y, double mS, double mZprime, void * pointer)
{
	double mu_s  = mS/mZprime;
	double alpha = mu_s*mu_s/(1.0-mu_s*mu_s);

	double invnorm;

	if(alpha < 0.01)
	{
		double invnorm_perturb_0 = (1.0/(1.0+alpha))*(7.0/4.0 + 41.0/60.0*alpha); 

		double invnorm_perturb_rest = -0.18333*pow(alpha,2.0)+0.22857*pow(alpha,3.0)-0.23274*pow(alpha,4.0)+0.22421*pow(alpha,5.0)-0.21190*pow(alpha,6.0)+0.19899*pow(alpha,7.0)-0.18662*pow(alpha,8.0)+0.17517*pow(alpha,9.0)-0.16475*pow(alpha,10.0)+0.15531*pow(alpha,11.0);
		invnorm = invnorm_perturb_0+invnorm_perturb_rest;
	}
	else 
	{		
		invnorm = (3.0/(2.0*pow(alpha,3.0)))*(2+alpha-3*pow(alpha,2.0))/(1.0+alpha) + (4.0*alpha*alpha-3.0)*(log(1+alpha)/log(exp(1.0)))/pow(alpha,4.0);
	}

	double ret = (1.0/invnorm)*x*(4-x*x)/((1.0+alpha*x)*(1.0+alpha*x));

//	printf("inv. norm: %.5lf\talpha: %.5lf\tret: %.5lf\n",invnorm,alpha,ret);

	if(ret<0){ ret = 0.0; }

return ret; 
}

double plotPdf(double mS, double mZprime, double (*pdf_function)(double, double, double, double, void *))
{
	double x;
	for(x=0.0;x<1.0+1e-5;x+=1e-2)
	{
		printf("%.5lf %.5lf %.5g\n",x,pdf_function(x,0.0,mS,mZprime,NULL),mS*mS/(mZprime*mZprime - mS*mS));

	}

return 0.0;
}

struct PDF_CHOICE choose_from_pdf(gsl_rng * r,double mS, double mZprime, double (*pdf)(double, double, double, double, void *))
{
	double mu_s  = mS/mZprime;
	double alpha = mu_s*mu_s/(1-mu_s*mu_s);

//	double PDF_MAX = 3.0/2.0; //pdf_function_test
	double PDF_MAX = 1.8; //pdf_function

	double x = gsl_rng_uniform(r);
	double y = -1.0 + 2.0*gsl_rng_uniform(r);
	double phi = 2*M_PI*gsl_rng_uniform(r);
	double z = (PDF_MAX+0.01)*gsl_rng_uniform(r);

	while(pdf(x,y,mS,mZprime,NULL)<=z)
	{
//		printf("I tried!\n");
		x = gsl_rng_uniform(r);
		y = -1.0 + 2.0*gsl_rng_uniform(r);
		z = (PDF_MAX+0.01)*gsl_rng_uniform(r);
	}

//	printf("I succeeded!\n");

	//printf("%.5lf %.5lf %.5lf %.5lf\n",x,y,z,pdf(x,y,mS,mZprime,NULL));
	struct PDF_CHOICE package;
	package.Enu = mS*x/2.0;
	package.cosThnu = y;
	package.Phinu = phi;

return package;
}

double getEvents(double mS, double mZprime, double events[][2])
{
	FILE *ptr_file;
    	
	char buf[3000];

	char * pch;

	int n = 1;
	int m = 0;
	char s[100];
	char filename[500] = "../MC/MC_\0";
	sprintf(s,"%.4lf_%.4lf.dat", mS, mZprime);
	strcat(filename,s);
//	printf("Filename: %s\n",filename);
	ptr_file =fopen(filename,"r");

    	if (!ptr_file)
       	{			
		printf("ERROR LOADING MC EVENTS\n");
		exit(1);
	}
    	while (fgets(buf,3000, ptr_file)!=NULL)
	{
		pch = strtok(buf,"\t");
		n=1;
 		while (pch != NULL)
		{
			if(n==3){	//printf("%.7g, ",strtof(pch,NULL));
					events[m][0] = strtof(pch,NULL);	 
				}
			if(n==4){	//printf("%.7g\n",strtof(pch,NULL));
					events[m][1] = strtof(pch,NULL);	 
				}
			pch = strtok(NULL,"\t");
		n++;	
		}
		m++;
	}
	fclose(ptr_file);

//	printf("%s\n",flux_temp);

//printf("Total lines: %d\n", m);

return 0;
}


int printEvents(double events[][2])
{
	int i;
	for (i=0;i<=NUMEVENTS-1;i++)
	{ 
		printf("%.9g\t%.9g\n",events[i][0],events[i][1]);
	}
return 0;
}


int drawRestFrameDist(gsl_rng * r, double mS, double mZprime, double output[3])
{
	struct PDF_CHOICE choice=choose_from_pdf(r,mS,mZprime,pdf_function);
	output[0]=choice.Enu;
	output[1]=choice.cosThnu;
	output[2]=choice.Phinu;
return 0;
}

int computeLabFrameVariables(OBSERVABLES * output, double mS, double Es, double costhS, double phiS, double restFrameParams[3])
{
	double Enu = restFrameParams[0];
	double Th = acos(restFrameParams[1]);
	double Phi = restFrameParams[2];

	double me = 0.0;// THIS CAUSES ERRORS! 5.11e-6; //GeV
	
	double Ee = (mS - Enu)/2.0;
	double Pe = sqrt(Ee*Ee-me*me);
	double beta = sqrt(1-mS*mS/(Es*Es));
	double gamma = 1.0/sqrt(1.0-beta*beta);

//printf("%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n", Enu, Ee, Pe, beta, gamma, mS); 

	double alpha = 2.0*acos( Enu/(2.0*Pe) );
	double theta_plus = M_PI - Th - alpha/2.0;
	double theta_minus = M_PI - Th + alpha/2.0;

	double Pplus_E = gamma*(Ee + beta*Pe*cos(theta_plus));
	double Pminus_E = gamma*(Ee + beta*Pe*cos(theta_minus));
	double Pplus_x = Pe*sin(theta_plus)*cos(Phi);
	double Pminus_x = Pe*sin(theta_minus)*cos(Phi);
	double Pplus_y = Pe*sin(theta_plus)*sin(Phi);
	double Pminus_y = Pe*sin(theta_minus)*sin(Phi);
	double Pplus_z = gamma*(Pe*cos(theta_plus) + beta*Ee);
	double Pminus_z = gamma*(Pe*cos(theta_minus) + beta*Ee);

	double Pee[] = {(Pplus_x + Pminus_x)/2.0, (Pplus_y + Pminus_y)/2.0, (Pplus_z + Pminus_z)/2.0};
	double Pplus[] = {Pplus_x, Pplus_y, Pplus_z};
	double Pminus[] = {Pminus_x, Pminus_y, Pminus_z};
	rotor(acos(costhS),phiS,Pee);
	rotor(acos(costhS),phiS,Pplus);
	rotor(acos(costhS),phiS,Pminus);


	output->E_sum = Pplus_E + Pminus_E; // energy is unaffected by rotation.
	output->Th_sum = (180.0/M_PI)*acos(Pee[2]/sqrt(Pee[0]*Pee[0] + Pee[1]*Pee[1] + Pee[2]*Pee[2] )); 
	output->AngSep = (180.0/M_PI)*acos((Pplus_x*Pminus_x + Pplus_y*Pminus_y + Pplus_z*Pminus_z)/(sqrt(Pplus_x*Pplus_x + Pplus_y*Pplus_y + Pplus_z*Pplus_z)*sqrt(Pminus_x*Pminus_x + Pminus_y*Pminus_y + Pminus_z*Pminus_z))); // opening angle is unaffected by rotation


//OLD BROKEN FS_ANGSEP.
//	if(Pminus[2] > 0 && Pplus[2] > 0)
//	{
//		output->FS_AngSep = (180.0/M_PI)*fabs(atan(Pplus[0]/Pplus[2]) - atan(Pminus[0]/Pminus[2]));
//	}
//	else 
//	{
//		output->FS_AngSep = 180- (180.0/M_PI)*fabs(atan(Pplus[0]/Pplus[2]) - atan(Pminus[0]/Pminus[2]));
//	}


	if(Pminus[2]*Pplus[2] >= 0 && Pplus[0]*Pplus[0] >= 0)
	{
		output->FS_AngSep = (180.0/M_PI)*fabs(fabs(atan(Pminus[0]/Pminus[2])) - fabs(atan(Pplus[0]/Pplus[2])));
	}
	else if(Pminus[2]*Pplus[2] >= 0 && Pplus[0]*Pplus[0] < 0)
	{
		output->FS_AngSep = (180.0/M_PI)*fabs(atan(Pminus[0]/Pminus[2]) + atan(Pplus[0]/Pplus[2]));
	}
	else if(Pminus[2]*Pplus[2] < 0 && Pplus[0]*Pplus[0] >= 0)
	{
		output->FS_AngSep = 180.0- (180.0/M_PI)*fabs(atan(Pminus[0]/Pminus[2]) + atan(Pplus[0]/Pplus[2]));
	}
	else if(Pminus[2]*Pplus[2] < 0 && Pplus[0]*Pplus[0] < 0)
	{
		output->FS_AngSep = 180.0- (180.0/M_PI)*fabs(atan(Pminus[0]/Pminus[2]) - atan(Pplus[0]/Pplus[2]));
	}




	if(Pplus_E < Pminus_E)
	{
		output->E_low = Pplus_E;  // E_low
		output->E_high = Pminus_E; // E_high

		output->Th_low = (180.0/M_PI)*acos(Pplus[2]/sqrt(Pplus[0]*Pplus[0] +Pplus[1]*Pplus[1] +Pplus[2]*Pplus[2]));
		output->Th_high = (180.0/M_PI)*acos(Pminus[2]/sqrt(Pminus[0]*Pminus[0] +Pminus[1]*Pminus[1] +Pminus[2]*Pminus[2]));
	}
	else 
	{
		output->E_low = Pminus_E; // E_low
		output->E_high = Pplus_E;  // E_high

		output->Th_low = (180.0/M_PI)*acos(Pminus[2]/sqrt(Pminus[0]*Pminus[0] +Pminus[1]*Pminus[1] +Pminus[2]*Pminus[2]));
		output->Th_high = (180.0/M_PI)*acos(Pplus[2]/sqrt(Pplus[0]*Pplus[0] +Pplus[1]*Pplus[1] +Pplus[2]*Pplus[2]));
	}

return 0;
}

int rotor(double theta, double phi, double vector[3])
{
	double x=vector[0];
	double y=vector[1];
	double z=vector[2];

	double rdotn1 = x*cos(phi) + y*sin(phi);
	double rdotn2 = z;
	double rdotn3 = -x*sin(phi) + y*cos(phi);

	vector[0]=(cos(theta)*rdotn1 + sin(theta)*rdotn2)*cos(phi) - rdotn3*sin(phi);
	vector[1]=(cos(theta)*rdotn1 + sin(theta)*rdotn2)*sin(phi) + rdotn3*cos(phi);
	vector[2]=-sin(theta)*rdotn1 + cos(theta)*rdotn2;

return 0;
}

/* ########## Main function ############### */
int main(int argc, char* argv[])
{

double mS = strtof(argv[1],NULL);
double mZprime = strtof(argv[2],NULL); 

static double events[NUMEVENTS][2]; //define the storage for all the events, [0] = E_s, [1] =cos\theta_s
getEvents(mS,mZprime,events); 

//plotPdf(mS,mZprime,pdf_function);
//printEvents(events);

//We enter the main loop over events. For each one, computing the relevant observables.
const gsl_rng_type * T;
gsl_rng * r;
gsl_rng_env_setup();

T = gsl_rng_default;
r = gsl_rng_alloc (T);

static OBSERVABLES Obs;
double phiS=0.0;
double restFrameParams[3];

int m;
for(m=0;m<=NUMEVENTS-1;m++)
{
	phiS = 2.0*M_PI*gsl_rng_uniform(r);

	//These are drawn from the rest frame decay.
	drawRestFrameDist(r,mS,mZprime,restFrameParams); //this will populate the doubles.

	//Below: populates Obs.
	computeLabFrameVariables(&Obs,mS,events[m][0],events[m][1],phiS,restFrameParams);

	Obs.E_sterile = events[m][0];
	Obs.Th_sterile = events[m][1];

	// E_sum / Th_sum / AngSep. / E_sterile / Th_sterile / E_high / Th_high/ E_low / Th_low
	printf("%.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf\n",Obs.E_sum, Obs.Th_sum, Obs.AngSep, Obs.E_sterile, Obs.Th_sterile, Obs.E_high, Obs.Th_high, Obs.E_low, Obs.Th_low, Obs.FS_AngSep);

}

//plotPdf(mS,mZprime,pdf_function_test);

/*	double vec[] = {0.0,0.0,1.0};
	double test_theta;
	//for(test_theta = 0.0; test_theta<2*M_PI; test_theta+=M_PI/100)
	//{
		rotor(M_PI/2,M_PI/2,vec);
		printf("%.5lf %.5lf %.5lf\n", vec[0], vec[1], vec[2]);
	//}
*/
	
gsl_rng_free(r);

return 0;
}
