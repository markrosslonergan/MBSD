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
	return ans1-(c1+c2*Q2+c3*Q2*Q2+c4*pow(Q2,3));
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

double histogrammer_indiv(CL_input in, double chiU, double cutEff, const double events[][NUM_EVENT_OBS], double Gram[], int which_var)
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

	for(i=0;i<NUMEVENTS;i++)
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
double boundUpeaky( double ms, double mz){
        std::vector<double  > peakBounds = {0.01238242210078086,0.0012535309071088556,0.00004547996544027682,0.00003486400316563211,0.000048498500851545515,0.000043481154180279925,0.00003814588076274375,0.000032810607345207566,   0.00002747533392767138,0.000022140060510135193,0.00001848513484898103,0.000016648710422498732,0.000014812285996016437,0.000012975861569534138,0.000011139437143051847,9.30301271656955e-6,   0.000011722952207256827,0.000022348900710157802,0.000014199121671679938,0.000015074668173337534,0.000020663287575263602,0.000025736549869991898,0.000020780763938661232,   0.000012406536764883298,9.2148124995948e-6,0.000010398349461906154,0.000013692168512906922,0.000010611689039650524,0.000014181454464076391,0.00003452605035229107,0.00009980899398935939,   0.00010519659404478175,0.0001287558277072979,0.002567828076615406,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.06000974457459983,0.0005491388780786994,   0.00008692580899813938,0.00007546678867820339,0.0000670313055282429,0.00006149850811247785,0.000056060882219892045,0.000051026402899493555,0.00004599192357909507,0.00004104215865575888,   0.00003610485398698345,0.00003213214118430577,0.000028244726729530042,0.000025169261508564636,0.000022141273340815894,0.000020562787108851245,0.0000189843008768866,0.000017852308552278194,   0.000016917644673245147,0.000015983782288916364,0.000015148034253223484,0.000014312286217530604,0.000013505417644784022,0.00001275916261975577,0.000012012907594727518,   0.000011338142722287882,0.000010704064528463935,0.00001006998633463999,9.526190122643592e-6,8.986956234247341e-6,8.467789113472166e-6,8.003854169017238e-6,7.53991922456231e-6,   7.160113098379005e-6,6.840955099379955e-6,6.5217971003809065e-6,6.3318009105642275e-6,6.156846509179875e-6,5.985745424164471e-6,5.829114967200497e-6,5.672484510236525e-6,   5.521599525386795e-6,5.37591519009579e-6,5.230230854804784e-6,5.087619257316122e-6,4.9455951893717475e-6,4.805406514441987e-6,4.6752942018612466e-6,4.545181889280506e-6,   4.417519134905032e-6,4.2926333217711365e-6,4.167747508637242e-6,4.049179574916104e-6,3.9323597050098456e-6,3.816522894308477e-6,3.709746907510374e-6,3.6029709207122713e-6,   3.499494812756753e-6,3.4007190291538117e-6,3.30194324555087e-6,3.2070640684738203e-6,3.1136475954414644e-6,3.020568557785226e-6,2.9352498238875896e-6,2.849931089989953e-6,   2.7672257597970884e-6,2.6892425193466815e-6,2.611259278896275e-6,2.537195845154522e-6,2.4650555080396757e-6,2.3929151709248293e-6,2.330525315023786e-6,2.268279280424445e-6,   2.214568314856435e-6,2.180766612833341e-6,2.1469649108102475e-6,2.4198725176186026e-6,2.7216475579144377e-6,2.0966538412613143e-6,2.0253156984318495e-6,1.9961258269669026e-6,   1.9673640166810762e-6,1.9399295294714553e-6,1.912495042261834e-6,1.8852879607601946e-6,1.8582606792909906e-6,1.8312333978217863e-6,1.8051073001928576e-6,1.7791122410725753e-6,   1.7534105274634778e-6,1.7289787269396175e-6,1.7045469264157571e-6,1.6802972597673258e-6,1.6562280965818638e-6,1.632158933396402e-6,1.6088406528305205e-6,1.5856907359328586e-6,   1.5627238101488947e-6,1.5409660542164883e-6,1.519208298284082e-6,1.4975945902102746e-6,1.4761597811967042e-6,1.4547249721831336e-6,1.4339128498842612e-6,1.4132966762533873e-6,   1.3929102467418403e-6,1.3753476443561706e-6,1.357785041970501e-6,1.3399414361535207e-6,1.3216580086911777e-6,1.3033745812288345e-6,1.2870494525040906e-6,1.2715459864491707e-6,   1.2560363895453765e-6,1.2402107856410494e-6,1.2243851817367222e-6,1.208606292505019e-6,1.1929207081162453e-6,1.1772351237274717e-6,1.1619023787167372e-6,1.1467611082781001e-6,   1.1316198378394632e-6,1.1173391202070031e-6,1.103091628191469e-6,1.0891755749541034e-6,1.0761227856741227e-6,1.0630699963941423e-6,1.049245380636674e-6,1.0348888779993191e-6,   1.0205323753619645e-6,1.0078534419944996e-6,9.953472956818846e-7,9.828612409601408e-7,9.704459174192359e-7,9.58030593878331e-7,9.458402983494787e-7,9.338449477363838e-7,   9.218495971232891e-7,9.099445015299672e-7,8.980552816649877e-7,8.865691437366877e-7,8.771226750527884e-7,8.67676206368889e-7,8.604063361970309e-7,8.554981404876014e-7,8.505899447781719e-7,   8.45965467421851e-7,8.414145063002377e-7,8.36870668157861e-7,8.323851723641769e-7,8.278996765704929e-7,8.229777591882111e-7,8.174615167119161e-7,8.119452742356211e-7,8.068986428948552e-7,   8.020187287402743e-7,7.97170820728625e-7,7.928993966947452e-7,7.886279726608653e-7,7.844750742472907e-7,7.805264348198531e-7,7.765777953924155e-7,7.723140844312032e-7,7.679033245258074e-7,   7.634925646204116e-7,7.583921238025582e-7,7.53289186433596e-7,7.49616095778962e-7,7.491108669201028e-7,7.486056380612436e-7,7.484164748860745e-7,7.484164748860745e-7,7.484164748860745e-7,   7.484164748860745e-7,7.484164748860745e-7,7.484164748860745e-7,7.484164748860745e-7,7.484164748860745e-7,7.484164748860745e-7,7.484164748860745e-7,7.484164748860745e-7,   7.484164748860745e-7,7.484164748860745e-7,7.484164748860745e-7,7.484164748860745e-7,7.484164748860745e-7,7.484164748860745e-7,7.484164748860745e-7,7.484164748860745e-7,   7.484164748860745e-7,7.484164748860745e-7,7.486326590018433e-7,7.501510949021517e-7,7.5166953080246e-7,7.5507081708941e-7,7.607095383894657e-7,7.663482596895214e-7,7.718708673676847e-7,   7.773590983725809e-7,7.828564088070003e-7,7.884500101962749e-7,7.940436115855495e-7,7.99790138024983e-7,8.05765322474927e-7,8.117405069248709e-7,8.178017044525785e-7,8.238971259124898e-7,   8.29996222249347e-7,8.362143002158433e-7,8.424323781823396e-7,8.487939530180891e-7,8.554284966675317e-7,8.620630403169741e-7,8.705623709518909e-7,8.800259643610491e-7,8.894895577702071e-7,   9.01033569021355e-7,9.126336532681836e-7,9.243508473482959e-7,9.363570657300174e-7,9.483632841117388e-7,9.60222145954891e-7,9.719839719334927e-7,9.837457979120944e-7,9.964954041190646e-7,   1.0093337728480576e-6,1.0222764886696772e-6,1.0355643454942937e-6,1.04885220231891e-6,1.0619921752723847e-6,1.0750095513212539e-6,1.0880269273701229e-6,1.1017430642494136e-6,   1.1155716865453261e-6,1.129621624355504e-6,1.1447081469858667e-6,1.159794669616229e-6};
      
	int which = floor(ms/0.001)-1;        
	return sqrt(peakBounds[which]);
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

double boundUtau(double ms){

	double result= 5.684603155097915*exp(-224.36106435213696*ms) +0.1538995418716346*exp(-62.349043302533644*ms)+ 0.0005635403465693567*exp(-6.096293940788484*ms);

	if(ms<0.009){ result = 1; }
	return sqrt(result);
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

double nuisFuncQE(const std::vector<double> &x, std::vector<double> &grad, void *my_data){
        
        nuisStruct *d = reinterpret_cast<nuisStruct*>(my_data);
        std::vector<double > qeGram = d->egram;
        double sigma_zeta = d->Sigma_Zeta;


        double zeta_b = x[0];
                                 //{204, 280, 214, 99, 83, 59, 51, 33, 37, 23, 19, 21, 12, 16, 4, 9, 4, 7, 3}
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

