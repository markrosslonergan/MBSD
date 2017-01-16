#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <string>
#include <sstream>
#include "flux.h"
#include "stat.h"
#include "flight.h"
#include "transit.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <ctime>
#include <globes/globes.h>







static void show_usage(std::string name)
{
    std::cerr << "\n Usage: " << name << " <options>\n"
	      << "Minimizer for finding ''local'' minima in NON unitary PMNS given experimental inputs\n"
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
              << "\t-f,--flux \t \tfun on flux generation mode\n"
              << "\t-c,--contour \t \tRun to generate .ddat's for pythoncontour plotting\n"
              << std::endl;
}


double f (double x, void * params) {
  double alpha = *(double *) params;
  double f = log(alpha*x) / sqrt(x);
  return f;
}



using namespace std;

int main(int argc, char* argv[]){
double extern MsMin;
double extern MsMax;
double extern MsStep;
double extern Umin;
double extern Umax;
double extern Ustep;
double extern GdMin; //These numbers made up for now
double extern GdMax;
double extern GdStep;
double extern eBinMin; //These numbers made up for now
double extern eBinMax;
double extern eBinStep;

vector<int > extern visdat;
vector<double > extern visback;
vector<double > extern visblank;


    bool testrun = false;
    bool fluxrun = false;
    bool contourrun = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_usage(argv[0]);
            return 0;
        } else if ((arg == "-t") || (arg == "--test")) {
	    testrun=true;
             
        } else if ((arg == "-f") || (arg == "--flux")) {
	    fluxrun=true;
             
        } else if ((arg == "-c") || (arg == "--contour")) {
	    contourrun=true;
             
        } 
    }

//#####################################################################################
//#####################################################################################
//################ Directory and File folder creation, and streams ####################
//#####################################################################################
//#####################################################################################

		time_t now = time(0);
		tm *ltm = localtime(&now); // time structure to read out mon:day:hour:min
	//cout << "Year: "<< 1900 + ltm->tm_year << endl;

		 std::vector<string > num2mon = {"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
 		 double nowmonth =  1 + ltm->tm_mon;
		 double nowday =  ltm->tm_mday;
   		 double nowhour =  ltm->tm_hour;
		 double nowmin =  ltm->tm_min;

		std::stringstream date_hour_min1; date_hour_min1.str( std::string() );	date_hour_min1.clear();
		std::stringstream date_hour_min2; date_hour_min2.str( std::string() );	date_hour_min2.clear();
		std::stringstream date_hour_min3; date_hour_min3.str( std::string() );	date_hour_min3.clear();
		std::stringstream cpstream; cpstream.str( std::string() );	cpstream.clear();

	 	date_hour_min1 << "mkdir -p  data/"<<num2mon[nowmonth]<<nowday<<"_"<<nowhour<<"_"<<nowmin<<"/flux";
		system(date_hour_min1.str().c_str());
	 	date_hour_min2 << "mkdir -p  data/"<<num2mon[nowmonth]<<nowday<<"_"<<nowhour<<"_"<<nowmin<<"/flight";
		system(date_hour_min2.str().c_str());
	 	date_hour_min3 << "mkdir -p  data/"<<num2mon[nowmonth]<<nowday<<"_"<<nowhour<<"_"<<nowmin<<"/transit";
		system(date_hour_min3.str().c_str());
	
	cpstream<< "cp globals.cxx  data/"<<num2mon[nowmonth]<<nowday<<"_"<<nowhour<<"_"<<nowmin<<"/";
	system(cpstream.str().c_str()); //Copy accross the config file. Has to be a faster way than this eh!


// Opening up all necessary streams
	std::ofstream teststream;  
	teststream.open("test.log");

//#####################################################################################
//#####################################################################################
//################ Next bit ####################
//#####################################################################################
//#####################################################################################

//std::vector<double> res1 =  kinconregions(0.9999,1.0,1,0.001,MKao) ;
//std::vector<double> res2 =  removeNAN(kinconregions(0.9999,1.0,1,0.001,MKao) );

//cout<<res1[0]<<"   "<<res1[1]<<"  "<<res1[2]<<"   "<<res1[3]<<endl;
//cout<<res2[0]<<"   "<<res2[1]<<"  "<<res2[2]<<"   "<<res2[3]<<endl;




//int gsl_integration_qag (const gsl_function * f, double a, double b, double epsabs, double epsrel, size_t limit, int key, gsl_integration_workspace * workspace, double * result, double * abserr) 


//#################################################################################
//                    Read a file in for fluxfile vector 
//#################################################################################
std::vector< std::vector<double > > fluxfile;

if (true){
	// fluxfile = readmein("table.dat");
	 fluxfile = readmein("table_long.dat");
} else {

	fluxfile.reserve(387);
	double mm = 0;

	for(int im = 0; im < 387; im=im+1) {
		mm = im*MsStep+MsStep;
		//	cout<<"im: "<<im<<"  mm: "<<mm<<endl;
		for (double ee=eBinMin+eBinStep; ee<eBinMax; ee=ee+eBinStep) {
			double fluxsum = lint(mm,ee,MKao)+lint(mm,ee,MP);
			cout<<fluxsum<<"  ";
			fluxfile[im].push_back(fluxsum);
		}
	cout<<endl;
	}

}




//#################################################################################
//                    Kaon Flux Practice Area 
//#################################################################################
//for(double ee=0.01; ee<=10; ee=ee+0.05){
//cout<<ee<<"  "<<Kflux_Asaka(ee)<<endl;
//}


//#################################################################################
//         decayflight returns a visbined result for a given M,U,Gdec
//#################################################################################
//dispvec(decayflight(0.018,pow(10,-2.65),1e-15,fluxfile)); // Works tick tick
//cout<<endl<<ebin[2]<<"  "<<ebin[3]<<endl;
//std::cout<<"Probability to decay is "<<probdec(1e-15,550,0.018,2)<<std::endl;
	

//#################################################################################
//         decayTransit Stuff
//#################################################################################
// double diffdecay_transit (double Ms, double Es, double Ev, double U
// double decaytransitwithE (double Ms, double Ev, double U, double Gdec, const std::vector<std::vector<double > > &fluxptr)  Seems to work fine
// std::vector<double > transitspectrum (double Ms, double U, double Gdec, const std::vector<std::vector<double > > &fluxptr) {
	 glbInit(argv[0]); 
//To normalise, we send the neutrino muon flux (estimated as fluxfile[0]) though globes and normalise to obtain 191 866  CCQE events.
//cout<<"Globes Normalization, m=0 flux, numu CCQE"<<endl;
//dispvec(globes_norming(fluxfile[0]));
//cout<<"Globes Spectrum, m=33 flux, nue CCQE"<<endl;
//globes_spectrum(fluxfile[0]);
//cout<<"Globes Spectrum, m=33 Decay in transit flux, nue CCQE"<<endl;
//globes_spectrum(transitspectrum(0.001,0.01,1e-15,fluxfile));


//std::cout<<"Probability to decay is "<<probdec(1e-15,550,0.03,2)<<std::endl;
//for(double ee=0.001;ee<=0.05; ee=ee+0.001) {
//cout<<"Es: "<<ee<<"  transit :"<<transitspectrum_withE (0.033, ee, 1, 1e-15, fluxfile)<<"   Probability to decay is "<<probdec(1e-15,550,0.033,ee)<<endl;
//}


//dispvec(decayflight(0.001,0.1,1e-15,fluxfile));
//cout<<"---------------------------------------------"<<probdec(1e-15,540,0.001,2)<<endl;
//dispvec(decaytransit(0.001,0.1,1e-15,fluxfile));

//decayflightwithE (0.3, 0.01, 0.1, 1e-15, fluxfile);
//decayflightwithE (0.3, 0.1, 0.1, 1e-15, fluxfile);
//decayflightwithE (0.3, 1, 0.1, 1e-15, fluxfile);
//decayflightwithE (0.3, 5, 0.1, 1e-15, fluxfile);
//dispvec(decayflight(0.01,0.1,1e-15,fluxfile));

//dispvec(fluxfile[32]);
//dispvec(transitspectrum(0.033,1,1e-15,fluxfile));

//cout<<transitspectrum_withE(0.04,2,1,1e-15,fluxfile)<<endl;



//dispvec(decayflight(0.250,pow(10,-5.0),1e-20,fluxfile));
//cout<<"---------------------------------------------"<<endl;
//decaytransit(0.101,0.01,1e-15,fluxfile);
//cout<<"---------------------------------------------"<<endl;
//dispvec(decaytransit(0.001,0.1,1e-15, fluxfile));

dispvec(decaytransit(0.101,0.01,1e-15, fluxfile));




//for (double m = 0.001; m<0.38;m=m+0.002){
//std::cout<<"M: "<<m<<" /0.38"<<std::endl;
if(false){
double m=0.08;
	for (double u=-8;u<=-0.05;u=u+0.05){
std::cout<<"U: "<<u<<std::endl;
		for (double gam = -30;gam <=-14;gam=gam+0.1){

			std::vector<double > kflight  = decayflight(m,pow(10,u),pow(10,gam),fluxfile);
			std::vector<double > ktransit = decaytransit(m,pow(10,u),pow(10,gam),fluxfile);
			teststream<<" "<<m<<" "<<u<<" "<<pow(10,gam)<<" "<<ifnan(vectorsum(kflight))<<"  "<<ifnan(vectorsum(ktransit))<<"  ";
			for(int i=0; i<kflight.size()-1;i++) {
				teststream<<ifnan(kflight[i])<<" ";
			}
			teststream<<ifnan(kflight[kflight.size()-1])<<"  ";
			for(int i=0; i<ktransit.size()-1;i++) {
				teststream<<ifnan(ktransit[i])<<" ";
			}
			teststream<<ifnan(ktransit[ktransit.size()-1])<<std::endl;
		}
	}
}//end if
//}

//#################################################################################
//         Log Liklihood stuff
//#################################################################################


//cout<<LogLikli(visdat, visback, decayflight(0.028,pow(10,-3.8),1e-25,fluxfile))<<endl;
//double likli_no_sig = LogLikli(visdat, visback, visblank);
//for (double ui=Umin; ui > Umax ;ui=ui+Ustep){
//	for (double mi=MsMin; mi<MsMax; mi=mi+MsStep){
//	cout<<ui<<"  "<<mi<<"  "<<LogLikli(visdat, visback, decayflight(mi,pow(10,ui),1e-25,fluxfile))<<endl;;
//	}
//}



teststream.close();
return 0;
} 
