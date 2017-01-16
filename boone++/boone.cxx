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
#include <boost/program_options.hpp>
#include <vector>
#include <iomanip>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <ctime>
#include <globes/globes.h>
#include "zprimecrosssection.h"


using namespace std;
namespace po = boost::program_options;
//test for 

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

vector<double > batchlist;

    bool testrun = false;
    bool fluxrun = false;
    bool r1run   = false;
    bool batchrun = false;
int selector = 0;

//#####################################################################################
//				Parameter Options Readining In!
//#####################################################################################

// Declare the supported options.
po::options_description desc("Usage Options for boone++");
desc.add_options()
    ("help,h", "produce help message")
    ("flux,f", "Run old decay in flight flux calculations")
    ("sterile,s", po::value<double>(),"Set Sterile Mass (GeV). Takes a double.")
    ("zprime,z",  po::value<double>(),"Set Z Prime Mass (GeV). Takes a double")
    ("runone,r",  po::value<int>(), "Run a single job. Takes a Int.")
    ("batch,b",   po::value<std::vector<double > >()->multitoken(), "Runs a batch montecarlo. Takes 7 doubles. Produces HIST and MC")
    ("polarity,p", po::value<int>(),"set negative or positive polarity")
;

po::variables_map vm;
po::store(po::parse_command_line(argc, argv, desc), vm);
po::notify(vm);    

if (vm.count("help")) {
    std::cerr << desc << "\n";
    return 1;
}

if (vm.count("flux")) {
   fluxrun = true;
}

if (vm.count("sterile")) {
    std::cout << "Sterile Mass (GeV) : " 
 << vm["sterile"].as<double>() << ".\n";
} else {
   //std::cout << "Sterile Mass Not Set.\n";
}

if (vm.count("zprime")) {
    std::cout << "Z Prime Mass (GeV) : " << vm["zprime"].as<double>() << ".\n";	
} else {
   //std::cout << "Z Prime Mass Not Set.\n";
}

if (vm.count("batch")) {
    std::cout << "Reading in Batch parameters"<<std::endl;
    batchlist = vm["batch"].as<std::vector< double>>() ;
    batchrun = true;

	if (batchlist.size()!=7){
		std::cerr<<"ERROR: PO::BATCH. Batchparameters do not number 6"<<std::endl;
	} 
} 

if (vm.count("runone")) {
	std::cout<<" RunOne Starting"<<std::endl;
    	r1run = true;
}

if(vm.count("polarity")){
	selector = vm["polarity"].as<int>();
	if (selector < 0 || selector > 3){
		std::cerr<<"ERROR PO:POLARITY must select 0,1,2,3 polarity corresponding to posnumu,posnumubar,negnumu, negnumubar"<<std::endl;
	}
}

//#####################################################################################
//#####################################################################################
//################ Directory and File folder creation, and streams ####################
//#####################################################################################
//#####################################################################################
if(false){
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
}

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

if (!fluxrun){
	// fluxfile = readmein("table.dat");
	 fluxfile = readmein("table_long.dat");
} else {

	fluxfile.reserve(387);
	double mm = 0;

	for(int im = 0; im < 387; im=im+1) {
		mm = im*MsStep+MsStep;
		//	cout<<"im: "<<im<<"  mm: "<<mm<<endl;
		for (double ee=eBinMin+eBinStep; ee<eBinMax; ee=ee+eBinStep) {  //WARNINGS CHANGE eBinMax to 10!!
			double fluxsum = lint(mm,ee,MKao)+lint(mm,ee,MP);
			cout<<fluxsum<<"  ";
			fluxfile[im].push_back(fluxsum);
		}
	cout<<endl;
	}

}

//dispvec(fluxfile[10]);



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

//dispvec(decaytransit(0.101,0.01,1e-15, fluxfile));


//dispvec(decayflight(0.250,pow(10,-5.0),1e-18,fluxfile));



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

//cout<<LogLikli(visdat, visback, visblank)<<endl;
//for (double ui=Umin; ui > Umax ;ui=ui+Ustep){
//	for (double mi=MsMin; mi<MsMax; mi=mi+MsStep){
//	cout<<ui<<"  "<<mi<<"  "<<LogLikli(visdat, visback, decayflight(mi,pow(10,ui),1e-25,fluxfile))<<endl;;
//	}
//}

//#################################################################################
//         		NEW z prime stuff Initilising the Random Number Generator
//#################################################################################

	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();
	
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

//#################################################################################
//         -r --runone	  r1: 1 Run of the data for aforred param
//#################################################################################
if(r1run){
	std::cout<<"Starting Single Run"<<std::endl;
	MCmassrun(vm["zprime"].as<double>(),vm["sterile"].as<double>(), vm["runone"].as<int>(),r, selector);
}
//#################################################################################
//         -b --batch	  batch Mode: Batch mode
//#################################################################################
if(batchrun){
	std::cout<<"BATCHMODE: Begining Batch"<<std::endl;
	double zstart = batchlist[0];
	double zstep  = batchlist[2];

	for( double is =batchlist[3]; is <= batchlist[4]; is+=batchlist[5]){
			zstart=batchlist[0];
		if(zstart <= is){
			zstart = ceil( static_cast<double>(is/zstep))*zstep;
			//std::cout<<ceil(0.11/0.02)*0.02<<" "<<ztep<<"  "<<zstart<<"  "<<zstart*batchlist[1]<<std::endl;
		}			
		for( double iz =zstart; iz <=batchlist[1]; iz+=zstep){
			std::cout<<"BatchRun~ Ms: "<<is<<" Mz: "<<iz<<std::endl;
			MCmassrun(iz,is,batchlist[6],r,selector) ;
		}
	}
}


//#################################################################################


//for(double ev=1.9;ev<=2;ev=ev+1e-6){
//	std::cout<<"ev: "<<ev<<" Flux: "<<fluxTempNormed(ev)<<std::endl;//returns negative value!!!!
//	std::cout<<"ev: "<<ev<<" cs: "<<totalCS (ev, 0.1, 0.001, 1.0,1.0)<<std::endl;//returns negative value!!!!
//}
//for(double izz=0.01; izz<=100; izz=izz+0.01){
//std::cout<<izz<<" "<<totalevents(1,1,izz,0.2)<<std::endl;
//}

//for(double iss=0.001; iss<=0.6; iss=iss+0.001){
//	for(double izz=0.1; izz<=1; izz=izz+0.01){
//		std::cout<<iss<<" "<<izz<<" "<<totalevents(1,1,izz,iss)<<std::endl;
//	}
//}

//double ppp= CSffactor(2, 1, 0.2, 0.1, 1, 0.1);
//double Tppp = totalCS(2, 2, 0.1, 0.2, 1);
//double teve = totalevents(1e-2,1e-2, 0.4, 0.2);


//double is = 0.021;
//for( double iz = 0.1; iz <= 1; iz=iz+0.01){
//	for( double is = 0.001; is <= min(0.2,iz);is=is+0.005){
//		std::cout<<"Ms/0.2:"<<is<<"Mz/1:"<<iz<<std::endl;
//		MCmassrun(iz,is,500000,r) ;
//	}
//}


//MCmassrun(0.1,0.02,2000,r);

//for( double iz = 0.5; iz <= 2.0; iz=iz+0.1){
//std::cout<<iz<<std::endl;
//MCmassrun(0.1,0.05,20000,r);

//for(double Evi =LABevMin(0.05,0.0,Mnuc); Evi<=3;Evi=Evi+0.01){
//double Evi=2;
//	for(double it = mandelTmax(0.05,Evi,0.0,1.0); it <= mandelTmin(0.05,Evi,0.0,1.0); it=it + 0.001){ 
//		std::cout<<Evi<<"  "<<it<<" "<<CSffactor(-it, Evi,  0.1, 0.05,  1.0,  1.0)<<" "<<CSffactor2(-it, Evi,  0.1, 0.05,  1.0,  1.0)<<std::endl;
//	}
//}

//std::cout<<totalCS(0.2, 0.3, 0.09, 0.0031, 0.01)<<"  "<<totalCS(0.2,91,0.00,1.0, 1.0)<<std::endl;
//std::cout<<totalCS(1, 0.3, 0.09, 0.0031, 0.01)<<"  "<<totalCS(1,91,0.00,1.0, 1.0)<<std::endl;
//std::cout<<totalCS(2, 0.3, 0.09, 0.0031, 0.01)<<"  "<<totalCS(2,91,0.00,1.0, 1.0)<<std::endl;

//std::cout<<totalCS(0.2, 0.6, 0.03, 0.01, 0.031)<<"  "<<totalCS(0.2,91,0.00,1.0, 1.0)<<std::endl;
//std::cout<<totalCS(1, 0.6, 0.03, 0.01, 0.031)<<"  "<<totalCS(1,91,0.00,1.0, 1.0)<<std::endl;
//std::cout<<totalCS(2, 0.6, 0.03, 0.01, 0.031)<<"  "<<totalCS(2,91,0.00,1.0, 1.0)<<std::endl;


//double pref=pow(10e9*5.068*10e4,-2)/(10e-38);
//double eee=10e-2;
//for(double le=-2;le<=2;le=le+0.02){
//eee=pow(10,le);
//std::cout<<eee<<"  "<<totalCS(eee, 0.3,0.09, 1,pow(0.000632659,2))*pref/eee<<"  "<<totalCS(eee,0.6, 0.03, 1,pow(1.84469e-06,2))*pref/eee<<"  "<<totalCS(eee,0.45, 0.05, 1,pow/( 0.00153338,2))*pref/eee<<"  "<<totalCS(eee, 91, 0.001, 1,1)*pref/eee<<std::endl;
//}

//double eee=pow(10,-2);
//for(double mz=-2;mz<=2;mz=mz+0.02){
//eee=pow(10,mz);
//std::cout<<eee<<"  "<<totalCS(1, eee,0.01, 1,1)/totalCS(1,91,0.01,1,1)<<std::endl;
//}


//std::vector<double > aop = MaxTotalCs(0.1,0.05,1,1);
//std::vector<double > aop2 = MaxDiffCs(0.1,0.05,1,1);
//dispvecHor(aop);
//dispvecHor(aop2);
//}
//################################# int.dat for KINK anaysis ######################
/*for(double eee=0.001;eee<=5;eee=eee+0.02){
std::cout<<eee<<"  "<<totalCS(eee, 0.3,0.05, 1.0, 1.0)<<"  "<<totalCS(eee,0.3, 0.001, 1.0, 1.0)<<"  "<<totalCS(eee, 0.3, 0.1, 1.0, 1.0)<<" "<<totalCS(eee, 0.8, 0.05, 1.0, 1.0)<<"  "<<totalCS(eee, 0.8, 0.001, 1.0, 1.0)<<"  "<<totalCS(eee,0.8,0.1,1.0, 1.0)<<std::endl;
}
*/
//################################# What Am I Sampling ######################
//double Evin = 6.0;
//for(double it = mandelTmax(0.05,Evin,0.0,1.0); it<=mandelTmin(0.05,Evin,0.0,1.0); it=it + 0.001){ 
//std::cout<<it<<" "<<CSffactor(-it, Evin,  91, 0.05,  1.0,  1.0)<<std::endl;
//}


//################################# out a given CS total ######################
//double smm=totalCS(2,91, 0.05, 1, 1);
//for(double mmz=0.1; mmz <= 1; mmz=mmz+0.01){
//std::cout<<mmz<<" "<<totalCS(2, mmz, 0.05, 1e-4, 1)/smm<<std::endl;
//}
//for( double ei= 0.01; ei<=3;ei=ei+0.01){
//std::cout<<ei<<" "<<totalCS(ei, 91, 0.05, 1, 1)<<std::endl;
//}
//std::vector<double > aop = MaxTotalCs(2,0.4,1,1);
//dispvecHor(aop);

//################################  Total events for fixed _ (ms or mz)  ##################

//std::cout<<CSProton(-0.2,2,0.5,0.08,1,1)<<"\t"<<CSNeutron(-0.2,2,0.5,0.08,1,1)<<std::endl;
//std::cout<<totalCSProton(2,0.5,0.08,1,1)<<" "<<totalCSNeutron(2,0.5,0.08,1,1)<<std::endl;
//std::cout<<"old "<<totalevents(1,1,0.5,0.08)<<" new "<<totalevents2(1,1,0.5,0.08)<<std::endl;
//double mz =0.5;
//	std::cout<<mz<<"\t"<<totalevents(1,1, mz, 0.08)<<"\t"<<totalevents2(1,1,mz,0.08)<<"\t"<<totalevents(1,1,mz,0.08)/totalevents2(1,1,mz,0.08)<<std::endl;
//for( double mz = 0.1; mz<=1; mz=mz+0.01){
//	std::cout<<mz<<"\t"<<totalevents(1,1, mz, 0.08)<<"\t"<<totalevents2(1,1,mz,0.08)<<"\t"<<totalevents(1,1,mz,0.08)/totalevents2(1,1,mz,0.08)<<std::endl;
//}
//std::cout<<totalevents(1,1, 91, 0.0000001)<<std::endl;

//MCmassrun(0.7,0.05,20000,r);
//double ee=2; double mus=0.08; double mtz=0.5;
//for(double ti=mandelTmax(mus,ee,0.0,1); ti<=mandelTmin(mus,ee,0.0,1.0); ti=ti+0.0001){
//	std::cout<<ti<<" "<<CSProton(ti,ee,mtz,mus,1,1)<<" "<<CSNeutron(ti,ee,mtz,mus,1,1)<<std::endl;
//}
//std::cout<<"First test of CCextion: "<<ppp<<std::endl;
//std::cout<<"First test of TotalCS: "<<Tppp<<std::endl;
//std::cout<<"First test of TeventsS: "<<teve<<std::endl;
//for(int p = 0; p<=200000;p++){
//	dispvecHor(MCkin (0.1,0.001, r));
//}

//for(double ee=0.01; ee<10; ee=ee+0.1){
//
//std::cout<<ee<<" "<<fluxTemp(ee)<<std::endl;
//}


gsl_rng_free(r);
teststream.close();
return 0;
} 
