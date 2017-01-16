/*globals.cxx*/

#include <vector>
#include <math.h>

double m2g (double m) {
return (100*1e5*1e9*m)/1.973;
}

std::vector<double > makerange (double start, double stop, double step) {
	std::vector<double > result;
	for(double i = start; i < stop; i=i+step) {
	   result.push_back(i);
	}
return result;
}
//##################  Generic physical constats #################
double MP =0.13957018;
double MMu = 0.1056583715;
double Me = 0.0005109989;
double Mdau = MMu;
double Mpar = MP;
double MKao = 0.493677;

double Mw = 80.385;
double Mz = 91.1876;
double Gf = 1.16637e-5;
double Piontotal = 2.5286208221283132e-17;
double FPI = 0.1307;
double FKA = 0.1561;
double VUS = 0.2252;
double VUD = 0.97425;
double TWein = acos(Mw/Mz);
double PI = 3.14159;

double theta12 = atan(sqrt(0.44)); //asin(sqrt(0.311));
double theta13 = asin(0.18); //asin(sqrt(0.0255));
double theta23 = M_PI/4;
double sdm = 7.7e-5;
double ldm_N = 2.4e-3;
double ldm_I = -2.4e-3;
double deltacp = M_PI/2;

//##################  Steps size for the entire run #################
double MsMin = 0.001;
double MsMax = 0.033;
double MsStep = 0.001;
double Umin = -0.2;
double Umax = -4.0;
double Ustep = -0.05;
double GdMin = -32.0; //These numbers made up for now
double GdMax = -1.0;
double GdStep = 0.5;
double eBinMin=0.00;
double eBinMax=9.95;  //WARNING CHANGED FROM 4.95 to 9.95
double eBinStep=0.05;
std::vector<double > ebin = {0.0,0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2., 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3., 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4., 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5, 4.55, 4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 4.9, 4.95, 5.};


//##################  Integration cosint bits, can be replaced with on-off axis #################
double CosIntBoundLow = 0.9999;
double CosIntBoundHigh = 1.0;

//##################### Detector and Experimental details #######################
double BASELINE = m2g(540.0);
double DECAYTUNNEL = m2g(50.0);
double DECvol = 4.818934563976785e22 ;
double DECpot = 6.46e20;
double DECeff= 0.9;

//##################### SW parameters and kaon flux parameters #######################
double TFUDGE = 0.1;
double PFUDGE = 0.3;
double SWNORM = 2.1325332782450387e-12;

//#################### MiniBooNE specific vectors of visdata ..etc.. #################
std::vector<double > visbin = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2,1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9};
std::vector<int > visdat = {204, 280, 214, 99, 83, 59, 51, 33, 37, 23, 19, 21, 12, 16, 4, 9, 4, 7, 3};
std::vector<double > visback = {151.5, 218.8, 155.6, 108.7, 72.5, 57.6, 45, 38.5, 31.4,22.2, 20.4, 17.2, 14.1, 10.2, 9.1 , 8.2 , 5.6 , 5.7 , 2.9};
std::vector<double > visdaterr = {14.3, 16.7, 14.6, 9.9 , 9.1 , 7.7, 7.1 , 5.7, 6.1,   4.8 , 4.4 , 4.6 , 3.5 , 4, 2, 3, 2 , 2.6, 1.7};
std::vector<double > visblank = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};



//##################### Statistics related variable #################################



//##################### Z Prime Related Variable #################################
double Mnuc = 1;

