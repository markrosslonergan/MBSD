#ifndef   BASEFNC_H
#define   BASEFNC_H

// Parameters as defined in the globals.cxx file.
extern double Mdau;
extern double Mpar;
extern double MKao;
extern double Gf;
extern double Piontotal;
extern double FPI;
extern double TWein;
extern double PI;
extern double VUD;
extern double MP;

extern double DECAYTUNNEL;

extern double CosIntBoundLow; 		//eventually be a function of detector.
extern double CosIntBoundHigh;

extern double TFUDGE ;
extern double PFUDGE ;
extern double SWNORM;

extern std::vector<double > ebin;
extern double BASELINE;
extern double DECvol;
extern double DECpot;
extern double DECeff;
extern double Gf;
extern double FPI;
extern double FKA;
extern double VUS;
extern double TWein;
extern double PI;
extern double VUD;
extern double Me;
extern double MMu;
extern double MP;


double gamma (double Ms, double Es){ return Es/Ms; }
double beta (double Ms,double Es) { return sqrt(1-pow(gamma(Ms,Es),-2) );}

double Gvvv (double M, double U)   {      return Gf*Gf*U*U*pow(M,5)/(96.0*PI*PI*PI);    }
double Gsg (double M, double U)    {      return Gvvv(M,U)*27.0/(32*137*PI);    }
double Gvee (double M, double U)   {      return Gvvv(M,U)*(0.25-pow(sin(TWein),2)+2*pow(sin(TWein),4) );    }
double Gpil (double M, double U)   {      return Gf*Gf*FPI*FPI*VUD*VUD*U*U*M*M*M/(16.0*PI);    }
double Gpiv (double M, double U)   {      return Gf*Gf*FPI*FPI*VUD*VUD*U*U*M*M*M/(64.0*PI);    }

double decaytotal (double Ms, double U, double Gdec) {
	double result = Gdec+Gvvv(Ms,U)+Gsg(Ms,U);
	if (Ms > 2*Me) {      result = result + Gvee(Ms,U); 	}
	if (Ms > MP)      {   result = result + Gpiv(Ms,U);      }
	if (Ms > (MP+Me) ) {  result = result + Gpil(Ms,U);      }
	return result;
} 

double M2G (double m) {
return (100*1e5*1e9*m)/1.973;
}

void printhere() {
	// printhere
	// just helps me print here
	std::cout<<"HERE!"<<std::endl;
}

void dispvec (std::vector<double > inputlist) {
	// prints a vector

	for(int i=0; i<inputlist.size();i++){
		std::cout<<inputlist[i]<<std::endl;
	}
}

void dispvecHor (std::vector<double > inputlist) {
	// prints a vector horizontally
	for(int i=0; i<inputlist.size();i++){
		std::cout<<inputlist[i]<<" ";
	}
	std::cout<<std::endl;
}

std::vector<std::vector<double > > readmein (std::string filename) {
	// readmein
	// func takes filename, and reads in a Ms times Ebin matrix into a vec<vec<double>>

	std::ifstream inputFile(filename); //inputfile 

	 std::string line;
	 std::vector< std::vector<double> > readflux;
	   while ( std::getline( inputFile, line ) ) {
	      std::istringstream is( line );
	      readflux.push_back(     std::vector<double>( std::istream_iterator<double>(is)   , std::istream_iterator<double>() ) );
	   }
return readflux;
}

double vectorsum (std::vector< double> input) {
	double ans = 0;
		for(int i=0;i<input.size();i++){
			ans=ans+input[i];
		}
	return ans;
}



// whichbin finds which bin in ebins, that a certain gamma energy is above.
int whichbin (double Ms, double Egam ){
	double tmp = Egam*(1+ Ms*Ms/(4.0*Egam*Egam));
	int result=0;
	
	for (int ie = 0; ie < ebin.size(); ie++) {
		if (ebin[ie] > Egam) {result = ie-1; break;}
	}
return result;	
}

int whichbinM (double M){
	int resulter=0;
	
	for (int ie = 0; ie < ebin.size(); ie++) {
		if (ebin[ie] > M) {resulter = ie-1; break;}
	}
return resulter;	
}

double probdec(double Gdec, double len, double Ms, double Es){

return 1-std::exp(-Gdec*M2G(len)/(beta(Ms,Es)*gamma(Ms,Es)));
}

std::vector<double> removeNAN (std::vector<double >inputlist) {
	// removes nan's from a list and sorts remaing reals
	// uses fact NAN != NAN 	
	// Used in sorting the kinematic boundaries of 2-body decays
	// Remove negative numbers too.
	
	std::vector<double> result;


	for(int i=0; i<inputlist.size();i++){
		if( inputlist[i]==inputlist[i]&&inputlist[i]>0.0){
				result.push_back(inputlist[i]);
		}
	}
	std::sort(result.begin(), result.end()) ; //Doesnt return, sorts at a pointer level
return result;
}

double ifnan(double inp){
if(inp!=inp){return 0;}
return inp;
}







#endif
