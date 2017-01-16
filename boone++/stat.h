// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005  

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/
#include <vector>
#include <math.h>
#include <algorithm>
#include <gsl/gsl_sf_gamma.h>

// Need log here for computational reasons
long double LogPoissonProb (int n, double s, double b) {
	return       log(s+b)*n -(s+b) ;//- gsl_sf_lnfact(n);
}


double LogLikli (std::vector<int > data, std::vector<double > backg, std::vector<double > signal) {
	double result=0.0;
	if (data.size()!=signal.size() ) {std::cout<<"ERROR data vector and signal vector are not equal in size. LogLikli fucn stat.h"<<std::endl;}
	int datasize = signal.size();
	
	for (int i=0; i<datasize; i++){
		result = result + LogPoissonProb(data[i],signal[i],backg[i]);
	}

	return result;
}
