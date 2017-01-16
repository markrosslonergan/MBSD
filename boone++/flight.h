// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005  

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/
#include <vector>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <algorithm>
#include "basefunc.h"
#ifndef  BASEFNC_H  
#include "basefunc.h"
#endif


// ############################################################################################
// ##########################   decay functions				#######################
// ############################################################################################

// todo: make a list for this to read from, do not call lint on the fly each time. almost done.
double decayflightwithE (double Ms, double Egam, double U, double Gdec, const std::vector<std::vector<double > > &fluxptr) {

	double result = 0.0;
	double binf = 0.025;

	for (int ie = whichbin (Ms, Egam) ; ie < ebin.size(); ie++) {
		double Es = ebin[ie]+binf;
		double G = gamma(Ms,Es);
		double B = beta(Ms,Es);
		//std::cout<<"Es: "<<Es<<" egam: "<<Egam<<" which: "<<whichbin (Ms, Egam) <<" exp: "<<exp(  -BASELINE*decaytotal(Ms,U,Gdec)/(G*B) )<<" flux: "<<fluxptr[(int)floor(Ms/0.001)-1][ie]<< " otherbit: "<<0.05*(U*U*Ms)/(pow(ebin[ie]+binf,2)*beta(Ms,ebin[ie]+binf))<<" addition: "<<0.05*(U*U*Ms)/(pow(ebin[ie]+binf,2)*beta(Ms,ebin[ie]+binf))*fluxptr[(int)floor(Ms/0.001)-1][ie]*exp(  -BASELINE*decaytotal(Ms,U,Gdec)/(G*B) )<<"  result: "<<result<<std::endl;
		result = result + DECvol*DECpot*DECeff*Gdec*0.05*(U*U*Ms)/(pow(ebin[ie]+binf,2)*beta(Ms,ebin[ie]+binf))*fluxptr[(int)floor(Ms/0.001)-1][ie]*exp(  -BASELINE*decaytotal(Ms,U,Gdec)/(G*B) );
	} 

return result; 
}

std::vector<double > decayflight (double Ms, double U, double Gdec, const std::vector<std::vector<double > > &fluxptr) {
	std::vector<double > result;
	for (int ee=2; ee<=38; ee=ee+2) {
		result.push_back(decayflightwithE(Ms,ebin[ee],U,Gdec, fluxptr) + decayflightwithE(Ms,ebin[ee+1],U,Gdec, fluxptr));
	}
return result;
}
