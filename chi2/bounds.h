#ifndef BOUNDS_H_
#define BOUNDS_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

// PLACE ALL BOUNDS ON U and chi not U^2 and Chi^2  
// All bounds on U valid 0.001 to 0.3
// All bounds on Chi valid 0.1 to 1

double GammaNeeded2DecayBefore(double mn, double en, double L, double L0, double assumedU);

double Gvee(double U, double ms);
double M2G (double m) ;

double boundBASEu(double ms);
double boundBASEzp(double mzp);
double boundPS191u(double ms);
double boundNOMADt(double ms);
double boundBABARzp( double mzp);

bool bound_is_legit_tau(double up, double ud, double chi, double ms, double mzp );
bool bound_is_legit_order1(double up, double chi, double ms, double mzp );

#endif

