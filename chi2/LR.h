#ifndef LR_H_
#define LR_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include "Minuit2/FCNBase.h"

#define NUMEVENTS 200000
#define EBINS 19
#define COSBINS 10

typedef struct CL_input { double mS; double mZprime; double eCut; double thCut; } CL_input;

double getTotalNumEvents(CL_input in);

int EtoBin(double E);
double BintoCentralE(int bin);
int CostoBin(double C);
double BintoCentralCos(int b);

double decayProb(CL_input input, double chiU, double Es);
double histogrammer(CL_input in, double chiU, double cutEff, const double events[][4], double eGram[], double cosGram[]);
int printEGram(double eGram[]);
int printCosGram(double cosGram[]);


namespace ROOT { namespace Minuit2 {

//class log_likelihood : public ROOT::Minuit2::FCNBase {
class E_log_likelihood : public FCNBase {
	
public:

	E_log_likelihood(double input_evs[NUMEVENTS][4]) 
	{
		int i,j;
		for(i=0;i<NUMEVENTS;i++)
		{
			for(j=0;j<4;j++)
			{
				events[i][j] = input_evs[i][j];
			}
		}
		
	}

	~E_log_likelihood() {}

	virtual double Up() const { return 0.5; }
	double operator()(const std::vector<double>& p) const;

private:

	double events[NUMEVENTS][4];
};

}} // closing namespaces
#endif

