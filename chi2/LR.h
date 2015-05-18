#ifndef LR_H_
#define LR_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#define NUMEVENTS 200000
#define NUM_EVENT_OBS 7 //Number of obervables stored for each event. E_sum/Th_sum/AngSep/E_sterile/E_high/E_low/Th_high

#define EBINS 19
#define COSBINS 10

typedef struct CL_input { double mS; double mZprime; double eCut; double thCut; double eFloor; double eRatio;double SysOn; double Sigma_Zeta; } CL_input;

double getTotalNumEvents(CL_input in);

int EtoBin(double E);
double BintoCentralE(int bin);
int CostoBin(double C);
double BintoCentralCos(int b);

double decayProb(CL_input input, double chiU, double Es);
double histogrammer(CL_input in, double chiU, double cutEff, const double events[][NUM_EVENT_OBS], double eGram[], double cosGram[]);
int printEGram(double eGram[]);
int printCosGram(double cosGram[]);

double boundU( double ms);
double boundChi(double ms);
double boundChiU(double ms, double mz);
double boundUtau(double ms);

typedef struct {std::vector<double >  egram; double Sigma_Zeta; } nuisStruct;
double nuisFuncE(const std::vector<double> &x, std::vector<double> &grad, void *my_data);
double nuisFuncA(const std::vector<double> &x, std::vector<double> &grad, void *my_data);
double nuisMarginalize(std::vector<double > * bf_zeta_b, double * chi, std::vector<double > * eVGram,int whi, double SIGMAZETA);


#endif

