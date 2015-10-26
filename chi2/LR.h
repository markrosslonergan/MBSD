#ifndef LR_H_
#define LR_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#define NUMEVENTS 40000
#define NUM_EVENT_OBS 7 //Number of obervables stored for each event. E_sum/Th_sum/AngSep/E_sterile/E_high/E_low/Th_high

#define EBINS 19
#define COSBINS 10
#define QEBINS 11

#define MCHI 80 

typedef struct CL_input { double mS; double mZprime; double eCut; double thCut; double eFloor; double eRatio;double SysOn; double Sigma_Zeta; int which_var; int pos_selector;} CL_input;
typedef struct BF_RESULT { double E_bf; double A_bf; double QE_bf; double E_bf_Up; double E_bf_Ud; double E_bf_Chi; double A_bf_Up; double A_bf_Ud; double A_bf_Chi; double QE_bf_Up; double QE_bf_Ud; double QE_bf_Chi; std::vector<double > E_spectrum;  std::vector<double > A_spectrum;  std::vector<double >QE_spectrum; double stats_bf;} BF_RESULT;

double getTotalNumEvents(CL_input in);

double getTotalNumEvents_NU(CL_input in);
double getTotalNumEvents_NUBAR(CL_input in);

int EtoBin(double E);
double BintoCentralE(int bin);
int CostoBin(double C);
double BintoCentralCos(int b);
int QEtoBin(double QE);
double BintoCentalQE(int b);
int which_BINS(int which);

double QEfromEandCos(double Evis,double costh);

double decayProb(CL_input input, double chiU, double Es);
double histogrammer(CL_input in, double chiU, double cutEff, const double events[][NUM_EVENT_OBS], double eGram[], double cosGram[], double qeGram[]);
double histogrammer_indiv(CL_input in, double chiUp, double cutEff, const double events[][NUM_EVENT_OBS], double Gram[], int which_var);
double histogrammer_indiv2(CL_input in, double Up, double Ud, double chi, double cutEff, const double events[][NUM_EVENT_OBS], double Gram[], int which_var, double finalScale);

int printEGram(double eGram[]);
int printCosGram(double cosGram[]);
int printQeGram(double qeGram[]);

int mcPrintVar(std::vector< double > list);

typedef struct {std::vector<double >  egram; double Sigma_Zeta; int pos_selector;} nuisStruct;
double nuisFuncE(const std::vector<double> &x, std::vector<double> &grad, void *my_data);
double nuisFuncA(const std::vector<double> &x, std::vector<double> &grad, void *my_data);
double nuisFuncQE(const std::vector<double> &x, std::vector<double> &grad, void *my_data);
double nuisMarginalize(std::vector<double > * bf_zeta_b, double * chi, std::vector<double > * eVGram,int whi, double SIGMAZETA,int pos_selector);


double GoF1(std::vector<double > best_spectrum,std::vector<double > spec_obs, std::vector<double >bkg, double bkzeta);
double GoF2(std::vector<double > best_spectrum,std::vector<double > spec_obs, std::vector<double >bkg, double bkzeta);
double pval(double chi2, double ndof);

double intpow( double base, int exponent );

double nuisMarginalize_dual(std::vector<double > * bf_zeta_b, double * chi, std::vector<double > * eVGram,int whi, double SIGMAZETA);
double histogrammer_indiv2_dual(CL_input in, double Up, double Ud, double chi, double cutEff_NU,double cutEff_NUBAR, const double events_NU[][NUM_EVENT_OBS], const double events_NUBAR[][NUM_EVENT_OBS], double Gram[], int which_var, double finalScale_NU, double finalScale_NUBAR);

#endif

