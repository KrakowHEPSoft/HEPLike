//   HEPLike: High Energy Physics Likelihoods
//
//   Header for PROFLIKELIHOOD class
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//////////////////////////////////////////////////
#ifndef HL_PROFLIKELIHOOD_H
#define HL_PROFLIKELIHOOD_H

//C++ headers
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>


//HEPLike headers
#include "HL_Stats.h"
#include "HL_Data.h"

//external:
#include "yaml-cpp/yaml.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_exp.h"
#include "TFile.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TMath.h"
#include "TH1D.h"

class HL_ProfLikelihood: public HL_Data
{

 public:

  explicit HL_ProfLikelihood() :  HL_Data() {};
  explicit HL_ProfLikelihood(std::string s) :  HL_Data(s) { };
  
  
  void read();
  double GetChi2(double theory, double theory_err);
  double GetLogLikelihood(double theory, double theory_err);
  double GetLikelihood(double theory, double theory_err);  
  
 private:
  
  double xmin;
  double xmax;
  double central_mes_val;

  
  std::string HL_RootFile;
  std::string HL_PATH;
  TGraph *likelihood;
  
  TFile *f;
  
  

  
};


#endif
