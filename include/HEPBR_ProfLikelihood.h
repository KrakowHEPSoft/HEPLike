//   HEPLike: High Energy Physics Likelihoods
//
//   Header for PROFLIKELIHOOD class
//
//   author: Marcin Chrzaszcz
//////////////////////////////////////////////////
#ifndef HEPBR_PROFLIKELIHOOD_H
#define HEPBR_PROFLIKELIHOOD_H

//C++ headers
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>


//HEPLike headers
#include "HEPStats.h"
#include "HEPData.h"

//external:
#include "yaml-cpp/yaml.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_exp.h"
#include "TFile.h"
#include "TGraph.h"
#include "TAxis.h"

class HEPBR_ProfLikelihood: public HEPData
{

 public:

  explicit HEPBR_ProfLikelihood() :  HEPData() {};
  explicit HEPBR_ProfLikelihood(std::string s) :  HEPData(s) { };
  
  
  void read();
  //double GetChi2(double theory, double theory_err=0);
  double GetLogLikelihood(double theory);//, double theory_err=0);
  //double GetLogLikelihood(double theory, double theory_err=0);  

 private:
  
  double xmin;
  double xmax;
  
  std::string HEPRootFile;
  std::string HEPPATH;
  TGraph *likelihood;
  
  TFile *f;
  
  

  
};


#endif
