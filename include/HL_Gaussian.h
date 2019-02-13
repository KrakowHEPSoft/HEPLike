//   HEPLike: High Energy Physics Likelihoods
//
//   Header for Gaussian likelihood
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//////////////////////////////////////////////////
#ifndef HL_GAUSSIAN_H
#define HL_GAUSSIAN_H

//C++ headers
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>


//HL_Like headers
#include "HL_Stats.h"
#include "HL_Data.h"

//external:
#include "yaml-cpp/yaml.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_exp.h"

class HL_Gaussian: public HL_Data
{

 public:

  explicit HL_Gaussian() :  HL_Data() {};
  explicit HL_Gaussian(std::string s) :  HL_Data(s) { };
  
  
  void read();
  double GetChi2(double theory, double theory_err=0);
  double GetLikelihood(double theory, double theory_err=0);
  double GetLogLikelihood(double theory, double theory_err=0);  

 private:
  
  double HL_Central;
  double HL_Sigma_stat;
  double HL_Sigma_syst;
  std::string ObsName;
  
  

  
};


#endif
