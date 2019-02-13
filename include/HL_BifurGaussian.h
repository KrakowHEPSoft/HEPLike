//   HL_Like: High Energy Physics Likelihoods
//
//   Header for construck likelihoods for Bifurcated gaussian distribution
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//////////////////////////////////////////////////
#ifndef HL_BIFURGAUSSIAN_H
#define HL_BIFURGAUSSIAN_H

//C++ headers
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>


//HL_Like headers
#include "HL_Stats.h"
#include "HL_Gaussian.h"

//external:
#include "yaml-cpp/yaml.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_exp.h"

class HL_BifurGaussian: public HL_Gaussian
{

 public:

  explicit HL_BifurGaussian() :  HL_Gaussian() {};
  explicit HL_BifurGaussian(std::string s) :  HL_Gaussian(s) { };
  
  
  void read();
  double GetChi2(double theory, double theory_err=0);
  double GetLikelihood(double theory, double theory_err=0);
  double GetLogLikelihood(double theory, double theory_err=0);  

 private:
  
  double HL_Central;
  double HL_Sigma_stat_R;
  double HL_Sigma_stat_L;
  double HL_Sigma_syst_R;
  double HL_Sigma_syst_L;

  
    

  

  
};


#endif
