//   HEPLike: High Energy Physics Likelihoods
//
//   Header for construck likelihoods for Bifurcated gaussian distribution
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//////////////////////////////////////////////////
#ifndef HEPBR_BIFURGAUSSIAN_H
#define HEPBR_BIFURGAUSSIAN_H

//C++ headers
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>


//HEPLike headers
#include "HEPStats.h"
#include "HEPBR_Gaussian.h"

//external:
#include "yaml-cpp/yaml.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_exp.h"

class HEPBR_BifurGaussian: public HEPBR_Gaussian
{

 public:

  explicit HEPBR_BifurGaussian() :  HEPBR_Gaussian() {};
  explicit HEPBR_BifurGaussian(std::string s) :  HEPBR_Gaussian(s) { };
  
  
  void read();
  double GetChi2(double theory, double theory_err=0);
  double GetLikelihood(double theory, double theory_err=0);
  double GetLogLikelihood(double theory, double theory_err=0);  

 private:
  
  double HEPCentral;
  double HEPSigma_stat_R;
  double HEPSigma_stat_L;
  double HEPSigma_syst_R;
  double HEPSigma_syst_L;

  
    

  

  
};


#endif
