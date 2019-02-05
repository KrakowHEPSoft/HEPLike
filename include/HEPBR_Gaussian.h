//   HEPLike: High Energy Physics Likelihoods
//
//   Header for limit class
//
//   author: Marcin Chrzaszcz
//////////////////////////////////////////////////
#ifndef HEPBR_GAUSSIAN_H
#define HEPBR_GAUSSIAN_H

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

class HEPBR_Gaussian: public HEPData
{

 public:

  explicit HEPBR_Gaussian() :  HEPData() {};
  explicit HEPBR_Gaussian(std::string s) :  HEPData(s) { };
  
  
  void read();
  double GetLogLikelihood(double);
  double GetLikelihood(double);   

 private:
  
  double HEPCentral;
  double HEPSigma_stat;
  double HEPSigma_syst;
    

  

  
};


#endif
