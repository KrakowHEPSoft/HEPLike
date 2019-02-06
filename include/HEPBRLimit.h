//   HEPLike: High Energy Physics Likelihoods
//
//   Header for limit class
//
//   author: Marcin Chrzaszcz
//////////////////////////////////////////////////
#ifndef HEPBRLIMIT_H
#define HEPBRLIMIT_H

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
class HEPBRLimit: public HEPData
{

 public:

  explicit HEPBRLimit() :  HEPData() {};
  explicit HEPBRLimit(std::string s) :  HEPData(s) { };
  
  
  void read();
  double GetChi2(double);     
  double GetLogLikelihood(double);
  double GetLikelihood(double);   

 protected:
  std::vector<double> CLs;
  std::vector<double> BR;
  
 private:
  double getCLs(double val);//, std::vector<double> br, std::vector<double> cls);
  

  
};


#endif
