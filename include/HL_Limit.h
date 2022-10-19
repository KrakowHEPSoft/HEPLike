//   HL_Like: High Energy Physics Likelihoods
//
//   Header for limit class
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//////////////////////////////////////////////////
#ifndef HL_LIMIT_H
#define HL_LIMIT_H

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
#include "boost/math/distributions/chi_squared.hpp"


class HL_Limit: public HL_Data
{

 public:

  explicit HL_Limit() :  HL_Data() {};
  explicit HL_Limit(std::string s) :  HL_Data(s) { };


  void Read();
  double GetChi2(double);
  double GetLogLikelihood(double);
  double GetLikelihood(double);
  double GetCLs(double val);

 protected:
  std::vector<double> CLs;
  std::vector<double> BR;
  double UL90CL;
  double UL95CL;
  bool useUL;

  //private:
  //double GetCLs(double val);//, std::vector<double> br, std::vector<double> cls);

};


#endif
