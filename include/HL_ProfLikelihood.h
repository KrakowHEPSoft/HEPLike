//   HEPLike: High Energy Physics Likelihoods
//
//   Header for PROFLIKELIHOOD class
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//   author: Tomas Gonzalo
//////////////////////////////////////////////////
#ifndef HL_PROFLIKELIHOOD_H
#define HL_PROFLIKELIHOOD_H

//C++ headers
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>


//HEPLike headers
#include "HL_Data.h"
#include "HL_Root.h"
#include "HL_Interpolator.h"
#include "HL_Function.h"
#include "HL_Minimizer.h"

//external:
#include "yaml-cpp/yaml.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_exp.h"


class HL_ProfLikelihood: public HL_Data
{

 public:

  explicit HL_ProfLikelihood() : HL_Data(), likelihood(nullptr), gmin(nullptr), fun(nullptr) {};
  explicit HL_ProfLikelihood(std::string s) :  HL_Data(s), likelihood(nullptr), gmin(nullptr), fun(nullptr){};

  HL_ProfLikelihood(const HL_ProfLikelihood &);
  HL_ProfLikelihood &operator=(const HL_ProfLikelihood &);

  ~HL_ProfLikelihood();

  void Read();
  double GetChi2(double theory);
  double GetChi2(double theory, double theory_err);

  double GetLogLikelihood(double theory);
  double GetLogLikelihood(double theory, double theory_err);

  double GetLikelihood(double theory);
  double GetLikelihood(double theory, double theory_err);

  int nxbins;
  double xmin;
  double xmax;
  double central_mes_val;
  std::string ObsName;

  std::string HL_RootFile;
  std::string HL_PATH;

  HL_Interpolator1D *likelihood;

  HL_Minimizer *gmin;

  HL_Function1D *fun;


};


#endif
