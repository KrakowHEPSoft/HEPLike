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
#include "HL_Root.h"

//external:
#include "yaml-cpp/yaml.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_exp.h"


class MyFunction: public ROOT::Math::IBaseFunctionOneDim{

 public:
  double DoEval(double theory_nuisance) const{

    double loglike=likelihood->Eval(theory_nuisance,0);
    double like=exp(loglike);


    double gauss_systematic=HL_Stats::gauss(theory_nuisance, theory_mean, theory_err);

    return loglike-log(gauss_systematic);// here the logligek is -\Delta LL so no minus before
  }
  ROOT::Math::IBaseFunctionOneDim* Clone() const{
    return new MyFunction();
  }
  void SetLikelihood(TGraph *l)
    {
      likelihood=l;
    };
  void SetTheory(double mean, double err)
  {
    theory_mean=mean;
    theory_err=err;
  };

 private:
  double theory_mean;
  double theory_err;
  TGraph *likelihood;

};




class HL_ProfLikelihood: public HL_Data
{

 public:

  explicit HL_ProfLikelihood() :  HL_Data() {};
  explicit HL_ProfLikelihood(std::string s) :  HL_Data(s) { };


  void Read();
  double GetChi2(double theory);
  double GetChi2(double theory, double theory_err);

  double GetLogLikelihood(double theory);
  double GetLogLikelihood(double theory, double theory_err);

  double GetLikelihood(double theory);
  double GetLikelihood(double theory, double theory_err);





 private:

  double xmin;
  double xmax;
  double central_mes_val;
  std::string ObsName;

  std::string HL_RootFile;
  std::string HL_PATH;
  TGraph *likelihood;

  // for minimaization
  ROOT::Math::Minimizer* gmin;


  TFile *f;

  MyFunction fun;


};


#endif
