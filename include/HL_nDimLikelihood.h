//   HEPLike: High Energy Physics Likelihoods
//
//   Header for ndim Likelihood
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//////////////////////////////////////////////////


#ifndef HL_NDIMLIKELIHOOOD_H
#define HL_NDIMLIKELIHOOOD_H

//C++ headers
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>


//HL_Like headers
#include "HL_Stats.h"
#include "HL_Data.h"

//external:
#include "TFile.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "yaml-cpp/yaml.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_exp.h"
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
//#include <boost/qvm/mat_operations.hpp>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"


/*
class MyFunction: public ROOT::Math::IBaseFunctionMultiDim{

 public:
  double DoEval(const double* theory_nuisance) const{

    bin=hist2D->FindBin(theory_nuisance[0], theory_nuisance[1]);    

    double log_likelihood=-hist2D->GetBinContent(bin); 
    

    
    double gauss_systematic=HL_Stats::gauss(theory_nuisance, theory_mean, theory_err);
    return loglike-log(gauss_systematic);// here the logligek is -\Delta LL so no minus before
  }
  ROOT::Math::IBaseFunctionOneDim* Clone() const{
    return new MyFunction();
  }
  void SetLikelihood(TH2D *l)
  {
    hist2D=l;
  };
  void SetTheory(double *mean, double *cov)
  {
    theory_mean=mean;
    theory_err=err;
  };

 private:
  double theory_mean;
  double theory_err;
  TH2D *hist2D

};


*/







class HL_nDimLikelihood: public HL_Data
{

 public:

  explicit HL_nDimLikelihood() :  HL_Data() {};
  explicit HL_nDimLikelihood(std::string s) :  HL_Data(s) { };
  
  
  void Read();
  double GetChi2( std::vector<double> theory) ;
  double GetLikelihood( std::vector<double> theory) ;  
  double GetLogLikelihood(  std::vector<double> theory) ;  
  //bool Restrict(std::vector<std::string>);
  double GetLogLikelihood(  std::vector<double> theory, std::vector<double> theory_error) ;
  

  
  void Profile();
  double GetChi2_profile( double theory, std::string);
  double GetLikelihood_profile( double theory, std::string axis) ;
  double GetLogLikelihood_profile(  double theory, std::string X);
  
  
  double loglikelihood_penalty;

  
 private:


  std::string HL_RootFile;
  std::string HL_PATH;
  std::vector<std::string> Observables;

  
  int NoOfObservables;
  int size_restricted;

  TH2D *hist2D;
  TH3D *hist3D;  

  TH1 *hist;
  
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;
    
  std::vector<double> central_mes_val;
  int dim;

  TH1D *hist_profileX;
  TH1D *hist_profileY;  
  TH1D *hist_profileZ;

    
  int n_binsX;  
  int n_binsY;
  int n_binsZ;  


  bool profiled;



  
  
};


#endif
