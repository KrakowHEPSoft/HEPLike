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
