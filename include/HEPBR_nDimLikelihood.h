//   HEPLike: High Energy Physics Likelihoods
//
//   Header for ndim Likelihood
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//////////////////////////////////////////////////


#ifndef HEPBR_NDIMLIKELIHOOOD_H
#define HEPBR_NDIMLIKELIHOOOD_H

//C++ headers
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>


//HEPLike headers
#include "HEPStats.h"
#include "HEPData.h"

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


class HEPBR_nDimLikelihood: public HEPData
{

 public:

  explicit HEPBR_nDimLikelihood() :  HEPData() {};
  explicit HEPBR_nDimLikelihood(std::string s) :  HEPData(s) { };
  
  
  void read();
  double GetChi2( std::vector<double> theory) ;
  double GetLikelihood( std::vector<double> theory) ;  
  double GetLogLikelihood(  std::vector<double> theory) ;  
  bool Restrict(std::vector<std::string>);
  
 private:


  std::string HEPRootFile;
  std::string HEPPATH;
  std::vector<std::string> Observables;

   
  
  int NoOfObservables;
  int size_restricted;

  TH2D *hist2D;
  TH3D *hist3D;  

  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;
    
  std::vector<double> central_mes_val;
  int dim;
  
  
};


#endif
