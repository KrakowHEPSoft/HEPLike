//   HEPLike: High Energy Physics Likelihoods
//
//   Header for Experimental Data Points class
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//////////////////////////////////////////////////
#ifndef HL_EXPPOINTS_H
#define HL_EXPPOINTS_H

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


class HL_ExpPoints: public HL_Data
{

 public:

  explicit HL_ExpPoints() :  HL_Data() {};
  explicit HL_ExpPoints(std::string s) :  HL_Data(s) { };


  void Read();
  double GetChi2(vector<double> theory);
  double GetLogLikelihood(vector<double> theory);
  double GetLikelihood(vector<double> theory);
  bool InitData();
  void SetFun( double FUN( vector<double>, vector<double>));


 private:

  TTree *HL_tree;
  double xmin;
  double xmax;
  double central_mes_val;
  std::string ObsName;

  std::string HL_RootFile;
  std::string HL_PATH;
  std::string HL_weight;
  vector<std::string> HL_obs;
  vector<TBranch*> HL_Branches;
  TBranch *HL_weight_branch;


  TFile *f;

  double (*fun)(vector<double>, vector<double>);

  vector<vector<double>> points;
  vector<double> weights;









};


#endif
