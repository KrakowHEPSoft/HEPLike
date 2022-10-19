//   HEPLike: High Energy Physics Likelihoods
//
//   Module to construck PROFLIKELIHOOD
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//////////////////////////////////////////////////


#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "HL_Stats.h"
#include "HL_Constants.h"

#include "HL_ProfLikelihood.h"

using namespace std;

bool debug = false;

void HL_ProfLikelihood::Read()
{

  // Note: amended to allow text file input
  if(! initialized)
   {
    std::cout << "HL_ProfLikelihood Warning, TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
    return;
  }

  read_standard();

  //Always use text input if it is available (avoids potential mem issues)
  if(config["TextData"])
  {
    if(debug) std::cout << "HL_ProfLikelihood is using text input" << std::endl;
    std::filename = find_filename(config["TextData"].as<std::string>);
    if(debug) std::cout << "Opening file " << filename << std::endl;
    std::ifstream in(filename.c_str());
    int nxbins;
    in >> nxbins >> xmin >> xmax;

    double x[nxbins];
    double y[nxbins];

    for(int point=0;point<nxbins+1;point++)
    {
      double point_,x_,y_;
      in >> point_ >> x_>> y_;
      x[point]=x_;
      y[point]=y_;
    }
    in.close();

    likelihood = LikelihoodInterpolator(nxbins,x,y);
    likelihood.SetLimits(xmin,xmax);
  }
  else if( config["ROOTData"])
  {
    #ifdef USE_ROOT
      // the HL_RootFile is something like data/RD/... but we need the absolute path:
      HL_RootTile = find_path(config["ROOTData"].as<std::string>();
      if(config["TGraphPath"]) HL_PATH=config["TGraphPath"].as<std::string>();

      // now opening files
      Tfile *f= new TFile(HL_RootFile.c_str(), "READ");
      TGraph *tmp=dynamic_cast<TGraph*>(f->Get(HL_PATH.c_str()));
      likelihood=dynamic_cast<TGraph*>(tmp->Clone());

      tmp->Delete();
      f->Close();
      delete f;

      xmin=likelihood->GetXaxis()->GetXmin () ;
      xmax=likelihood->GetXaxis()->GetXmax () ;

    #else
      throw std::runtime_error("Requested ROOT file but ROOT is not enabled. Please enable ROOT or provide data in text format.");
    #endif
  }
  else
  {
    throw std::runtime_error("You didn't provide neither a text nor a ROOT file!!! HL_ProfLikelihood class is protesting!");
  }

  if(config["Observables"] )
  {
    YAML::Node node  = config["Observables"];
    ObsName=node[0][0].as<std::string>();
    central_mes_val=node[0][1].as<double>();
  }

  //gmin=ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "ConjugateFR");
  //gmin->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
  //gmin->SetMaxIterations(1000000);  // for GSL
  //gmin->SetTolerance(0.0001);
  //gmin->SetPrintLevel(3);
  gmin = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, 1);
  double niter = 1000000;
  double tolerance = 0.0001;

  fun=MyFunction();
  fun.SetLikelihood(likelihood);

}

double HL_ProfLikelihood::GetLogLikelihood(double theory)
{
  if(theory < xmin || theory > xmax) return -1.e10;
  double loglikelihood=(-1)*likelihood->Eval(theory,0);

  return loglikelihood;
}

double HL_ProfLikelihood::GetLogLikelihood(double theory, double theory_variance)
{
  double theory_err=sqrt(theory_variance);

  if(theory < xmin || theory > xmax) return -1.e10;
  fun.SetTheory(theory,theory_err);
  ROOT::Math::Functor  f1(fun,1);


  gmin->SetFunction(f1);
  double step[1] = {0.01*theory_err};
  double variable[1]={theory-5.*theory_err};
  gmin->SetVariable(0,"x",variable[0], step[0]);
  gmin->SetVariableInitialRange(0,theory-5*theory_err, theory+5.*theory_err);
  gmin->SetVariableLimits(0,theory-5.*theory_err, theory+5.*theory_err);
  gmin->Minimize();
  const double *theory_nuisance = gmin->X();

  double loglikelihood=likelihood->Eval(theory_nuisance[0],0);

  return (-1.)*loglikelihood;
}


double HL_ProfLikelihood::GetChi2(double theory)
{
  double log_likelihood=GetLogLikelihood(theory);
  double chi2=-2.*log_likelihood;
  return chi2;
}

double HL_ProfLikelihood::GetChi2(double theory, double theory_err)
{
  double log_likelihood=GetLogLikelihood(theory, theory_err);
  double chi2=-2.*log_likelihood;
  return chi2;
}


double HL_ProfLikelihood::GetLikelihood(double theory)
{
  double log_likelihood=GetLogLikelihood(theory);
  return gsl_sf_exp(log_likelihood);
}

double HL_ProfLikelihood::GetLikelihood(double theory, double theory_err)
{
  double log_likelihood=GetLogLikelihood(theory, theory_err);
  return gsl_sf_exp(log_likelihood);
}
