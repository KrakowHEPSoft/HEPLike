//   HEPLike: High Energy Physics Likelihoods
//
//   Module to construct PROFLIKELIHOOD
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
#include "HL_Interpolator.h"
#include "HL_ProfLikelihood.h"

using namespace std;

HL_ProfLikelihood::~HL_ProfLikelihood()
{
  delete likelihood;
  delete gmin;
  delete fun;
}

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
    std::string filename = find_path(config["TextData"].as<std::string>());
    if(debug) std::cout << "Opening file " << filename << std::endl;
    std::ifstream in(filename.c_str());
    in >> nxbins >> xmin >> xmax;
    // Data has entries from 0 to nxbins, so a total of nxbins+1
    nxbins++;

    double x[nxbins];
    double y[nxbins];

    for(int point=0;point<nxbins;point++)
    {
      double point_,x_,y_;
      in >> point_ >> x_>> y_;
      x[point]=x_;
      y[point]=y_;
    }
    in.close();

    likelihood = new HL_Interpolator1D(nxbins,x,y);
    likelihood->SetLimits(xmin,xmax);
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
      likelihood = new HL_Interpolator1D(tmp->Clone());

      tmp->Delete();
      f->Close();
      delete f;

      xmin=likelihood->TG->GetXaxis()->GetXmin();
      xmax=likelihood->TG->GetXaxis()->GetXmax();
      likelihood->SetLimits(xmin,xmax);

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

  gmin = new HL_Minimizer("ConjugateFR", 1);
  gmin->SetMaxIterations(1000000);
  gmin->SetTolerance(0.0001);

  fun = new HL_Function1D(likelihood);
}

double HL_ProfLikelihood::GetLogLikelihood(double theory)
{
  if(theory < xmin || theory > xmax) return -1.e10;
  double loglikelihood=(-1)*likelihood->Eval(theory);

  return loglikelihood;
}

double HL_ProfLikelihood::GetLogLikelihood(double theory, double theory_variance)
{
  double theory_err=sqrt(theory_variance);

  if(theory < xmin || theory > xmax) return -1.e10;
  fun->SetTheory(theory,theory_err);
  //ROOT::Math::Functor  f1(fun,1);

  //gmin->SetFunction(f1);
  //double step[1] = {0.01*theory_err};
  //double variable[1]={theory-5.*theory_err};
  //gmin->SetVariable(0,"x",variable[0], step[0]);
  //gmin->SetVariableInitialRange(0,theory-5*theory_err, theory+5.*theory_err);
  //gmin->SetVariableLimits(0,theory-5.*theory_err, theory+5.*theory_err);
  //gmin->Minimize();
  //const double *theory_nuisance = gmin->X();

  //double loglikelihood=likelihood->Eval(theory_nuisance[0],0);

  gmin->SetFunction(fun);
  double step[1] = {0.01*theory_err};
  double variable[1]={theory-5.*theory_err};
  gmin->SetVariable(0, variable[0], step[0]);
  gmin->Minimize();
  std::vector<double> theory_nuisance = gmin->X();

  double loglikelihood = likelihood->Eval(theory_nuisance[0]);

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
