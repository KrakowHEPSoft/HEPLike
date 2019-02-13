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

void HL_ProfLikelihood::read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }

  read_standard();

  //cout<<HL_BibEntry<<endl;

  if( config["ROOTData"])  HL_RootFile=config["ROOTData"].as<std::string>();
  else
    {
      std::cout<<"You didn't profice a root file!!! HL_ProfLikelihood class is protesting!"<<std::endl;
    }
  
  if(config["TGraphPath"]) HL_PATH=config["TGraphPath"].as<std::string>();
  
  // now opening files
  f= new TFile(HL_RootFile.c_str(), "READ");
  likelihood=dynamic_cast<TGraph*>(f->Get(HL_PATH.c_str()));
  xmin=likelihood->GetXaxis()->GetXmin () ;
  xmax=likelihood->GetXaxis()->GetXmax () ;
  int N=likelihood->GetN();
  Double_t* x=likelihood->GetX();
  Double_t* y=likelihood->GetY();
  double minX=0;
  double minY=10e10;
  for(int i=2;i<N-2;i++) // exclude the last points
    {
      if(y[i]<minY)
        {
          minY=y[i];
          minX=x[i];
        }
    }
  central_mes_val=minX;


}

double HL_ProfLikelihood::GetChi2(double theory, double theory_err=-1.)
{
  double log_likelihood=GetLogLikelihood(theory,theory_err); 
  double chi2=-2.*log_likelihood; 
  return chi2;
}

double HL_ProfLikelihood::GetLogLikelihood(double theory, double theory_error=-1.)
{
  if(theory < xmin || theory > xmax) return -1.e10;
  if(theory_error<0.){
    double loglikelihood=(-1.)*likelihood->Eval(theory,0, "S" );
    return loglikelihood;
  }
  double loglikelihood=(-1.)*likelihood->Eval(theory,0, "S" );
  double chi2=-2.*loglikelihood;
  // now this is nasty, you should always profile over the theory error but if you want to inflate the chi2:
  double delta=theory-central_mes_val;
  double experimental_err2=chi2/(delta*delta);
  double err2=experimental_err2+theory_error*theory_error;
  chi2=delta*delta/(err2);
  loglikelihood=-0.5*chi2;

  return loglikelihood;
  
}
double HL_ProfLikelihood::GetLikelihood(double theory, double theory_err=-1.) 
{
  double log_likelihood=GetLogLikelihood(theory,theory_err);
  return gsl_sf_exp(log_likelihood);  
}




