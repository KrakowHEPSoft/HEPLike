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

#include "HEPStats.h"
#include "HEPConstants.h"
#include "HEPBR_ProfLikelihood.h"

using namespace std;

void HEPBR_ProfLikelihood::read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }

  read_standard();

  //cout<<HEPBibEntry<<endl;

  if( config["ROOTData"])  HEPRootFile=config["ROOTData"].as<std::string>();
  else
    {
      std::cout<<"You didn't profice a root file!!! HEPBR_ProfLikelihood class is protesting!"<<std::endl;
    }
  
  if(config["TGraphPath"]) HEPPATH=config["TGraphPath"].as<std::string>();
  
  // now opening files
  f= new TFile(HEPRootFile.c_str(), "READ");
  likelihood=dynamic_cast<TGraph*>(f->Get(HEPPATH.c_str()));
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

double HEPBR_ProfLikelihood::GetChi2(double theory, double theory_err=-1.)
{
  double log_likelihood=GetLogLikelihood(theory,theory_err); 
  double chi2=-2.*log_likelihood; 
  return chi2;
}

double HEPBR_ProfLikelihood::GetLogLikelihood(double theory, double theory_error=-1.)
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
double HEPBR_ProfLikelihood::GetLikelihood(double theory, double theory_err=-1.) 
{
  double log_likelihood=GetLogLikelihood(theory,theory_err);
  return gsl_sf_exp(log_likelihood);  
}




