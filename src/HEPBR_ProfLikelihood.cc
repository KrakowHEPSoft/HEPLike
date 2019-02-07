//   HEPLike: High Energy Physics Likelihoods
//
//   Module to construck PROFLIKELIHOOD
//
//   author: Marcin Chrzaszcz
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
  cout<<xmin<<" "<<xmax<<endl;
}
/*
double HEPBR_ProfLikelihood::GetChi2(double theory, double theory_err)
{
  double loglikelihood=(-1.)*likelihood->Eval

  double err2=HEPSigma_stat*HEPSigma_stat+ HEPSigma_syst*HEPSigma_syst+theory_err*theory_err;
  double chi2=(HEPCentral-theory)*(HEPCentral-theory)/err2;
  return chi2;
}
*/
double HEPBR_ProfLikelihood::GetLogLikelihood(double theory)
{
  if(theory < xmin || theory > xmax) return 1e10;
  
  double loglikelihood=(-1.)*likelihood->Eval(theory,0, "S" );
  return loglikelihood;

}
/*
double HEPBR_ProfLikelihood::GetLikelihood(double theory, double theory_err) 
{
  double log_likelihood=GetLogLikelihood(theory,theory_err);
  return gsl_sf_exp(log_likelihood);  
}

*/


