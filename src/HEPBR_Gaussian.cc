//   HEPLike: High Energy Physics Likelihoods
//
//   Module to construck likelihoods for gaussian distribution
//
//   author: Jihyun Bhom,  Marcin Chrzaszcz
//////////////////////////////////////////////////


#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "HEPStats.h"
#include "HEPConstants.h"
#include "HEPBR_Gaussian.h"

using namespace std;

void HEPBR_Gaussian::read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }

  read_standard();
  
  HEPCentral=config["Br"].as<double>();
  HEPSigma_stat=config["Stat"].as<double>();
  HEPSigma_syst=config["Syst"].as<double>();
  
}
double HEPBR_Gaussian::GetChi2(double theory, double theory_err)
{
  double err2=HEPSigma_stat*HEPSigma_stat+ HEPSigma_syst*HEPSigma_syst+theory_err*theory_err;
  double chi2=(HEPCentral-theory)*(HEPCentral-theory)/err2;
  return chi2;
}


double HEPBR_Gaussian::GetLogLikelihood(double theory, double theory_err)  
{

  double chi2=GetChi2(theory,theory_err);
  
  return -0.5*chi2;
}
double HEPBR_Gaussian::GetLikelihood(double theory, double theory_err) 
{
  double log_likelihood=GetLogLikelihood(theory,theory_err);
  return gsl_sf_exp(log_likelihood);  
}




