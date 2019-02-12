//   HEPLike: High Energy Physics Likelihoods
//
//   Module to construck likelihoods for Bifurcated gaussian distribution
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
#include "HEPBR_BifurGaussian.h"

using namespace std;

void HEPBR_BifurGaussian::read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }

  read_standard();
  
  HEPCentral=config["Br"].as<double>();
  HEPSigma_stat_L=config["Stat_L"].as<double>();
  HEPSigma_stat_R=config["Stat_R"].as<double>();
  if(config["Syst"] )
    {
      HEPSigma_syst_R=config["Syst"].as<double>();
      HEPSigma_syst_L=config["Syst"].as<double>();   
    }
  else
    {
      HEPSigma_syst_R=config["Syst_R"].as<double>();
      HEPSigma_syst_L=config["Syst_L"].as<double>();
    }
  cout<<HEPSigma_stat_R<<" "<<HEPSigma_stat_L<<" "<<HEPSigma_syst_R<<" "<<HEPSigma_syst_L<<endl;
}
double HEPBR_BifurGaussian::GetChi2(double theory, double theory_err)
{
  double delta=HEPCentral-theory;
  double chi2=0.;
  if(delta>0.)
    {
      double err2=HEPSigma_stat_L*HEPSigma_stat_L+ HEPSigma_syst_L*HEPSigma_syst_L+theory_err*theory_err;    
      chi2=delta*delta/err2;
    }
  else
    {
      double err2=HEPSigma_stat_R*HEPSigma_stat_R+ HEPSigma_syst_R*HEPSigma_syst_R+theory_err*theory_err;
      chi2=delta*delta/err2;
    }
  
  return chi2;

}





