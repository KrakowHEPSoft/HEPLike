//   HL_Like: High Energy Physics Likelihoods
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

#include "HL_Stats.h"
#include "HL_Constants.h"
#include "HL_BifurGaussian.h"

using namespace std;

void HL_BifurGaussian::read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }

  read_standard();
  
  HL_Central=config["Br"].as<double>();
  HL_Sigma_stat_L=config["Stat_L"].as<double>();
  HL_Sigma_stat_R=config["Stat_R"].as<double>();
  if(config["Syst"] )
    {
      HL_Sigma_syst_R=config["Syst"].as<double>();
      HL_Sigma_syst_L=config["Syst"].as<double>();   
    }
  else
    {
      HL_Sigma_syst_R=config["Syst_R"].as<double>();
      HL_Sigma_syst_L=config["Syst_L"].as<double>();
    }
  cout<<HL_Sigma_stat_R<<" "<<HL_Sigma_stat_L<<" "<<HL_Sigma_syst_R<<" "<<HL_Sigma_syst_L<<endl;
}
double HL_BifurGaussian::GetChi2(double theory, double theory_err)
{
  double delta=HL_Central-theory;
  double chi2=0.;
  if(delta>0.)
    {
      double err2=HL_Sigma_stat_L*HL_Sigma_stat_L+ HL_Sigma_syst_L*HL_Sigma_syst_L+theory_err*theory_err;    
      chi2=delta*delta/err2;
    }
  else
    {
      double err2=HL_Sigma_stat_R*HL_Sigma_stat_R+ HL_Sigma_syst_R*HL_Sigma_syst_R+theory_err*theory_err;
      chi2=delta*delta/err2;
    }
  
  return chi2;

}





