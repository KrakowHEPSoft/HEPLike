//   HEPLike: High Energy Physics Likelihoods
//
//   Module to read yaml files
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
#include "HEPBR_Gaussian.h"

using namespace std;

void HEPBR_Gaussian::read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }
  YAML::Node config = YAML::LoadFile(HFile);

  HEPDOI=config["DOI"].as<std::string>();
  HEPBibCite=config["BibCite"].as<std::string>();
  HEPBibEntry=config["BibEntry"].as<std::string>();
  HEPFileName=config["FileName"].as<std::string>();
  HEPHFLAV=config["HFLAV"].as<std::string>();
  HEPSource=config["Source"].as<std::string>();
  HEPYear=config["Year"].as<std::string>();
  HEPName=config["Name"].as<std::string>();
  HEPDecay=config["Decay"].as<std::string>();
  
  HEPCentral=config["BR"].as<double>();
  HEPSigma_stat=config["Stat"].as<double>();
  HEPSigma_syst=config["Syst"].as<double>();
  
  
}
/*
double HEPBR_Gaussian::GetLogLikelihood(double br)
{
  double cls=getCLs(br) ;
  //std::cout<<gsl_cdf_gaussian_P(1., 1)-gsl_cdf_gaussian_P(-1., 1.)<<std::endl;
  double nsigma=0.001;
  double dsigma=0.0001;

  double p=0;
  while (p<1.- cls) {
    p= gsl_sf_erf(nsigma/M_SQRT2);
    nsigma+=dsigma;
  }
  //std::cout<<"n of sigmas= "<<nsigma<<"  "<<cls<<std::endl;
  double chi2=nsigma*nsigma;
  //double loglikelihood=gsl_sf_exp(chi2);
  
  return -0.5*chi2;
}
double HEPBR_Gaussian::GetLikelihood(double br)
{
  double log_likelihood=GetLogLikelihood(br);
  return gsl_sf_exp(log_likelihood);  
}

*/



