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
#include "HEPBRLimit.h"

using namespace std;

void HEPBRLimit::read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }
  read_standard();
  
  YAML::Node node  = config["Cls"];
  for(YAML::const_iterator it = node.begin(); it != node.end();  ++it )
    {
      //std::cout << *it << std::endl;
      //std::cout << (*it)[0] <<  std::endl;   
      CLs.push_back( ((*it)[1]).as<double>()  );
      BR.push_back( ((*it)[0]).as<double>()  ); 
    }

}
double HEPBRLimit::GetChi2(double br)
{
  double cls=GetCLs(br) ;
  //std::cout<<gsl_cdf_gaussian_P(1., 1)-gsl_cdf_gaussian_P(-1., 1.)<<std::endl;
  //double nsigma=0.001;
  //double dsigma=0.0001;
  //double p=0;
  double nsigma=HEPStats::get_sigma_from_pval(1-cls);
  /*
  while (p<1.- cls) {
    p= gsl_sf_erf(nsigma/M_SQRT2);
    nsigma+=dsigma;
  }
  */
  std::cout<<"n of sigmas= "<<nsigma<<"  "<<cls<<std::endl;
  double chi2=nsigma*nsigma;
  //double loglikelihood=gsl_sf_exp(chi2);
  
  return chi2;
}
double HEPBRLimit::GetLogLikelihood(double br)
{
  double chi2=GetChi2(br);
  return -0.5*chi2;
  
}

double HEPBRLimit::GetLikelihood(double br)
{
  double log_likelihood=GetLogLikelihood(br);
  return gsl_sf_exp(log_likelihood);  
}
// this algorithm is O(log n) fast:
double  HEPBRLimit::GetCLs(double value)//, std::vector<double> a, std::vector<double> cls)
{
  int length =BR.size();
  if(value < BR[0]) {
    return CLs[0];
  }
  if(value > BR[length-1]) {
    return CLs[length-1];
  }
  int lo = 0;
  int hi = length - 1;
  while (lo <= hi) {
    int mid = (hi + lo) / 2;
    if (value < BR[mid]) {
      hi = mid - 1;
    }
    else if (value > BR[mid]) {
      lo = mid + 1;
    }
    else {// if we are spot on:
      return CLs[mid];
    }
    
  }// while loop
  // now the lo and hi are next to each other
  // let's use linear interpolation:
  return CLs[lo] + (value - BR[lo])*(CLs[hi]- CLs[lo])/(BR[hi]-BR[lo]); 
  //return (BR[lo] - value) < (value - BR[hi]) ? CLs[lo] : CLs[hi];
}
