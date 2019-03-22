//   HL_Like: High Energy Physics Likelihoods
//
//   Module to read yaml files
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
#include "HL_Limit.h"

using namespace std;

void HL_Limit::Read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }
  read_standard();
  
  YAML::Node node;
  if( config["Cls"]) node= config["Cls"];
  else if(  config["p-value"] )  node= config["p-value"] ;

  for(YAML::const_iterator it = node.begin(); it != node.end();  ++it )
    {
      //std::cout << *it << std::endl;
      //std::cout << (*it)[0] <<  std::endl;   
      CLs.push_back( ((*it)[1]).as<double>()  );
      BR.push_back( ((*it)[0]).as<double>()  ); 
    }

};
double HL_Limit::GetChi2(double br)
{
  double cls=GetCLs(br) ;
  //std::cout<<gsl_cdf_gaussian_P(1., 1)-gsl_cdf_gaussian_P(-1., 1.)<<std::endl;
  //double nsigma=0.001;
  //double dsigma=0.0001;
  //double p=0;
  double chi2=gsl_cdf_chisq_Pinv(1-cls,1);

  return chi2;
}
double HL_Limit::GetLogLikelihood(double br)
{
  double chi2=GetChi2(br);
  return -0.5*chi2;
  
}
double HL_Limit::GetLikelihood(double br)
{
  double log_likelihood=GetLogLikelihood(br);
  return gsl_sf_exp(log_likelihood);  
}
// this algorithm is O(log n) fast:
double  HL_Limit::GetCLs(double value)//, std::vector<double> a, std::vector<double> cls)
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
