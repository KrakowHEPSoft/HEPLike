//   HL_Like: High Energy Physics Likelihoods
//
//   Module to read yaml files
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//   author: Tomas Gonzalo
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
    throw std::runtime_error("HL_Limit, warninig, TRYING TO READ WITHOUT GIVING ANY FILE!");
  }
  read_standard();
  useUL=false;
  UL90CL=-1.;
  UL95CL=-1.;
  YAML::Node node;
  // case where we have a p-value scan:
  if( config["Cls"]) node= config["Cls"];
  else if(  config["p-value"] )  node= config["p-value"] ;
  if(node)
  {
    for(YAML::const_iterator it = node.begin(); it != node.end();  ++it )
    {
      CLs.push_back( ((*it)[1]).as<double>()  );
      BR.push_back( ((*it)[0]).as<double>()  );
    }
  }
  // case where we have just UL:
  if(config["UL90CL"])
  {
    UL90CL=config["UL90CL"].as<double>();
    useUL=true;
  }
  if(config["UL95CL"])
  {
    UL95CL=config["UL95CL"].as<double>();
    useUL=true;
  }

};

double HL_Limit::GetChi2(double br)
{
  double chi2;

  if(useUL== false)
  {
    double cls=GetCLs(br) ;
    chi2=gsl_cdf_chisq_Pinv(1-cls,1);
  }
  else
  {
    double error;
    if(UL95CL>0.)
    {
      error=fabs(0. - UL95CL)/1.96;
    }
    else if(UL90CL>0.)
    {
      error = fabs(0. - UL90CL)/1.64;
    }
    double loglikelihood=HL_Stats::gaussian_upper_limit(br, 0., 0., error, false);

    chi2=loglikelihood*(-2.);

  }
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
  if(value < BR[0])
    return CLs[0];
  if(value > BR[length-1])
    return CLs[length-1];
  int lo = 0;
  int hi = length - 1;
  while (lo <= hi)
  {
    int mid = (hi + lo) / 2;
    if (value < BR[mid])
      hi = mid - 1;
    else if (value > BR[mid])
      lo = mid + 1;
    else // if we are spot on:
      return CLs[mid];

  }// while loop
  // now the lo and hi are next to each other
  // let's use linear interpolation:
  return CLs[lo] + (value - BR[lo])*(CLs[hi]- CLs[lo])/(BR[hi]-BR[lo]);
  //return (BR[lo] - value) < (value - BR[hi]) ? CLs[lo] : CLs[hi];
}
