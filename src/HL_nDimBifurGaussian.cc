//   HL_Like: High Energy Physics Likelihoods
//
//   Module to construck likelihoods for ndim bifurcated gaussian distribution
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
#include "HL_nDimBifurGaussian.h"

using namespace std;

static bool ndbg_debug = false;

void HL_nDimBifurGaussian::Read()
{
  if(! initialized)
  {
    throw std::runtime_error("HL_nDimBifurGaussian Warninig, TRYING TO READ WITHOUT GIVING ANY FILE!");
    return;
  }
  read_standard();
  restricted=false;
  if(config["Observables"])
  {
    YAML::Node node  = config["Observables"];

    for(YAML::const_iterator it = node.begin(); it != node.end();  ++it )
    {
      Observables.push_back( ((*it)[0]).as<std::string>()  );
      central.push_back( ((*it)[1]).as<double>()  );

      stat_error_right.push_back( ((*it)[2]).as<double>()  );
      stat_error_left.push_back( (-1.)* ((*it)[3]).as<double>()  );
      syst_error_right.push_back( ((*it)[4]).as<double>()  );
      if( (*it).size()>5 )
      {
        syst_error_left.push_back( (-1.)*((*it)[5]).as<double>()  );

      }
      else syst_error_left.push_back( ((*it)[4]).as<double>()  );

    }
  }// read the errors and cenral vaules, now correlation


  NoOfObservables=Observables.size();

  for(int i=0 ; i<  NoOfObservables; i++)
  {
    error_right.push_back( sqrt(stat_error_right[i]*stat_error_right[i] + syst_error_right[i]*syst_error_right[i]));
    error_left.push_back( sqrt(stat_error_left[i]*stat_error_left[i] + syst_error_left[i]*syst_error_left[i]));
  }

  if(config["Correlation"])
  {
    YAML::Node node  = config["Correlation"];
    int row=0;

    HL_correlation= boost::numeric::ublas::matrix<double>(NoOfObservables,NoOfObservables);
    HL_cov = boost::numeric::ublas::matrix<double>(NoOfObservables,NoOfObservables);
    HL_cov_inv = boost::numeric::ublas::matrix<double>(NoOfObservables,NoOfObservables);
    for(YAML::const_iterator it = node.begin(); it != node.end();  ++it )
    {
      if(row==0)
      {
        for(int i=0;i<NoOfObservables; i++)
        {
          if(Observables[i] != ((*it)[i]).as<std::string>()  )
          {
            throw std::runtime_error("Error in HL_nDimBifurGaussian, wrong ordering of names in correlation matrix");
          }
        }
      }// row=0
      else
      {
        for(int i=0;i<NoOfObservables; i++)
        {
          HL_correlation(row-1,i)= ((*it)[i]).as<double>();
        }
      }

      row++;
    }
  }
  else // no correlation
  {
    HL_correlation= boost::numeric::ublas::matrix<double>(NoOfObservables,NoOfObservables);
    HL_cov = boost::numeric::ublas::matrix<double>(NoOfObservables,NoOfObservables);
    HL_cov_inv = boost::numeric::ublas::matrix<double>(NoOfObservables,NoOfObservables);

    for(int i=0;i<NoOfObservables; i++)
    {
      for(int j=0;j<NoOfObservables; j++)
      {
        if(i==j)
        {
          HL_correlation(i,j) =1.;
          HL_cov(i,j) =1.;
          HL_cov_inv(i,j)=1.;
        }//diagonal
        else
        {
          HL_correlation(i,j) =0.;
          HL_cov(i,j) =0.;
          HL_cov_inv(i,j)=0.;
        }//ofdiagonal
      }// j
    }// i

  }//no correlation case
  size_restricted=NoOfObservables;

}

bool HL_nDimBifurGaussian::Restrict(std::vector<std::string> names)
{
  size_restricted=  names.size();
  std::vector<int> indexes;
  for(int i=0 ; i<NoOfObservables ; i++)
  {
    std::string element=Observables[i];
    if(std::find(names.begin(), names.end(), element) != names.end())
    {
      // names have the element:
      indexes.push_back(i);
    }
  }
  HL_cov_restricted= boost::numeric::ublas::matrix<double>(size_restricted, size_restricted);
  HL_correlation_restricted=  boost::numeric::ublas::matrix<double>(size_restricted, size_restricted);
  HL_cov_inv_restricted=  boost::numeric::ublas::matrix<double>(size_restricted, size_restricted);

  for(int i=0; i< size_restricted; i++)
  {
    central_restricted.push_back(central [ indexes[i] ] );
    for(int j=0; j< size_restricted; j++)
    {
      int i_index=indexes[i];
      int j_index=indexes[j];
      HL_correlation_restricted(i,j)=HL_correlation(i_index, j_index);
    }
  }
  restricted=true;

  return true;
}


double HL_nDimBifurGaussian::GetChi2(std::vector<double> theory)  //, double theory_err)
{
  if(!restricted) // if we don't resctric and use whole matrix
  {
    HL_correlation_restricted=HL_correlation;
    //HL_Stats::inverse(HL_correlation_restricted, size_restricted);
    size_restricted=NoOfObservables;
    HL_Stats::InvertMatrix(HL_cov_restricted,HL_cov_inv_restricted);
    //      restricted=true;
    central_restricted=central;
  }

  // now calculating chi2
  double chi2=0;
  vector<double> diff;
  if(theory.size() !=  central_restricted.size())
  {
    if(ndbg_debug) std::cout<<"Theory: "<<theory.size()<<"  Exp:"<< central_restricted.size()<<"  "<< central.size() <<std::endl;
    throw std::runtime_error("Error in HL_nDimBifurGaussian::GetChi2, you had different dimensions in theory and experiment");
    return -1e10;
  }
  for(size_t i=0;i<central_restricted.size(); i++) { diff.push_back(central_restricted[i] - theory[i] );}


  HL_cov_restricted=  boost::numeric::ublas::matrix<double>(size_restricted, size_restricted);
  HL_cov_inv_restricted=  boost::numeric::ublas::matrix<double>(size_restricted, size_restricted);

  for(int i=0; i< size_restricted; i++)
  {
    for(int j=0; j< size_restricted; j++)
    {
      if(diff[i] >= 0. && diff[j] >= 0.)      HL_cov_restricted(i,j)=error_right[i]*error_right[j]*HL_correlation_restricted(i,j);
      if(diff[i] >= 0. && diff[j] < 0.)      HL_cov_restricted(i,j)=error_right[i]*error_left[j]*HL_correlation_restricted(i,j);
      if(diff[i] < 0. && diff[j] >= 0.)      HL_cov_restricted(i,j)=error_left[i]*error_right[j]*HL_correlation_restricted(i,j);
      if(diff[i] < 0. && diff[j] < 0.)      HL_cov_restricted(i,j)=error_left[i]*error_left[j]*HL_correlation_restricted(i,j);
    }
  }

  HL_Stats::InvertMatrix(HL_cov_restricted,HL_cov_inv_restricted);
  for (size_t i=0; i < HL_cov_inv_restricted.size1(); ++i)
  {
    for (size_t j=0; j<HL_cov_inv_restricted.size2(); ++j)
    {
      chi2+= diff[i] * HL_cov_inv_restricted(i,j)*diff[j] ;
    }
  }
  if (ndbg_debug) std::cout<<"Returning chi2: "<<chi2<<std::endl;
  return chi2;
}

double HL_nDimBifurGaussian::GetLogLikelihood(std::vector<double> theory)
{
  double chi2=GetChi2(theory);

  return -0.5*chi2;
}

double HL_nDimBifurGaussian::GetLikelihood(std::vector<double> theory)
{
  double log_likelihood=GetLogLikelihood(theory);
  return gsl_sf_exp(log_likelihood);
}


double HL_nDimBifurGaussian::GetChi2(std::vector<double> theory,  boost::numeric::ublas::matrix<double> theory_cov)
{
  if(theory_cov.size1() != theory.size() )
  {
    throw std::runtime_error("Error in HL_nDimBifurGaussian::GetChi2, you had different dimensions in theory and cov matrix");
  }
  if(theory_cov.size2() != theory_cov.size1() )
  {
    throw std::runtime_error("Error in HL_nDimBifurGaussian::GetChi2, your theory cov matrix is not square!");
  }


  if(!restricted) // if we don't restrict and use whole matrix
  {
    HL_correlation_restricted=HL_correlation;
    size_restricted=NoOfObservables;
    HL_Stats::InvertMatrix(HL_cov_restricted,HL_cov_inv_restricted);
      //restricted=true;
    central_restricted=central;
  }
  // now calculating chi2
  double chi2=0;
  vector<double> diff;
  if(theory.size() !=  central_restricted.size())
  {
    if(ndbg_debug) std::cout<<"Theory: "<<theory.size()<<"  Exp:"<< central_restricted.size()<<" "<<  central.size()<<std::endl;
    throw std::runtime_error("Error in HL_nDimBifurGaussian::GetChi2, you had different dimensions in theory and experiment");
    return -1e10;
  }
  for(size_t i=0;i<central_restricted.size(); i++) { diff.push_back(central_restricted[i] - theory[i] );}


  HL_cov_restricted=  boost::numeric::ublas::matrix<double>(size_restricted, size_restricted);
  HL_cov_inv_restricted=  boost::numeric::ublas::matrix<double>(size_restricted, size_restricted);

  for(int i=0; i< size_restricted; i++)
  {
    for(int j=0; j< size_restricted; j++)
    {
      if(diff[i] >= 0. && diff[j] >= 0.)      HL_cov_restricted(i,j)=error_right[i]*error_right[j]*HL_correlation_restricted(i,j);
      if(diff[i] >= 0. && diff[j] < 0.)      HL_cov_restricted(i,j)=error_right[i]*error_left[j]*HL_correlation_restricted(i,j);
      if(diff[i] < 0. && diff[j] >= 0.)      HL_cov_restricted(i,j)=error_left[i]*error_right[j]*HL_correlation_restricted(i,j);
      if(diff[i] < 0. && diff[j] < 0.)      HL_cov_restricted(i,j)=error_left[i]*error_left[j]*HL_correlation_restricted(i,j);
    }
  }

  //adding theory error:
  HL_cov_restricted+=theory_cov;

  HL_Stats::InvertMatrix(HL_cov_restricted,HL_cov_inv_restricted);

  for (size_t i=0; i < HL_cov_inv_restricted.size1(); ++i)
  {
    for (size_t j=0; j<HL_cov_inv_restricted.size2(); ++j)
    {
      chi2+= diff[i] * HL_cov_inv_restricted(i,j)*diff[j] ;
    }
  }
  return chi2;
}

double HL_nDimBifurGaussian::GetLogLikelihood(std::vector<double> theory, boost::numeric::ublas::matrix<double> theory_cov)
{
  double chi2=GetChi2(theory, theory_cov);
  return -0.5*chi2;
}

double HL_nDimBifurGaussian::GetLikelihood(std::vector<double> theory, boost::numeric::ublas::matrix<double> theory_cov)
{
  double log_likelihood=GetLogLikelihood(theory,  theory_cov);
  return gsl_sf_exp(log_likelihood);
}

