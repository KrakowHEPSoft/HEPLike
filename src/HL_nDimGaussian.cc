//   HEPLike: High Energy Physics Likelihoods
//
//   Module to construck likelihoods for ndim gaussian distribution
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
#include "HL_nDimGaussian.h"


using namespace std;



void HL_nDimGaussian::Read()
{
  if(! initialized)
    {
      std::cout << "HL_nDimGaussian, TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }
  read_standard();
  if(config["Observables"])
    {
      YAML::Node node  = config["Observables"];

      for(YAML::const_iterator it = node.begin(); it != node.end();  ++it )
        {
          Observables.push_back( ((*it)[0]).as<std::string>()  );
          central.push_back( ((*it)[1]).as<double>()  );
          stat_error.push_back( ((*it)[2]).as<double>()  );
          //cout<<((*it)[0]).as<std::string>()<<" "<<((*it)[1]).as<double>()<<" "<< ((*it)[2]).as<double>()<<endl;
          if( (*it).size()>3 )
            {
              syst_error.push_back( ((*it)[3]).as<double>()  );

            }
          else
            {
              syst_error.push_back( 0.);
            }
        }
    }// read the errors and cenral vaules, now correlation
  NoOfObservables=Observables.size();

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
                      cout<<"Error in HL_nDimGaussian, wrong ordering of names in correlation matrix"<<endl;
                      return;
                    }}}// row=0
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

    }//no corrlation case


  //calculating cov:
  for(int i=0; i<NoOfObservables; i++)
    {
      for(int j=0;j<NoOfObservables; j++)
        {
          HL_cov(i,j)=(stat_error[i]*stat_error[j]+syst_error[i]*syst_error[j]) * HL_correlation(i,j);
        } // j
    }// i

  restricted=false;
  //cout<<HL_correlation<<endl;
  //cout<<HL_cov<<endl;
}
bool HL_nDimGaussian::Restrict(std::vector<std::string> names)
{
  size_restricted=  names.size();
  std::vector<int> indexes;
  for(int i=0 ; i<NoOfObservables ; i++)
    {
      std::string element=Observables[i];
      if(std::find(names.begin(), names.end(), element) != names.end()) {
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
          HL_cov_restricted(i,j)=HL_cov(i_index, j_index);
          HL_cov_inv_restricted(i,j)=HL_cov(i_index, j_index);
          HL_correlation_restricted(i,j)=HL_correlation(i_index, j_index);
        }
    }
  restricted=true;
  HL_Stats::InvertMatrix(HL_cov_restricted,HL_cov_inv_restricted);


  return true;
}


double HL_nDimGaussian::GetChi2(std::vector<double> theory)  //, double theory_err)
{
  if(!restricted) // if we don't resctric and use whole matrix
    {
      HL_cov_restricted=HL_cov;
      HL_correlation_restricted=HL_correlation;
      HL_cov_inv_restricted=HL_cov;
      central_restricted=central;

      //HL_Stats::inverse(HL_correlation_restricted, size_restricted);
      size_restricted=NoOfObservables;
      HL_Stats::InvertMatrix(HL_cov_restricted,HL_cov_inv_restricted);
      restricted=true;
    }
  // now calculating chi2
  double chi2=0;
  vector<double> diff;
  if(theory.size() !=  central_restricted.size())
    {
      std::cout<<"Error in HL_nDimGaussian::GetChi2, you had different dimensions in theory and experiment"<<std::endl;
      return -1e10;
    }
  for(int i=0;i<central_restricted.size(); i++) { diff.push_back(central_restricted[i] - theory[i] );}

  //cout<<HL_cov_restricted<<endl;
  //cout<<HL_cov_inv_restricted<<endl;
  for (int i=0; i < HL_cov_inv_restricted.size1(); ++i)
    {
      for (int j=0; j<HL_cov_inv_restricted.size2(); ++j)
        {
          chi2+= diff[i] * HL_cov_inv_restricted(i,j)*diff[j] ;
        }
    }

  return chi2;
}
double HL_nDimGaussian::GetLogLikelihood(std::vector<double> theory)
{

  double chi2=GetChi2(theory);

  return -0.5*chi2;
}
double HL_nDimGaussian::GetLikelihood(std::vector<double> theory)
{
  double log_likelihood=GetLogLikelihood(theory);
  return gsl_sf_exp(log_likelihood);
}




double HL_nDimGaussian::GetChi2(std::vector<double> theory , boost::numeric::ublas::matrix<double> theory_cov)
{
  if(theory_cov.size1() != theory.size() )
    {
      std::cout<<"Error in HL_nDimGaussian::GetChi2, you had different dimensions in theory and cov matrix"<<std::endl;
     }
  if(theory_cov.size2() != theory_cov.size1() )
    {
      std::cout<<"Error in HL_nDimGaussian::GetChi2, your theory cov matrix is not square!"<<std::endl ;
      
    }
  
  if(!restricted) // if we don't resctric and use whole matrix
    {
      HL_cov_restricted=HL_cov;
      HL_correlation_restricted=HL_correlation;
      HL_cov_inv_restricted=HL_cov;
      central_restricted=central;

      size_restricted=NoOfObservables;

      HL_cov_restricted+=theory_cov;
      cout<<"experimental cov1: "<<HL_cov_restricted<<endl;
      HL_Stats::InvertMatrix(HL_cov_restricted,HL_cov_inv_restricted);
      
      
      restricted=true;
    }
  else
      {
        cout<<"experimental cov2: "<<HL_cov_restricted<<endl;       
        HL_cov_restricted+=theory_cov;
        HL_Stats::InvertMatrix(HL_cov_restricted,HL_cov_inv_restricted);
      }

    // now calculating chi2
    double chi2=0;
    vector<double> diff;
    if(theory.size() !=  central_restricted.size())
      {
        std::cout<<"Error in HL_nDimGaussian::GetChi2, you had different dimensions in theory and experiment"<<std::endl  ;
        return -1e10;
      }
    for(int i=0;i<central_restricted.size(); i++) { diff.push_back(central_restricted[i] - theory[i] );}


    for (int i=0; i < HL_cov_inv_restricted.size1(); ++i)
      {
        for (int j=0; j<HL_cov_inv_restricted.size2(); ++j)
          {
            chi2+= diff[i] * HL_cov_inv_restricted(i,j)*diff[j] ;
          }
      }
    
    return chi2;
}
double HL_nDimGaussian::GetLogLikelihood(std::vector<double> theory, boost::numeric::ublas::matrix<double> theory_cov)
{
  double chi2=GetChi2(theory, theory_cov);
  return -0.5*chi2;
}
double HL_nDimGaussian::GetLikelihood(std::vector<double> theory, boost::numeric::ublas::matrix<double> theory_cov)
{
  double log_likelihood=GetLogLikelihood(theory,  theory_cov);
  return gsl_sf_exp(log_likelihood);
}
  
