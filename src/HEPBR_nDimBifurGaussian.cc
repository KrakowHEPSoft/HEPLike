//   HEPLike: High Energy Physics Likelihoods
//
//   Module to construck likelihoods for ndim bifurcated gaussian distribution
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
#include "HEPBR_nDimBifurGaussian.h"


using namespace std;



void HEPBR_nDimBifurGaussian::read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
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
      
      HEP_correlation= boost::numeric::ublas::matrix<double>(NoOfObservables,NoOfObservables);
      HEP_cov = boost::numeric::ublas::matrix<double>(NoOfObservables,NoOfObservables);
      HEP_cov_inv = boost::numeric::ublas::matrix<double>(NoOfObservables,NoOfObservables);
      for(YAML::const_iterator it = node.begin(); it != node.end();  ++it )
        {
          if(row==0)
            {
              for(int i=0;i<NoOfObservables; i++)
                {
                  if(Observables[i] != ((*it)[i]).as<std::string>()  )
                    {
                      cout<<"Error in HEPBR_nDimBifurGaussian, wrong ordering of names in correlation matrix"<<endl;
                      return;
                    }}}// row=0
          else
            {
              for(int i=0;i<NoOfObservables; i++)
                {
                  HEP_correlation(row-1,i)= ((*it)[i]).as<double>();
                }
            }
          
          row++;
        }
    }
  else // no correlation
    {
      HEP_correlation= boost::numeric::ublas::matrix<double>(NoOfObservables,NoOfObservables);
      HEP_cov = boost::numeric::ublas::matrix<double>(NoOfObservables,NoOfObservables);
      HEP_cov_inv = boost::numeric::ublas::matrix<double>(NoOfObservables,NoOfObservables);
      
      
      for(int i=0;i<NoOfObservables; i++)
        {
          for(int j=0;j<NoOfObservables; j++)
            {
              if(i==j)
                {
                  HEP_correlation(i,j) =1.;
                  HEP_cov(i,j) =1.;  
                  HEP_cov_inv(i,j)=1.;   
                }//diagonal
              else
                {
                  HEP_correlation(i,j) =0.;
                  HEP_cov(i,j) =0.;
                  HEP_cov_inv(i,j)=0.;
                }//ofdiagonal
            }// j
        }// i
      
    }//no corrlation case

  cout<<HEP_correlation<<endl;
  //cout<<HEP_cov<<endl;
}
bool HEPBR_nDimBifurGaussian::Restrict(std::vector<std::string> names)
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
  HEP_cov_restricted= boost::numeric::ublas::matrix<double>(size_restricted, size_restricted);
  HEP_correlation_restricted=  boost::numeric::ublas::matrix<double>(size_restricted, size_restricted);
  HEP_cov_inv_restricted=  boost::numeric::ublas::matrix<double>(size_restricted, size_restricted);

  for(int i=0; i< size_restricted; i++)
    {
      central_restricted.push_back(central [ indexes[i] ] );
      for(int j=0; j< size_restricted; j++)
        {
          int i_index=indexes[i];
          int j_index=indexes[j];
          HEP_correlation_restricted(i,j)=HEP_correlation(i_index, j_index); 
        }
    }
  cout<<HEP_correlation_restricted<<endl;
  restricted=true;
  
  return true;
}


double HEPBR_nDimBifurGaussian::GetChi2(std::vector<double> theory)  //, double theory_err)
{
  if(!restricted) // if we don't resctric and use whole matrix
    {
      HEP_correlation_restricted=HEP_correlation;
      //HEPStats::inverse(HEP_correlation_restricted, size_restricted);   
      size_restricted=NoOfObservables;
      HEPStats::InvertMatrix(HEP_cov_restricted,HEP_cov_inv_restricted);
      restricted=true; 
    }
  // now calculating chi2
  double chi2=0;
  vector<double> diff;
  if(theory.size() !=  central_restricted.size())
    {
      std::cout<<"Error in HEPBR_nDimBifurGaussian::GetChi2, you had different dimensions in theory and experiment"<<std::endl;
      return -1e10;
    }
  for(int i=0;i<central_restricted.size(); i++) { diff.push_back(central_restricted[i] - theory[i] );}


  HEP_cov_restricted=  boost::numeric::ublas::matrix<double>(size_restricted, size_restricted);
  HEP_cov_inv_restricted=  boost::numeric::ublas::matrix<double>(size_restricted, size_restricted);

  for(int i=0; i< size_restricted; i++)
    {
      for(int j=0; j< size_restricted; j++)
        {
          cout<<diff[i]<<" "<<diff[j]<<endl;
          if(diff[i] >= 0. && diff[j] >= 0.)      HEP_cov_restricted(i,j)=error_right[i]*error_right[j]*HEP_correlation_restricted(i,j);
          if(diff[i] >= 0. && diff[j] < 0.)      HEP_cov_restricted(i,j)=error_right[i]*error_left[j]*HEP_correlation_restricted(i,j);
          if(diff[i] < 0. && diff[j] >= 0.)      HEP_cov_restricted(i,j)=error_left[i]*error_right[j]*HEP_correlation_restricted(i,j);
          if(diff[i] < 0. && diff[j] < 0.)      HEP_cov_restricted(i,j)=error_left[i]*error_left[j]*HEP_correlation_restricted(i,j);
        }
    }
  cout<<HEP_cov_restricted<<endl;
  HEPStats::InvertMatrix(HEP_cov_restricted,HEP_cov_inv_restricted);
  cout<<HEP_cov_restricted<<endl;
  cout<<HEP_cov_inv_restricted<<endl;
  for (int i=0; i < HEP_cov_inv_restricted.size1(); ++i)
    {
      for (int j=0; j<HEP_cov_inv_restricted.size2(); ++j)
        {
          chi2+= diff[i] * HEP_cov_inv_restricted(i,j)*diff[j] ;
        }
    }
  return chi2;
}
double HEPBR_nDimBifurGaussian::GetLogLikelihood(std::vector<double> theory)
{

  double chi2=GetChi2(theory);
  
  return -0.5*chi2;
}
double HEPBR_nDimBifurGaussian::GetLikelihood(std::vector<double> theory)
{
  double log_likelihood=GetLogLikelihood(theory);
  return gsl_sf_exp(log_likelihood);  
}






