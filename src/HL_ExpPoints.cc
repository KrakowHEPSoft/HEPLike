//   HEPLike: High Energy Physics Likelihoods
//
//   Module to construct Experimental Data Points class 
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

#include "HL_ExpPoints.h"

using namespace std;




void HL_ExpPoints::Read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }

  read_standard();
  


  if( config["ROOTData"])  HL_RootFile=config["ROOTData"].as<std::string>();
  else
    {
      std::cout<<"You didn't provide a root file!!! HL_ExpPoints class is protesting!"<<std::endl;
    }
  
  if(config["TTreePath"]) HL_PATH=config["TTreePath"].as<std::string>();
  // now opening files
  f= new TFile(HL_RootFile.c_str(), "READ");
  HL_tree=dynamic_cast<TTree*>(f->Get(HL_PATH.c_str()));
  if(config["Weight"]) HL_weight=config["Weight"].as<std::string>();   
  else{
    std::cout<<"You didn't provide a weight name!!! HL_ExpPoints class is protesting!"<<std::endl;
  }
  if(config["Observables"])
    {
      YAML::Node node  = config["Observables"];         
      for(YAML::const_iterator it = node.begin(); it != node.end();  ++it )
        {
          HL_obs.push_back( ((*it)[0]).as<std::string>()  );
        }
    }// if Observables exist
  else{
    std::cout<<"You didn't provide a observable name!!! HL_ExpPoints class is protesting!"<<std::endl;
  }
  bool success=InitData();
  if(!success ) std::cout<<"HL_ExpPoints couldn't read data points!!! HL_ExpPoints class is protesting!"<<std::endl;
  


}
bool HL_ExpPoints::InitData()
{
  //initializing tree:
  int entries=HL_tree->GetEntries();
  int nVars=HL_obs.size();
  vector<TBranch*> branches(nVars);
  double vars[nVars];
  double weight_tmp;
  for(unsigned i =0; i < nVars ; ++i)
    {
      HL_Branches.push_back( new TBranch);
      //HL_tree->Branch(HL_obs[i].c_str(), &vars[i], (HL_obs[i]+"/D").c_str()); 
      HL_tree->SetBranchAddress( HL_obs[i].c_str(), &vars[i], &HL_Branches[i] );
      
    }
  HL_tree->SetBranchAddress( HL_weight.c_str() , &weight_tmp, &HL_weight_branch);    
  
  
  // storing data points in memory:
  cout<<entries<<endl;
  for(unsigned i=0; i < entries; ++i)
    {
      HL_tree->GetEntry(i);
      weights.push_back(weight_tmp);
      vector<double> tmp;
      for(int j=0; j < nVars ; ++j)
        {
          //cout<<vars[j]<<" "<<weights[i]<<endl;
          tmp.push_back(vars[j]);
        }
      points.push_back(tmp);
    }  
  return true;
}
void HL_ExpPoints::SetFun(  double FUN( vector<double> , vector<double> ) )
{
  fun=FUN;
}


double HL_ExpPoints::GetChi2(vector<double> theory)
{
  double log_likelihood=GetLogLikelihood(theory);
  double chi2=-2.*log_likelihood; 
  return chi2;
}

double HL_ExpPoints::GetLogLikelihood(vector<double> theory)
{
  double loglikelihood=0.;
  for(unsigned i=0; i <points.size(); ++i)
    {
      loglikelihood += weights[i] * log( (*fun )( theory, points[i] ));
    }
  return loglikelihood;
}
double HL_ExpPoints::GetLikelihood(vector<double> theory )
{
  double log_likelihood=GetLogLikelihood(theory);
  return gsl_sf_exp(log_likelihood);  
}




