//   HEPLike: High Energy Physics Likelihoods
//
//   Module to construck likelihoods for ndim likelihood
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
#include "HEPBR_nDimLikelihood.h"


using namespace std;



void HEPBR_nDimLikelihood::read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }
  read_standard();
  if( config["ROOTData"])  HEPRootFile=config["ROOTData"].as<std::string>();
  else
    {
      std::cout<<"You didn't profice a root file!!! HEPBR_ProfLikelihood class is protesting!"<<std::endl;
    }
  if(config["TH2Path"])
    {
      HEPPATH=config["TH2Path"].as<std::string>();
      dim=2;
    }
  else if(config["TH3Path"])
    {
      HEPPATH=config["TH3Path"].as<std::string>(); 
      dim=3;
    }
      
  TFile *f= new TFile(HEPRootFile.c_str(), "READ");
  if(dim==2) hist2D=dynamic_cast<TH2D*>(f->Get(HEPPATH.c_str()));
  else if(dim==3)  hist3D=dynamic_cast<TH3D*>(f->Get(HEPPATH.c_str()));
  if(config["Observables"])
    {
      YAML::Node node  = config["Observables"];
      for(YAML::const_iterator it = node.begin(); it != node.end();  ++it )
        {
          Observables.push_back( ((*it)[0]).as<std::string>()  );
          central_mes_val.push_back( ((*it)[1]).as<double>()  );
        }
      
      
    }
  cout<<central_mes_val[0]<<" "<<central_mes_val[1]<<endl;

  
}
double HEPBR_nDimLikelihood::GetChi2(std::vector<double> theory)  //, double theory_err)
{
  double log_likelihood=GetLogLikelihood(theory);
  return -2.*log_likelihood;

}
double HEPBR_nDimLikelihood::GetLogLikelihood(std::vector<double> theory)
{
  int bin;
  if(theory.size() ==2)
    {
      bin=hist2D->FindBin(theory[0], theory[1]);
      
    }
  double log_likelihood=hist2D->GetBinContent(bin);
  return -log_likelihood;

  
}
double HEPBR_nDimLikelihood::GetLikelihood(std::vector<double> theory)
{
  double log_likelihood=GetLogLikelihood(theory);
  return gsl_sf_exp(log_likelihood);  
}





