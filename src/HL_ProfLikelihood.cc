//   HEPLike: High Energy Physics Likelihoods
//
//   Module to construck PROFLIKELIHOOD
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

#include "HL_ProfLikelihood.h"

using namespace std;

void HL_ProfLikelihood::Read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }

  read_standard();

  if( config["ROOTData"])
    {
      HL_RootFile=config["ROOTData"].as<std::string>();
      // the HL_RootFile is something like data/RD/... but we need the absolute path:
      int pos=HFile.find("/data/");
      //string path=HFile.substr (0, pos);
      if(pos<0) std::cout<<"Error in HL_nDimLikelihood, didn't find 'data'"<<std::endl;
      while(true)
        {
          int pos_new=HFile.find("/data/", pos+1);
          if(pos_new>0)
            {
              pos=pos_new;
            }
          else break;
        }
      string path=HFile.substr (0, pos);


      path=path+"/"+HL_RootFile;
      HL_RootFile=path;
    }
  else
    {
      std::cout<<"You didn't profice a root file!!! HL_ProfLikelihood class is protesting!"<<std::endl;
    }

  if(config["TGraphPath"]) HL_PATH=config["TGraphPath"].as<std::string>();

  // now opening files
  f= new TFile(HL_RootFile.c_str(), "READ");
  likelihood=dynamic_cast<TGraph*>(f->Get(HL_PATH.c_str()));
  likelihood->Print();



  xmin=likelihood->GetXaxis()->GetXmin () ;
  xmax=likelihood->GetXaxis()->GetXmax () ;
  std::cout<<xmin<<" "<<xmax<<std::endl;
  if(config["Observables"] )
    {
      YAML::Node node  = config["Observables"];
      ObsName=node[0][0].as<std::string>();
      central_mes_val=node[0][1].as<double>();

    }

}
double HL_ProfLikelihood::GetChi2(double theory, double theory_err=-1.)
{
  double log_likelihood=GetLogLikelihood(theory,theory_err);
  double chi2=-2.*log_likelihood;
  return chi2;
}

double HL_ProfLikelihood::GetLogLikelihood(double theory, double theory_error=-1.)
{
  if(theory < xmin || theory > xmax) return -1.e10;
  double loglikelihood=(-1.)*likelihood->Eval(theory,0, "S" );
  std::cout<<"Likelihood: "<<loglikelihood<<" for "<<theory<<std::endl;
  if(theory_error<0.){
      return loglikelihood;
  }
  // we inflate the likelihood by theory error square:
  loglikelihood*=(1.+theory_error)*(1.+theory_error);


  return loglikelihood;

}
double HL_ProfLikelihood::GetLikelihood(double theory, double theory_err=-1.)
{
  double log_likelihood=GetLogLikelihood(theory,theory_err);
  return gsl_sf_exp(log_likelihood);
}
