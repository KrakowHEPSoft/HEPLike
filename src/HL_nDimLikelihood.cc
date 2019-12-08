//   HEPLike: High Energy Physics Likelihoods
//
//   Module to construck likelihoods for ndim likelihood
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
#include "HL_nDimLikelihood.h"


using namespace std;



void HL_nDimLikelihood::Read()
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

      
      path=path+HL_RootFile;
      HL_RootFile=path;
      cout<<"HEPLike testing the path: "<<HL_RootFile<<endl;
    }
  else
    {
      std::cout<<"You didn't profice a root file!!! HL_ProfLikelihood class is protesting!"<<std::endl;
    }
  if(config["TH2Path"])
    {
      HL_PATH=config["TH2Path"].as<std::string>();
      dim=2;
    }
  else if(config["TH3Path"])
    {
      HL_PATH=config["TH3Path"].as<std::string>();
      dim=3;
    }

  TFile *f= new TFile(HL_RootFile.c_str(), "READ");
  if(dim==2)
    {
      hist2D=dynamic_cast<TH2D*>(f->Get(HL_PATH.c_str()));
      n_binsX=hist2D->GetNbinsX();
      n_binsY=hist2D->GetNbinsY();
      n_binsZ=hist2D->GetNbinsZ();
      hist=hist2D;

    }
  else if(dim==3)
    {
      hist3D=dynamic_cast<TH3D*>(f->Get(HL_PATH.c_str()));

      n_binsX=hist3D->GetNbinsX();
      n_binsY=hist3D->GetNbinsY();
      n_binsZ=hist3D->GetNbinsZ();
      hist=hist3D;

    }
  xmin=hist->GetXaxis()->GetXmin();
  xmax=hist->GetXaxis()->GetXmax();
  ymin=hist->GetYaxis()->GetXmin();
  ymax=hist->GetYaxis()->GetXmax();
  zmin=0.;
  zmax=0.;
  if(dim==3)
    {
      zmin=hist->GetZaxis()->GetXmin();
      zmax=hist->GetZaxis()->GetXmax();
    }

  if(config["Observables"])
    {
      YAML::Node node  = config["Observables"];
      for(YAML::const_iterator it = node.begin(); it != node.end();  ++it )
        {
          Observables.push_back( ((*it)[0]).as<std::string>()  );
          central_mes_val.push_back( ((*it)[1]).as<double>()  );
        }


    }
  profiled=false;

  //cout<<central_mes_val[0]<<" "<<central_mes_val[1]<<endl;


}
double HL_nDimLikelihood::GetChi2(std::vector<double> theory)  //, double theory_err)
{
  double log_likelihood=GetLogLikelihood(theory);
  return -2.*log_likelihood;

}
double HL_nDimLikelihood::GetLogLikelihood(std::vector<double> theory)
{
  int bin;
  if(theory.size() ==2)
    {
      bin=hist2D->FindBin(theory[0], theory[1]);

    }
  double log_likelihood=hist2D->GetBinContent(bin);
  return -log_likelihood;
}
double HL_nDimLikelihood::GetLikelihood(std::vector<double> theory)
{
  double log_likelihood=GetLogLikelihood(theory);
  return gsl_sf_exp(log_likelihood);
}
double HL_nDimLikelihood::GetLogLikelihood(std::vector<double> theory , vector<double> theory_err)
{
  double log_likelihood=GetLogLikelihood(theory);
  double chi2=-2.*log_likelihood;
  
  for(unsigned i=0; i< theory.size(); i++)
    {
      chi2+=(theory[i]-central_mes_val[i])*(theory[i]-central_mes_val[i])/(  theory_err[i]*theory_err[i]);

    }
  
  return -0.5*chi2;
        
}



void HL_nDimLikelihood::Profile()
{
  //profiling over X:
  hist_profileX=new TH1D("profX", "profX", n_binsX,xmin, xmax);
  
  profiled=true;

  if(dim==2){
    for(int ix=1 ; ix < n_binsX ; ++ix)
      {
        double min=1.e10;

        for(int iy=1 ; iy < n_binsY ; ++iy)
          {
            if(hist->GetBinContent(ix,iy)<min) min=hist->GetBinContent(ix,iy);

          }//
        hist_profileX->SetBinContent(ix, min);
      }
  }
  else if(dim==3)
    {
      for(int ix=1 ; ix < n_binsX ; ++ix)
        {
          double min=1.e10;

          for(int iy=1 ; iy < n_binsY ; ++iy)
            {
              for(int iz=1 ; iz< n_binsZ ; ++iz)
                {
                  if(hist->GetBinContent(ix,iy,iz)<min) min=hist->GetBinContent(ix,iy,iz);
                }
            }
          hist_profileX->SetBinContent(ix, min);
        }
    }
  // profile Y:
  hist_profileY=new TH1D("profY", "profY", n_binsY,ymin, ymax);

  if(dim==2){
    for(int iy=1 ; iy < n_binsY ; ++iy)
      {
        double min=1.e10;

        for(int ix=1 ; ix < n_binsX ; ++ix)
          {
            if(hist->GetBinContent(ix,iy)<min) min=hist->GetBinContent(ix,iy);

          }//
        hist_profileY->SetBinContent(iy, min);
      }
  }
  else if(dim==3)
    {
      for(int iy=1 ; iy < n_binsY ; ++iy)
        {
          double min=1.e10;

          for(int ix=1 ; ix < n_binsX ; ++ix)
            {
              for(int iz=1 ; iz< n_binsZ ; ++iz)
                {
                  if(hist->GetBinContent(ix,iy,iz)<min) min=hist->GetBinContent(ix,iy,iz);
                }
            }
          hist_profileY->SetBinContent(iy, min);
        }
    }
  cout<<"?"<<endl;
  TFile *ftmp=new TFile("tmp.root", "RECREATE");
  hist_profileY->Write();
  hist_profileX->Write();
  if(dim==2) return;

  hist_profileZ=new TH1D("profZ", "profZ", n_binsZ,zmin, zmax);       
  for(int iz=1 ; iz < n_binsZ ; ++iz)
    {
      double min=1.e10;

      for(int iy=1 ; iy < n_binsY ; ++iy)
        {
          for(int ix=1 ; ix< n_binsX ; ++ix)
            {
              if(hist->GetBinContent(ix,iy,iz)<min) min=hist->GetBinContent(ix,iy,iz);
            }
        }
      hist_profileZ->SetBinContent(iz, min);
    }
  

}
double HL_nDimLikelihood::GetLogLikelihood_profile(  double theory, std::string X)
{

  if(X==Observables[0])
    {
      int bin= hist_profileX->FindBin(theory);
      double log_likelihood=hist_profileX->GetBinContent(bin);
      return -log_likelihood;
    }
  else if(X==Observables[1])
    {
      int  bin= hist_profileY->FindBin(theory);
      double log_likelihood=hist_profileY->GetBinContent(bin);
      return -log_likelihood;
    }
  else if(X==Observables[2])
    {
      int bin= hist_profileZ->FindBin(theory);
      double log_likelihood=hist_profileZ->GetBinContent(bin);
      return -log_likelihood;
    }
  else
    {
      std::cout<<"WRONG observable NAME!! in HL_nDimLikelihood profiling"<<std::endl;
      return 1.e20;
    }

}
double HL_nDimLikelihood::GetLikelihood_profile(double theory, std::string X)
{
  double log_likelihood=GetLogLikelihood_profile(theory, X);
  return gsl_sf_exp(log_likelihood);
}
double HL_nDimLikelihood::GetChi2_profile(double theory, std::string X)
{
  double log_likelihood=GetLogLikelihood_profile(theory, X);
  return -2.*log_likelihood;

}



