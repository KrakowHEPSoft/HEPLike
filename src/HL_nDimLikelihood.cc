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

  // Note: amended to allow text file input for 2D arrays

  if(! initialized)
    {
      std::cout << "HL_nDimLikelihood Warning, TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }
  read_standard();
  loglikelihood_penalty=-1.e6;

  //Always use text input if it is available (avoids potential mem issues)

  if( config["TextData"])
    {
    std::cout << "HL_nDimLikelihood is using text input" << std::endl;
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

	std::string filename=path+"/"+config["TextData"].as<std::string>();
	std::cout << "Opening file " << filename << std::endl;
	std::ifstream in(filename.c_str());
	in >> n_binsX >> xmin >> xmax >> n_binsY >> ymin >> ymax;
	hist2D = new TH2D("HLhist","HLhist",n_binsX,xmin,xmax,n_binsY,ymin,ymax);

	for(int x=0;x<n_binsX+1;x++){
	  for(int y=0;y<n_binsY+1;y++){
	    double bx,by,binc;
	    in >> bx >> by>> binc;
	    hist2D->SetBinContent(bx,by,binc);
	  }
	}
	in.close();
	hist=hist2D;

    }
  else {

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

	std::cout<<"You didn't provide a root file for a text file!!! HL_nDimLikelihood class is protesting!"<<std::endl;
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
    cout<<"here"<<endl;
    TFile *f= new TFile(HL_RootFile.c_str(), "READ");
    if(dim==2)
      {
	TH2D *hist2D_tmp=dynamic_cast<TH2D*>(f->Get(HL_PATH.c_str()));
	hist2D=dynamic_cast<TH2D*>(hist2D_tmp->Clone());
	hist2D->SetDirectory(0);
	hist2D_tmp->Delete();

	//delete hist2D_tmp;

	n_binsX=hist2D->GetNbinsX();
	n_binsY=hist2D->GetNbinsY();
	n_binsZ=hist2D->GetNbinsZ();
	hist=hist2D;

      }
    else if(dim==3)
      {
	TH3D *hist3D_tmp=dynamic_cast<TH3D*>(f->Get(HL_PATH.c_str()));
	hist3D=dynamic_cast<TH3D*>(hist3D_tmp->Clone());
	hist3D->SetDirectory(0) ;
	hist3D_tmp->Delete();
	//delete hist3D_tmp;

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
    f->Close();
    delete f;
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

  //profiled=false;
  /*gmin=ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "ConjugatePR");
  //gmin=ROOT::Math::Factory::CreateMinimizer("Minuit2");
  gmin->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
  gmin->SetMaxIterations(1000000);  // for GSL
  gmin->SetTolerance(0.000001);
  gmin->SetPrintLevel(3);*/

  fun=MyFunction2D();
  fun.SetLikelihood(hist2D);

  //cout<<central_mes_val[0]<<" "<<central_mes_val[1]<<endl;

  cout<<"Finished"<<endl;
}
double HL_nDimLikelihood::GetChi2(std::vector<double> theory)  //, double theory_err)
{
  double log_likelihood=GetLogLikelihood(theory);
  return -2.*log_likelihood;

}
double HL_nDimLikelihood::GetChi2(std::vector<double> theory,  boost::numeric::ublas::matrix<double> theory_cov)
{
  double log_likelihood=GetLogLikelihood(theory);
  return -2.*log_likelihood;

}

double HL_nDimLikelihood::GetLogLikelihood(std::vector<double> theory)
{
  int bin;
  double delta = 1e-2;
  double xlim = theory[0], ylim = theory[1];
  if(theory[0]>xmax) xlim = xmax*(1-delta);
  if(theory[0]<xmin) xlim = xmin*(1+delta);
  if(theory[1]>ymax) ylim = ymax*(1-delta);
  if(theory[1]<ymin) ylim = ymin*(1+delta);
  /*
  if(theory.size() ==2)
    {
      bin=hist2D->FindBin(theory[0], theory[1]);

    }
    double log_likelihood=hist2D->GetBinContent(bin);
  */
  if(theory[0]>xmax or theory[0]<xmin or theory[1]>ymax or theory[1]<ymin)
  {
    double sigma[] = {delta*(theory[0]+xlim) , delta*(theory[1]+ylim)};
    return -1. * hist2D->Interpolate(xlim,ylim) - pow((theory[0] - xlim)/sigma[0],2) - pow((theory[1] - ylim)/sigma[1],2);
  }
  else
  {
    double log_likelihood=hist2D->Interpolate (theory[0], theory[1]);
    return (-1.)*log_likelihood;
  }
}

double HL_nDimLikelihood::GetLogLikelihood(std::vector<double> theory, boost::numeric::ublas::matrix<double> theory_cov)
{
  if(theory_cov.size1() != theory.size() )
    {
      std::cout<<"Error in HL_nDimLikelihood::GetLogLikelihood  you had different dimensions in theory and cov matrix"<<std::endl;
    }
  if(theory_cov.size2() != theory_cov.size1() )
    {
      std::cout<<"Error in HL_nDimLikelihood::GetLogLikelihood, your theory cov matrix is not square!"<<std::endl ;

    }

  int bin;
  if(theory[0]>xmax) return loglikelihood_penalty;
  if(theory[0]<xmin) return loglikelihood_penalty;
  if(theory[1]>ymax) return loglikelihood_penalty;
  if(theory[1]<ymin) return loglikelihood_penalty;

  fun.SetTheory(theory,theory_cov);
  ROOT::Math::Functor  f1(fun, 2);


  ROOT::Math::Minimizer* gmin=ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "ConjugatePR");
  //gmin=ROOT::Math::Factory::CreateMinimizer("Minuit2");
  gmin->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
  gmin->SetMaxIterations(1000000);  // for GSL
  gmin->SetTolerance(0.000001);
  gmin->SetPrintLevel(3);
  gmin->SetFunction(f1);

  double step[2] = {0.05*sqrt(theory_cov(0,0)), 0.05*sqrt(theory_cov(1,1)) };
  double variable[2] = { theory[0], theory[1]};


  gmin->SetVariable(0,"x",variable[0]+0.5*sqrt(theory_cov(0,0)), step[0]);
  gmin->SetVariable(1,"y",variable[1]+0.5*sqrt(theory_cov(1,1)), step[1]);

  gmin->SetVariableInitialRange(0,theory[0]-5.*sqrt(theory_cov(0,0)), theory[0]+5.*sqrt(theory_cov(0,0)));
  gmin->SetVariableInitialRange(1,theory[1]-5.*sqrt(theory_cov(1,1)), theory[1]+5.*sqrt(theory_cov(1,1)));

  gmin->SetVariableLimits(0,theory[0]-5.*sqrt(theory_cov(0,0)), theory[0]+5.*sqrt(theory_cov(0,0)));
  gmin->SetVariableLimits(1,theory[1]-5.*sqrt(theory_cov(1,1)), theory[1]+5.*sqrt(theory_cov(1,1)));


  gmin->Minimize();
  const double *theory_nuisance = gmin->X();

  if(theory_nuisance[0]>xmax) return loglikelihood_penalty;
  if(theory_nuisance[0]<xmin) return loglikelihood_penalty;
  if(theory_nuisance[1]>ymax) return loglikelihood_penalty;
  if(theory_nuisance[1]<ymin) return loglikelihood_penalty;
  /*
  bin=hist2D->FindBin(theory_nuisance[0], theory_nuisance[1]);

  double log_likelihood=hist2D->GetBinContent(bin);
  */
  double log_likelihood=hist2D->Interpolate(theory_nuisance[0], theory_nuisance[1]);


  return (-1.)*log_likelihood;
}


double HL_nDimLikelihood::GetLikelihood(std::vector<double> theory)
{
  double log_likelihood=GetLogLikelihood(theory);
  return gsl_sf_exp(log_likelihood);
}

double HL_nDimLikelihood::GetLikelihood(std::vector<double> theory, boost::numeric::ublas::matrix<double> theory_cov)
{
  double log_likelihood=GetLogLikelihood(theory,theory_cov);
  return gsl_sf_exp(log_likelihood);
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

  TFile *ftmp=new TFile("tmp.root", "RECREATE");
  hist_profileY->Write();
  hist_profileX->Write();
  //hist_profileY->Delete();
  //hist_profileX->Delete();


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
  hist_profileZ->Write();
  //  hist_profileZ->Delete();
  ftmp->Close();
  ftmp->Delete();
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



