//   HEPLike: High Energy Physics Likelihoods
//
//   Module to construck likelihoods for ndim likelihood
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
#include "HL_nDimLikelihood.h"

using namespace std;

static bool ndim_debug = false;

HL_nDimLikelihood::~HL_nDimLikelihood()
{
  delete hist2D;
  // delete hist3D;
  delete hist_profileX;
  delete hist_profileY;
  // delete hist_profileZ;
  delete gmin;
  delete fun;
}

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
    if (ndim_debug) std::cout << "HL_nDimLikelihood is using text input" << std::endl;
    std::string filename = find_path(config["TextData"].as<std::string>());
    if (ndim_debug) std::cout << "Opening file " << filename << std::endl;
    std::ifstream in(filename.c_str());
    // TODO: This assumes it is always 2D, need to change this to make it work with 3D too
    dim=2;
    in >> n_binsX >> xmin >> xmax >> n_binsY >> ymin >> ymax;
    // Data has entries from 0 to nxbins, so a total of nxbins+1
    n_binsX++;
    n_binsY++;
    double **x; x = new double *[n_binsX];
    double **y; y = new double *[n_binsX];
    double **c; c = new double *[n_binsX];

    for(int i=0; i<n_binsX; i++)
    {
      x[i] = new double[n_binsY];
      y[i] = new double[n_binsY];
      c[i] = new double[n_binsY];
      for(int j=0; j<n_binsY; j++)
      {
        double bx,by,binc;
        in >> bx >> by >> binc;
        x[i][j] = bx;
        y[i][j] = by;
        c[i][j] = binc;
      }
    }
    in.close();

    hist2D = new HL_Interpolator2D(n_binsX, n_binsY, x, y, c);
    hist2D->SetLimits(xmin,xmax,ymin,ymax);
  }
  else if( config["ROOTData"])
  {
    #ifdef USE_ROOT
      HL_RootFile = find_path(config["ROOTData"].as<std::string>());
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
        TH2D *hist2D_tmp=dynamic_cast<TH2D*>(f->Get(HL_PATH.c_str()));
        hist2D_tmp->SetDirectory(0);
        hist2D = new HL_Interpolator2D(hist2D_tmp->Clone());
        hist2D_tmp->Delete();
        f->Close();
        delete f;

        n_binsX = hist2D->TH->GetNbinsX();
        n_binsY = hist2D->TH->GetNbinsY();

        xmin=hist->TH->GetXaxis()->GetXmin();
        xmax=hist->TH->GetXaxis()->GetXmax();
        ymin=hist->TH->GetYaxis()->GetXmin();
        ymax=hist->TH->GetYaxis()->GetXmax();

        hist2D->SetLimits(xmin,xmax,ymin,ymax);
      }
      else if(dim==3)
      {
        throw std::runtime_error("Dimensions greater than 2 are not implemented");
        // TODO: Implement Interpolator3D
        //TH3D *hist3D_tmp=dynamic_cast<TH3D*>(f->Get(HL_PATH.c_str()));
        //hist3D_tmp->SetDirectory(0) ;
        //hist3D = new HL_Interpolator3D(hist3D_tmp->Clone());
        //hist3D_tmp->Delete();

        //n_binsX=hist3D->TH->GetNbinsX();
        //n_binsY=hist3D->TH->GetNbinsY();
        //n_binsZ=hist3D->TH->GetNbinsZ();

        //xmin=hist->TH->GetXaxis()->GetXmin();
        //xmax=hist->TH->GetXaxis()->GetXmax();
        //ymin=hist->TH->GetYaxis()->GetXmin();
        //ymax=hist->TH->GetYaxis()->GetXmax();
        //zmin=hist->TH->GetZaxis()->GetXmin();
        //zmax=hist->TH->GetZaxis()->GetXmax();

        //hist3D->SetLimits(xmin,xmax,ymin,ymax,zmin,zmax);

      }
    #else
      throw std::runtime_error("Requested ROOT file but ROOT is not enabled. Please enable ROOT or provide data in text format.");
    #endif

  }
  else
  {
    throw std::runtime_error("You didn't provide neither a text file nor a ROOT file!!! HL_nDimLikelihood class is protesting!");
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

  if(dim == 2)
    fun = new HL_Function2D(hist2D);
  else if(dim == 3)
  {
    throw std::runtime_error("Dimensions greater than 2 are not implemented");
    //fun = new HL_Function3D(hist3D);
  }
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
  double delta = 1e-2;

  if(dim == 2)
  {
    double xlim = theory[0], ylim = theory[1];
    if(theory[0]>xmax) xlim = xmax*(1-delta);
    if(theory[0]<xmin) xlim = xmin*(1+delta);
    if(theory[1]>ymax) ylim = ymax*(1-delta);
    if(theory[1]<ymin) ylim = ymin*(1+delta);
    if(theory[0]>xmax or theory[0]<xmin or theory[1]>ymax or theory[1]<ymin)
    {
      double sigma[] = {delta*(theory[0]+xlim) , delta*(theory[1]+ylim)};
      return -1. * hist2D->Eval(xlim,ylim) - pow((theory[0] - xlim)/sigma[0],2) - pow((theory[1] - ylim)/sigma[1],2);
    }
    else
    {
      double log_likelihood=hist2D->Eval(theory[0], theory[1]);
      return (-1.)*log_likelihood;
    }
  }
  else
    throw std::runtime_error("Dimensions greater than 2 are not implemented");
}

double HL_nDimLikelihood::GetLogLikelihood(std::vector<double> theory, boost::numeric::ublas::matrix<double> theory_cov)
{
  if(theory_cov.size1() != theory.size() )
  {
    throw std::runtime_error("Error in HL_nDimLikelihood::GetLogLikelihood  you had different dimensions in theory and cov matrix");
  }
  if(theory_cov.size2() != theory_cov.size1() )
  {
    throw std::runtime_error("Error in HL_nDimLikelihood::GetLogLikelihood, your theory cov matrix is not square!");
  }

  if(dim == 2)
  {
    if(theory[0]>xmax) return loglikelihood_penalty;
    if(theory[0]<xmin) return loglikelihood_penalty;
    if(theory[1]>ymax) return loglikelihood_penalty;
    if(theory[1]<ymin) return loglikelihood_penalty;

    fun->SetTheory(theory,theory_cov);

    gmin = new HL_Minimizer("ConjugatePR", 2);
    gmin->SetMaxIterations(1000000);  // for GSL
    gmin->SetTolerance(0.000001);

    gmin->SetFunction(fun);

    double step[2] = {0.05*sqrt(theory_cov(0,0)), 0.05*sqrt(theory_cov(1,1)) };
    double variable[2] = { theory[0], theory[1]};

    gmin->SetVariable(0,variable[0]+0.5*sqrt(theory_cov(0,0)), step[0]);
    gmin->SetVariable(1,variable[1]+0.5*sqrt(theory_cov(1,1)), step[1]);

    //gmin->SetVariableInitialRange(0,theory[0]-5.*sqrt(theory_cov(0,0)), theory[0]+5.*sqrt(theory_cov(0,0)));
    //gmin->SetVariableInitialRange(1,theory[1]-5.*sqrt(theory_cov(1,1)), theory[1]+5.*sqrt(theory_cov(1,1)));

    //gmin->SetVariableLimits(0,theory[0]-5.*sqrt(theory_cov(0,0)), theory[0]+5.*sqrt(theory_cov(0,0)));
    //gmin->SetVariableLimits(1,theory[1]-5.*sqrt(theory_cov(1,1)), theory[1]+5.*sqrt(theory_cov(1,1)));


    gmin->Minimize();
    //const double *theory_nuisance = gmin->X();
    std::vector<double> theory_nuisance = gmin->X();

    if(theory_nuisance[0]>xmax) return loglikelihood_penalty;
    if(theory_nuisance[0]<xmin) return loglikelihood_penalty;
    if(theory_nuisance[1]>ymax) return loglikelihood_penalty;
    if(theory_nuisance[1]<ymin) return loglikelihood_penalty;

    double log_likelihood=hist2D->Eval(theory_nuisance[0], theory_nuisance[1]);

    return (-1.)*log_likelihood;
  }
  else
    throw std::runtime_error("Dimensions greater than 2 are not implemented");
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


void HL_nDimLikelihood::Profile(std::string mode)
{
  profiled=true;

  double newX[n_binsX];
  double newY[n_binsY];

  auto compare = [&](double val, double min, double max)
  {
    if(mode == "min") return val < min;
    else if(mode == "max") return val > max;
    else throw std::runtime_error("Unknown mode");
  };

  // Profile over X
  if(dim==2)
  {
    for(int ix=0 ; ix < n_binsX ; ++ix)
    {
      double min = 1.e10;
      double max = 0;

      for(int iy=0 ; iy < n_binsY ; ++iy)
      {
        double val = hist2D->z_data[ix][iy];
        if(compare(val, min, max))
        {
          min = val;
          max = val;
          newX[ix] = hist2D->x_data[ix][iy];
          newY[ix] = val;
        }
      }
    }
    hist_profileX = new HL_Interpolator1D(n_binsX, newX, newY);
  }
  else if(dim==3)
  {
    throw std::runtime_error("Dimensions greater than 2 are not implemented");
    // TODO: Implement 3D
    //for(int ix=1 ; ix < n_binsX ; ++ix)
    //{
    //  double min=1.e10;
    //  double max = 0;
    //
    //  for(int iy=1 ; iy < n_binsY ; ++iy)
    //  {
    //    for(int iz=1 ; iz< n_binsZ ; ++iz)
    //    {
    //      double val = hist3D->w_data[ix][iy][iz]);
    //      if(compare(val, min, max)
    //      {
    //        min = val;
    //        max = val;
    //        newX[ix] = hist->x_data[ix][iy][iz];
    //        newY[ix] = val;
    //      }
    //    }
    //  }
    //}
    //hist_profileX = new HL_Interpolator1D(n_binsX, newX, newY);
  }

  // Profile over Y
  if(dim==2)
  {
    for(int iy=0 ; iy < n_binsY ; ++iy)
    {
      double min = 1.e10;
      double max = 0;

      for(int ix=0 ; ix < n_binsX ; ++ix)
      {
        double val = hist2D->z_data[ix][iy];
        if(compare(val, min, max))
        {
          min = val;
          max = val;
          newX[iy] = hist2D->y_data[ix][iy];
          newY[iy] = val;
        }
      }
    }
    hist_profileY = new HL_Interpolator1D(n_binsY, newX, newY);
  }
  else if(dim==3)
  {
    throw std::runtime_error("Dimensions greater than 2 are not implemented");
    // TODO: Implement 3D
    //for(int iy=0 ; iy < n_binsY ; ++iy)
    //{
    //  double min=1.e10;
    //  double max = 0;
    //
    //  for(int ix=0 ; ix < n_binsX ; ++ix)
    //  {
    //    for(int iz=0 ; iz< n_binsZ ; ++iz)
    //    {
    //      double val = hist3D->w_data[ix][iy][iz];
    //      if(compare(val, min, max)
    //      {
    //        min = val;
    //        max = val;
    //        newX[iy] = hist3D->y_data[ix][iy][iz];
    //        newY[iy] = val;
    //      }
    //    }
    //  }
    //}
    //hist_profileY = new HL_Interpolator1D(n_binsY, newX, newY);
  }

  // Profile over Z
  if(dim ==3)
  {
    throw std::runtime_error("Dimensions greater than 2 are not implemented");
    // TODO: Implement 3D
    //for(int iz=0 ; iz < n_binsZ ; ++iz)
    //{
    //  double min=1.e10;
    //  double max = 0;
    //
    //  for(int ix=0 ; ix < n_binsX ; ++ix)
    //  {
    //    for(int iy=0 ; iy< n_binsY ; ++iy)
    //    {
    //      double val = hist3D->w_data[ix][iy][iz];
    //      if(compare(val, min, max)
    //      {
    //        min = val;
    //        max = val;
    //        newX[iz] = hist3D->z_data[ix][iy][iz];
    //        newY[iz] = val;
    //      }
    //    }
    //  }
    //}
    //hist_profileZ = new HL_Interpolator1D(n_binsZ, newX, newY);
  }

}

double HL_nDimLikelihood::GetLogLikelihood_profile(double theory, std::string X)
{
  if(X==Observables[0])
  {
    double log_likelihood = hist_profileX->Eval(theory);
    return -log_likelihood;
  }
  else if(X==Observables[1])
  {
    double log_likelihood = hist_profileY->Eval(theory);
    return -log_likelihood;
  }
  else if(X==Observables[2])
  {
    throw std::runtime_error("Dimensions greater than 2 are not implemented");
    // TODO: Implement 3D
    //double log_likelihood = hist_profileZ->Eval(theory);
    //return -log_likelihood;
  }
  else
    throw std::runtime_error("WRONG observable NAME in HL_nDimLikelihood profiling");
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



