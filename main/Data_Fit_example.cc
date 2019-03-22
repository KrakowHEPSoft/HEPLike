#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>


#include "HL_Stats.h"
#include "HL_Constants.h"
#include "HL_Limit.h"
#include "HL_Gaussian.h"
#include "HL_BifurGaussian.h"
#include "HL_ProfLikelihood.h"
#include "HL_nDimGaussian.h"
#include "HL_nDimBifurGaussian.h"
#include "HL_nDimLikelihood.h"
#include "HL_ExpPoints.h"

#include "yaml-cpp/yaml.h"

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"

using namespace std;


double FUN(vector<double> in, vector<double> point)
{
  double mean=in[0];
  double sigma=in[1];
  double x=point[0];
  double gauss= 1./(sqrt(2. * HL_Const::pi * sigma*sigma)) *pow(HL_Const::e, -(x-mean)*(x-mean)/(2.*sigma*sigma) );
  return gauss;
}


int main (int argc, char *argv[])
{
  TCanvas *c1 = new TCanvas("c1", "c1", 800,600);
  
  // profile likelihood test
  HL_ExpPoints *br = new HL_ExpPoints("/storage/github/HEPLike/data_toy/toy/dummy.yaml");
  br->Read();
  br->SetFun(FUN);

  double sigma=2.;
  double mean=-1.;

  vector<double> MEAN;
  vector<double> LL;
  double min=1.e10;
  while (mean<2.)
    {
      vector<double> point;
      point.push_back(mean);
      point.push_back(sigma);
     
      cout<<"Loglikelihood= "<<br->GetLogLikelihood(point)<<" "<<mean<<endl;
      MEAN.push_back(mean );
      double tmp=(-1.)*br->GetLogLikelihood(point);
      if(tmp<min) min=tmp;
      LL.push_back(tmp);
      
      
    mean+=0.05;
  }
  double LLA[MEAN.size()];
  double MEANA[MEAN.size()];

  for(int i=0; i< MEAN.size() ; i++)
    {
      LLA[i]=LL[i]-min;
      MEANA[i]=MEAN[i];

    }
  TGraph* gr = new TGraph(MEAN.size(), MEANA, LLA);
  gr->Draw("AC*");
  c1->SaveAs("GaussFit.pdf");
  c1->SaveAs("GaussFit.root");


  return 1;
}
