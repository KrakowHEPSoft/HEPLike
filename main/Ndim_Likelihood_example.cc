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

#include "yaml-cpp/yaml.h"

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "TH1D.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"



using namespace std;



int main (int argc, char *argv[])
{
  TCanvas *c1 = new TCanvas("c1", "c1", 800,600);
  TCanvas *c2 = new TCanvas("c2", "c2", 800,600);  
  if(argc !=2)
    {
      return 1;

    }
  string pwd=argv[1];
  
  HL_nDimLikelihood *br = new HL_nDimLikelihood(pwd+"/data/examples/b2mumu.yaml");
  br->Read();

  cout<<"AA"<<endl;
  
  TH2D *hist=dynamic_cast<TH2D*>(br->GetHist());
  cout<<"AA"<<endl;
  TH2D *hist_post=dynamic_cast<TH2D*>(hist->Clone("AA"));
  cout<<"AA"<<endl;
  boost::numeric::ublas::matrix<double> theory_cov(2,2);
  theory_cov(0,0)=(6.e-9)*(6.e-9)*0.1*0.1;
  theory_cov(1,1)=(4.e-10)*(4.e-10)*0.1*0.1;  
  theory_cov(0,1)=0.;
  theory_cov(1,0)=0.;
  //vector<double> theory={6.e-9, 4.e-10};

    
  int Nx=hist_post->GetNbinsX();
  int Ny=hist_post->GetNbinsY(); 

  
  for(unsigned i=1; i<=Nx; ++i)
    {
      for(unsigned j=1; j<=Ny; ++j) 
        {
          
          double x = hist->GetXaxis()->GetBinCenter(i);
          double y = hist->GetYaxis()->GetBinCenter(j);  
          vector<double> theory={x,y};
          double LL=br->GetLogLikelihood(theory,theory_cov);
          hist_post->SetBinContent(i,j,-LL);

        }
    }
  c1->cd();
  hist->Draw("COLZ");
  c1->SaveAs("oryginal.pdf");
  c2->cd();
  hist_post->Draw("COLZ");
  c2->SaveAs("new.pdf");
  

  return 0;
}
