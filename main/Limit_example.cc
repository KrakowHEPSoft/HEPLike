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
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"



using namespace std;



int main (int argc, char *argv[])
{
  HL_Limit *limit = new HL_Limit("data/examples/HFLAV_2019_180.yaml");
  limit->Read();
  double br=5.e-10;
  vector<double> BR;
  vector<double> like;
  vector<double> like_fake;

  HL_Limit *limit2 = new HL_Limit("data/examples/HFLAV_2019_180_dummy.yaml");
  limit2->Read();
  
  
  double max1=-1.e10;
  double max2=-1.e10;
  

  while(br<2e-8)
    {
      BR.push_back(br);
      double tmp_like_p=limit->GetLogLikelihood(br);
      like.push_back(tmp_like_p);
      double tmp_like=limit2->GetLogLikelihood(br); 
      like_fake.push_back(tmp_like);

      /*
      double error= fabs(9.9e-9)/1.64; 
      double tmp_ like=HL_Stats::gaussian_upper_limit(br, 0. ,0.,  error, false);
      cout<<br<<"  "<<   tmp_like_p<<"  "<<tmp_like<<endl; 
      like_fake.push_back(tmp_like);
      */
      if(tmp_like_p > max1) max1=tmp_like_p;
      if(tmp_like > max2) max2=tmp_like;

         
      br+=1e-10;
    }
  int bins=BR.size();
  double BRA[bins];
  double LL[bins];
  double LLF[bins];
  for(int i=0; i <bins; ++i)
    {
      BRA[i]=BR[i];
      LL[i]=like[i]-max1;
      LLF[i]=like_fake[i]-max2;

    }

  TGraph* gr = new TGraph(BR.size(), BRA, LL);
  TGraph* gr1 = new TGraph(BR.size(), BRA, LLF); 
  //############################################3
  TCanvas *c1= new TCanvas("c1", "c1", 800,600);

  gr->SetLineColor(kRed-3);
  gr1->SetLineColor(kGreen-3);
  gr->SetLineWidth(3);
  gr1->SetLineWidth(3);       

  gr->GetXaxis()->SetTitle("Br(#tau -> #mu #mu e)");
  gr1->GetXaxis()->SetTitle("Br(#tau -> #mu #mu e)");
  gr->GetYaxis()->SetTitle("#Delta log-likelihood");
  gr1->GetYaxis()->SetTitle("#Delta log-likelihood");    



  
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr);
  mg->Add(gr1); 

  mg->Draw("AL");   
  
  mg->GetXaxis()->SetTitle("Br(#tau -> #mu e e)");
  mg->GetYaxis()->SetTitle("#Delta log-likelihood");
  
  mg->Draw("AL");
  

  TLegend *leg = new TLegend(0.55,0.7,0.9,0.9);
  leg->AddEntry(gr, "Proper Likelihood", "l");
  leg->AddEntry(gr1, "Aprox. Proper Likelihood", "l");    

  leg->Draw();


  c1->SaveAs("dupa.png");

  
  return 0;
}
