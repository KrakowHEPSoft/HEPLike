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
#include "TGraph.h"
#include "TCanvas.h"

using namespace std;



int main (int argc, char *argv[])
{
  TCanvas *c1 = new TCanvas("c1", "c1", 800,600);
  
  // profile likelihood test
  HL_ProfLikelihood *br = new HL_ProfLikelihood("data/LHCb/RD/RKstar/CERN-EP-2017-100_q2_1.1_6.yaml");
  br->Read();
  
  // now let's see
  vector<double> BR;
  vector<double> LL;
  vector<double> LL2;
  
  double ibr=0.3;
  double dbr=0.005;
  double min1=1e10;
  double min2=1e10;
  while (ibr<1.3) {

    BR.push_back(ibr);
    double iLL=(-1.)*br->GetLogLikelihood(ibr, -2.);
    double iLL2=(-1.)*br->GetLogLikelihood(ibr,  0.1);// adding 10% error.
    LL.push_back(iLL);
    LL2.push_back(iLL2);
    
    if(iLL<min1) min1=iLL;
    if(iLL2<min2) min2=iLL2;
        
    ibr+=dbr;
  }
  int size=LL.size();
  double BR2[size];
  double LLA[size];
  double LL2A[size];

  
  for(int i=0; i< size ; i++)
    {
      BR2[i]=BR[i];
      LLA[i]=LL[i] -min1 ;
      LL2A[i]=LL2[i] - min2;;
      
    }
  TGraph* gr = new TGraph(size,BR2,LLA);
  TGraph* gr2 = new TGraph(size,BR2,LL2A);
  
  gr->Draw("AC*");
  c1->SaveAs("ProfLikelihood_example.pdf");
  gr2->SetLineColor(kBlue);
  gr2->Draw("AC* SAME");       
  c1->SaveAs("ProfLikelihood_example2.pdf");
  
  return 1;
}
