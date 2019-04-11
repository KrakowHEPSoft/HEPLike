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

using namespace std;



int main (int argc, char *argv[])
{

  HL_nDimBifurGaussian *br = new HL_nDimBifurGaussian("data/examples/test_3dimassym.yaml");
  br->Read();
  vector<string> a;//
  a.push_back("BR1");
  a.push_back("BR3");
  br->Restrict(a);
  //TH1D *hist;
  vector<double> theory;
  theory.push_back(-0.1);
  theory.push_back(0.9);
  cout<<"Chi2: "<<br->GetChi2(theory)<<endl;


  return 0;
}
