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
  HL_BifurGaussian *br = new HL_BifurGaussian("data/examples/test_bi.yaml");
  br->read();
  cout<<"likehood: "<<br->GetLikelihood(0.16)<<endl;
  
  
  /*
  HL_Limit *limit = new HL_Limit("data/HFLAV/2019/tau/HFLAV_2019_180.yaml");
  limit->read();
  double br=atof(argv[1]);
  cout<<"CLS: "<<limit->GetCLs(br)<<endl;
  cout<<"Chi2: "<<limit->GetChi2(br)<<endl;
  cout<<"likehood: "<<limit->GetLikelihood(br)<<endl;
  cout<<"Loglikehood: "<<limit->GetLogLikelihood(br)<<endl;
  */
  /*
  // profile likelihood test
    HL_ProfLikelihood *br = new HL_ProfLikelihood("data/LHCb/RD/RKstar_3fb/RKstar_lowq2.yaml");
    br->read();
  */
  /*
  HL_nDimGaussian *br = new HL_nDimGaussian("data/test_3dim.yaml");
  br->read();
  vector<string> a;//
  a.push_back("BR1");
  a.push_back("BR3");
  br->Restrict(a);
  //TH1D *hist;
  vector<double> theory;
  theory.push_back(0.1);
  theory.push_back(0.1);
  cout<<"Chi2: "<<br->GetChi2(theory)<<endl;
  */
  /*
  HL_nDimBifurGaussian *br = new HL_nDimBifurGaussian("data/test_3dimassym.yaml");
  br->read();
  vector<string> a;//
  a.push_back("BR1");
  a.push_back("BR3");
  br->Restrict(a);
  //TH1D *hist;
  vector<double> theory;
  theory.push_back(-0.1);
  theory.push_back(0.9);
  cout<<"Chi2: "<<br->GetChi2(theory)<<endl;
  */
  /*
  HL_nDimLikelihood *br = new HL_nDimLikelihood("data/LHCb/RD/Bs2mumu_5fb/b2mumu.yaml");
  br->read();  
  */

  return 0;
}
