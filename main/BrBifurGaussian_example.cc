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
  br->Read();
  cout<<"likehood: "<<br->GetLikelihood(0.16)<<endl;
  
  

  return 0;
}
