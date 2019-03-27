#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "HL_Data.h"
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

  if(argc != 2 )
    {
      std::cout<<"Wrong number of arguments, please use in the following way:"<<std::endl;
      std::cout<<"./Test_YAML <YAML FILE>"<<std::endl;
      return -1; 
    }


  HL_Data *br = new HL_Data(argv[1]);
  br->set_debug_yaml(true);
  br->Read();

  return 0;
}
