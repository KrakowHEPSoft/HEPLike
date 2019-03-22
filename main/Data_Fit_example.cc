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



int main (int argc, char *argv[])
{
  TCanvas *c1 = new TCanvas("c1", "c1", 800,600);
  
  // profile likelihood test
  HL_ExpPoints *br = new HL_ExpPoints("/storage/github/HEPLike/data_toy/toy/dummy.yaml");
  br->Read();
  


  return 1;
}
