//   HEPLike: High Energy Physics Likelihoods
//
//   Header importing ROOT headers or wrapper classes
//
//   author: Tomas Gonzalo
//////////////////////////////////////////////////
#ifndef HL_ROOT_H
#define HL_ROOT_H

// If using root, just include relevant headers
#ifdef USE_ROOT

  #include "TFile.h"
  #include "TTree.h"
  #include "TGraph.h"
  #include "TAxis.h"
  #include "TMath.h"
  #include "TH1D.h"
  #include "TH2D.h"
  #include "TH3D.h"
  #include "TBranch.h"
  #include "TRandom2.h"
  #include "TError.h"

  #include "Math/Minimizer.h"
  #include "Math/Factory.h"
  #include "Math/Functor.h"
  #include "Math/IFunction.h"

// If not, then just throw errors if types are used
#else


class TBranch
{
  public:

    TBranch()  { throw std::runtime_error("TBranch class is only available in ROOT, please enable ROOT to use"); }
    TBranch(const TBranch&)  { throw std::runtime_error("TBranch class is only available in ROOT, please enable ROOT to use"); }
};

class TTree
{
  public:

    TTree(const TTree&)  { throw std::runtime_error("TTree class is only available in ROOT, please enable ROOT to use"); }
    int GetEntries()  { throw std::runtime_error("TTree class is only available in ROOT, please enable ROOT to use"); }
    int GetEntry(long int) { throw std::runtime_error("TTree class is only available in ROOT, please enable ROOT to use"); }
    int SetBranchAddress(const char *, void *add, TBranch ** = 0) { throw std::runtime_error("TTree class is only available in ROOT, please enable ROOT to use"); }
};

class TFile
{
  public:

    TFile(std::string, std::string) { throw std::runtime_error("TFile class is only available in ROOT, please enable ROOT to use"); }
    TTree* Get(std::string) { throw std::runtime_error("TFile class is only available in ROOT, please enable ROOT to use"); }

};

class TGraph
{
  public:

    TGraph(const TGraph&)  { throw std::runtime_error("TGraph class is only available in ROOT, please enable ROOT to use"); }
};


#endif

#endif
