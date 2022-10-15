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
  #include "TBranch.h"
  #include "TRandom2.h"
  #include "TError.h"

  #include "Math/Minimizer.h"
  #include "Math/Factory.h"
  #include "Math/Functor.h"

// If not, then wrap usual methods in root classes
#else

// TODO: Expand these types to do sometthing

class TBranch
{

};

class TTree : std::vector<TBranch*>
{
  public:

    // Number of entries
    int GetEntries();

    // Get one entry
    int GetEntry(long int);

    // Get address of branch
    int SetBranchAddress(const char *, void *add, TBranch ** = 0);

};

class TFile
{
  public:

    // Constructor
    TFile(std::string, std::string);

    // Getter
    TTree* Get(std::string);

};

class TGraph
{

};

namespace ROOT
{

}

#endif

#endif
