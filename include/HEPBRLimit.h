//   HEPLike: High Energy Physics Likelihoods
//
//   Header for limit class
//
//   author: Marcin Chrzaszcz
//////////////////////////////////////////////////
#ifndef HEPBRLIMIT_H
#define HEPBRLIMIT_H

//C++ headers
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>


//HEPLike headers
#include "HEPStats.h"
#include "HEPData.h"

#include "yaml-cpp/yaml.h"


class HEPBRLimit: public HEPData
{

 public:

  explicit HEPBRLimit() :  HEPData() {};
  explicit HEPBRLimit(std::string s) :  HEPData(s) { };
  
  //HEPBRLimit();
  //HEPBRLimit(std::string);
  
  //std::string HFile;
  void read();
  
  

 protected:
  std::vector<double> CLs;
  std::vector<double> BR;
  

 /*
  std::string HEPDOI;
  std::string HEPBibCite;
  std::string HEPBibEntry;
  std::string HEPFileName;
  std::string HEPHFAG;
  std::string HEPSource;
  std::string HEPYear;
  std::string HEPName;

      
  bool initialized;
  */


  
};


#endif
