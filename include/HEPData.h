//   HEPLike: High Energy Physics Likelihoods
//
//   Header for limit class
//
//   author: Marcin Chrzaszcz
//////////////////////////////////////////////////
#ifndef HEPDATA_H
#define HEPDATA_H

#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>



#include "HEPStats.h"

#include "yaml-cpp/yaml.h"


class HEPData
{

 public:

  HEPData();
  HEPData(std::string);
  
  std::string HFile;
  void read();
  
  
  
 protected:

  std::string HEPDOI;
  std::string HEPBibCite;
  std::string HEPBibEntry;
  std::string HEPFileName;
  std::string HEPHFAG;
  std::string HEPSource;
  std::string HEPYear;
  std::string HEPName;

      
  bool initialized;
  
};


#endif
