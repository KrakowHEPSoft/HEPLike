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

  YAML::Node config;

  //HEPLIke info:
  std::string HEPFileName;
  std::string HEPName;
  
  // bibligraphic data:
  std::string HEPDOI;
  std::string HEPBibCite;
  std::string HEPBibEntry;
  std::string HEPArxiv;
  std::string HEPSubmissionYear;
  std::string HEPPublicationYear;
  
  
  
  // Collaborations numbers
  std::string HEPHFLAV;
  std::string HEPSource;
  std::string HEPCollaborations;

  
  // addition data
  std::string HEPExperimentalDataSet; // 1fb for example
  std::string HEPExperimentalDataYears;
  std::string HEPDecay;
  std::vector<std::string> HEPObservables;
  std::string HEPKinematics;
  




  bool initialized;

  void read_standard();

  
};


#endif
