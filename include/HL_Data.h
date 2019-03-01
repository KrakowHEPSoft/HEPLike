//   HEPLike: High Energy Physics Likelihoods
//
//   Header for baseic class
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//////////////////////////////////////////////////
#ifndef HL_DATA_H
#define HL_DATA_H

#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>



#include "HL_Stats.h"

#include "yaml-cpp/yaml.h"


class HL_Data
{

 public:

  HL_Data();
  HL_Data(std::string);
  
  std::string HFile;
  void Read();
  
  
  
 protected:

  YAML::Node config;

  //HL_LIke info:
  std::string HL_FileName;
  std::string HL_Name;
  
  // bibligraphic data:
  std::string HL_DOI;
  std::string HL_BibCite;
  std::string HL_BibEntry;
  std::string HL_Arxiv;
  std::string HL_SubmissionYear;
  std::string HL_PublicationYear;
    
  // Collaborations numbers
  std::string HL_Collaborations;
  std::string HL_Collaboration_number;
  
  // addition data
  std::string HL_ExperimentalDataLumi; // 1fb for example
  std::string HL_ExperimentalDataYears;
  std::string HL_Decay;
  std::string HL_Kinematics;
  
  




  bool initialized;
  bool HL_debug_yaml;

  void read_standard();
  void set_debug_yaml(bool);

  
};


#endif
