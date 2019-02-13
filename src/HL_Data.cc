//   HEPLike: High Energy Physics Likelihoods
//
//   Module to read yaml files
//
//   author: Jihyun Bhom,  Marcin Chrzaszcz
//////////////////////////////////////////////////


#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "HEPStats.h"
#include "HEPConstants.h"
#include "HL_Data.h"



HL_Data:: HL_Data()
{
  initialized=false;
  HFile="";

}
HL_Data:: HL_Data(std::string s)
{

  initialized=true;
  HFile=s;
}



void HL_Data::read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }
  read_standard();
}
    

  
void HL_Data::read_standard()  
{
  config = YAML::LoadFile(HFile);

  HL_FileName=config["FileName"].as<std::string>();
  if(config["Name"] ) HL_Name=config["Name"].as<std::string>();
  if(config["DOI"] ) HL_DOI=config["DOI"].as<std::string>();
  if(config["BibCite"]) HL_BibCite=config["BibCite"].as<std::string>();
  if(config["BibEntry"]) HL_BibEntry=config["BibEntry"].as<std::string>();
  if(config["arxiv"]) HL_Arxiv=config["arxiv"].as<std::string>();
  if(config["Arxiv"]) HL_Arxiv=config["Arxiv"].as<std::string>();
  if(config["SubmissionYear"]) HL_SubmissionYear =config["SubmissionYear"].as<std::string>();
  if(config["PublicationYear"]) HL_PublicationYear=config["PublicationYear"].as<std::string>();
  if(config["HFLAV"]) HL_HFLAV=config["HFLAV"].as<std::string>();
  if(config["Source"]) HL_Source= config["Source"].as<std::string>();
  if(config["Collaborations"]) HL_Collaborations=config["Collaborations"].as<std::string>();
  if(config["ExperimentalDataSet"]) HL_ExperimentalDataSet=config["ExperimentalDataSet"].as<std::string>();
  if(config["ExperimentalDataYears"]) HL_ExperimentalDataYears=config["ExperimentalDataYears"].as<std::string>();
  if(config["Decay"]) HL_Decay=config["Decay"].as<std::string>(); 
  /*
  if(config["Observables"])
    {
      YAML::Node node  = config["Observables"];
      for(YAML::const_iterator it = node.begin(); it != node.end();  ++it )
        {
          HL_Observables.push_back( (*it).as<std::string>()  ); 
        }
      
    }
  */
  if(config["Kinematics"]) HL_Kinematics=config["Kinematics"].as<std::string>(); 

  
}