//   HEPLike: High Energy Physics Likelihoods
//
//   Module to read yaml files
//
//   author: Marcin Chrzaszcz
//////////////////////////////////////////////////


#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "HEPStats.h"
#include "HEPConstants.h"
#include "HEPData.h"



HEPData:: HEPData()
{
  initialized=false;
  HFile="";

}
HEPData:: HEPData(std::string s)
{

  initialized=true;
  HFile=s;
}



void HEPData::read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }
  read_standard();
}
    

  
void HEPData::read_standard()  
{
  config = YAML::LoadFile(HFile);

  HEPFileName=config["FileName"].as<std::string>();
  if(config["Name"] ) HEPName=config["Name"].as<std::string>();
  if(config["DOI"] ) HEPDOI=config["DOI"].as<std::string>();
  if(config["BibCite"]) HEPBibCite=config["BibCite"].as<std::string>();
  if(config["BibEntry"]) HEPBibEntry=config["BibEntry"].as<std::string>();
  if(config["arxiv"]) HEPArxiv=config["arxiv"].as<std::string>();
  if(config["Arxiv"]) HEPArxiv=config["Arxiv"].as<std::string>();
  if(config["SubmissionYear"]) HEPSubmissionYear =config["SubmissionYear"].as<std::string>();
  if(config["PublicationYear"]) HEPPublicationYear=config["PublicationYear"].as<std::string>();
  if(config["HFLAV"]) HEPHFLAV=config["HFLAV"].as<std::string>();
  if(config["Source"]) HEPSource= config["Source"].as<std::string>();
  if(config["Collaborations"]) HEPCollaborations=config["Collaborations"].as<std::string>();
  if(config["ExperimentalDataSet"]) HEPExperimentalDataSet=config["ExperimentalDataSet"].as<std::string>();
  if(config["ExperimentalDataYears"]) HEPExperimentalDataYears=config["ExperimentalDataYears"].as<std::string>();
  if(config["Decay"]) HEPDecay=config["Decay"].as<std::string>(); 
  if(config["Observables"])
    {
      YAML::Node node  = config["Observables"];
      for(YAML::const_iterator it = node.begin(); it != node.end();  ++it )
        {
          HEPObservables.push_back( (*it).as<std::string>()  ); 
        }
      
    }
  if(config["Kinematics"]) HEPKinematics=config["Kinematics"].as<std::string>(); 

  
}
