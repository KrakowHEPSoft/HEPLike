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

#include "HL_Stats.h"
#include "HL_Constants.h"
#include "HL_Data.h"



HL_Data:: HL_Data()
{
  initialized=false;
  HL_debug_yaml=false;
  HFile="";

}
HL_Data:: HL_Data(std::string s)
{

  initialized=true;
  HL_debug_yaml=false;
  HFile=s;
}



void HL_Data::Read()
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
  else if (HL_debug_yaml) std::cout<<"HEPLike Warning, The <Name> is missing in the yaml file"<<std::endl; 

  if(config["DOI"] ) HL_DOI=config["DOI"].as<std::string>();
  else if (HL_debug_yaml) std::cout<<"HEPLike Warning, The <DOI> is missing in the yaml file"<<std::endl; 

  if(config["BibCite"]) HL_BibCite=config["BibCite"].as<std::string>();
  else if (HL_debug_yaml) std::cout<<"HEPLike Warning, The <BibCite> is missing in the yaml file"<<std::endl; 

  if(config["BibEntry"]) HL_BibEntry=config["BibEntry"].as<std::string>();
  else if (HL_debug_yaml) std::cout<<"HEPLike Warning, The <BibEntry> is missing in the yaml file"<<std::endl; 

  if(config["Arxiv"]) HL_Arxiv=config["Arxiv"].as<std::string>();
  else if (HL_debug_yaml) std::cout<<"HEPLike Warning, The <Arxiv> is missing in the yaml file"<<std::endl; 

  if(config["SubmissionYear"]) HL_SubmissionYear =config["SubmissionYear"].as<std::string>();
  else if (HL_debug_yaml) std::cout<<"HEPLike Warning, The <SubmissionYear> is missing in the yaml file"<<std::endl; 

  if(config["PublicationYear"]) HL_PublicationYear=config["PublicationYear"].as<std::string>();
  else if (HL_debug_yaml) std::cout<<"HEPLike Warning, The <PublicationYear> is missing in the yaml file"<<std::endl; 

  if(config["Collaborations"]) HL_Collaborations=config["Collaborations"].as<std::string>();
  else if (HL_debug_yaml) std::cout<<"HEPLike Warning, The <Collaborations> is missing in the yaml file"<<std::endl; 

  if(config["ExperimentalDataLumi"]) HL_ExperimentalDataLumi=config["ExperimentalDataLumi"].as<std::string>();
  else if (HL_debug_yaml) std::cout<<"HEPLike Warning, The <ExperimentalDataLumi> is missing in the yaml file"<<std::endl; 

  if(config["ExperimentalDataYears"]) HL_ExperimentalDataYears=config["ExperimentalDataYears"].as<std::string>();
  else if (HL_debug_yaml) std::cout<<"HEPLike Warning, The <ExperimentalDataYears> is missing in the yaml file"<<std::endl; 

  if(config["Decay"]) HL_Decay=config["Decay"].as<std::string>(); 
  else if (HL_debug_yaml) std::cout<<"HEPLike Warning, The <Decay> is missing in the yaml file"<<std::endl; 

  if(config["Kinematics"]) HL_Kinematics=config["Kinematics"].as<std::string>(); 
  else if (HL_debug_yaml) std::cout<<"HEPLike Warning, The <Kinematics> is missing in the yaml file"<<std::endl; 

  
  
}
void HL_Data::set_debug_yaml(bool b)
{
    HL_debug_yaml = b;

}
