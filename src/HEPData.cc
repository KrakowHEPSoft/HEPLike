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
  HEPDOI=config["DOI"].as<std::string>();
  HEPBibCite=config["BibCite"].as<std::string>();
  HEPBibEntry=config["BibEntry"].as<std::string>();
  HEPFileName=config["FileName"].as<std::string>();
  HEPHFLAV=config["HFLAV"].as<std::string>();
  HEPSource=config["Source"].as<std::string>();
  HEPYear=config["Year"].as<std::string>();
  HEPName=config["Name"].as<std::string>();
  HEPDecay=config["Decay"].as<std::string>();
}
