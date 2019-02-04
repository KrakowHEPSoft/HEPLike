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

#include "stats.h"
#include "constants.h"
#include "HepData.h"


void HEPData:: HEPData()
{
  initialized=false;
  HFile="";

}
void HEPData:: HEPData(std::string s)
{

  initialized=true;
  HFile=s;
}



void HEPData::read()
{
  if(! initialized)
    {
      std::std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }
  YAML::Node config = YAML::LoadFile(HFile);
  HEPDOI=config["DOI"];
  HEPBibCite=config["BibCite"];
  HEPBibEntry=config["BibEntry"];
  HEPFileName=config["FileName"];
  HEPHFAG=config["HFAG"];
  HEPSource=config["HFAG"];
  HEPYear=config["Year"];
  HEPName=config["Name"];
  
}
    

  
