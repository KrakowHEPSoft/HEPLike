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
#include "HEPBRLimit.h"



void HEPBRLimit::read()
{
  if(! initialized)
    {
      std::cout << "TRYING TO READ WITHOUT GIVING ANY FILE!" << std::endl;
      return;
    }
  YAML::Node config = YAML::LoadFile(HFile);
  HEPDOI=config["DOI"].as<std::string>();
  HEPBibCite=config["BibCite"].as<std::string>();
  HEPBibEntry=config["BibEntry"].as<std::string>();
  HEPFileName=config["FileName"].as<std::string>();
  HEPHFAG=config["HFAG"].as<std::string>();
  HEPSource=config["HFAG"].as<std::string>();
  HEPYear=config["Year"].as<std::string>();
  HEPName=config["Name"].as<std::string>();
  YAML::Node node  = config["Cls"];
  for(YAML::const_iterator it = node.begin(); it != node.end();  ++it )
    {
      std::cout << *it << std::endl;
      std::cout << (*it)[0] <<  std::endl;   
      CLs.push_back( ((*it)[0]).as<double>()  );
      BR.push_back( ((*it)[1]).as<double>()  ); 
    }

  
}
    

  
