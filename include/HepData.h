//   HEPLike: High Energy Physics Likelihoods
//
//   Header for limit class
//
//   author: Marcin Chrzaszcz
//////////////////////////////////////////////////


#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>



#include "stats.h"



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
  std::string HEPname;

      
  bool initialized;
  
}
