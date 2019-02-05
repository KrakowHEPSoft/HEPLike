//   HEPLike: High Energy Physics Likelihoods
//
//   Header with statistic functions
//
//   author: Marcin Chrzaszcz
//////////////////////////////////////////////////
#ifndef STATS_H
#define STATS_H

#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>

namespace HEPStats
{

  double gaussian_loglikelihood(double theory, double obs, double theoryerr, double obserr, bool profile_systematics);
  double lognormal_loglikelihood(double theory, double obs, double theoryerr, double obserr, bool profile_systematics);
  double lognormal_loglikelihood_relerror(double theory, double obs, double reltheoryerr, double relobserr, bool profile_systematics);
  double gaussian_upper_limit(double theory, double obs, double theoryerr, double obserr, bool profile_systematics);
  double gaussian_lower_limit(double theory, double obs, double theoryerr, double obserr, bool profile_systematics);
  double get_Likelihood_from_Cls(double cls);

}

#endif
