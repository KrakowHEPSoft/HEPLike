//   HEPLike: High Energy Physics Likelihoods
//
//   Module with statistic functions
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
namespace HEPStats
{
    
  /// Use a detection to compute a simple chi-square-like log likelihood, for the case when obs is Gaussian distributed. 
  double gaussian_loglikelihood(double theory, double obs, double theoryerr, double obserr)
  {
    if(theoryerr <0.) std::cerr<<"Negative theory error, in gaussian_loglikelihood"<<std::endl;
    if(obserr<0.) std::cerr<<"Negative error error, in gaussian_loglikelihood"<<std::endl; 

    double errsq = theoryerr*theoryerr + obserr*obserr;
    double chi2 = pow(theory-obs,2)/errsq;     
    double loglikelihood=-0.5*chi2;

    return loglikelihood;
  }
  /// Use a detection to compute a simple chi-square-like log likelihood, for the case when obs is log-normal distributed.
  /// Version that takes absolute errors.
  double lognormal_loglikelihood(double theory, double obs, double theoryerr, double obserr, bool profile_systematics)
  {
    return lognormal_loglikelihood_relerror(theory, obs, theoryerr/theory, obserr/obs, profile_systematics);
  }
  /*
  /// Use a detection to compute a simple chi-square-like log likelihood, for the case when obs is log-normal distributed.
  /// Version that takes fractional (relative) errors.
  double lognormal_loglikelihood_relerror(double theory, double obs, double reltheoryerr, double relobserr, bool profile_systematics)
  {
    if (reltheoryerr < 0) utils_error().raise(LOCAL_INFO, "Theory uncertainty cannot be negative.");
    if (relobserr <= 0) utils_error().raise(LOCAL_INFO, "Observational uncertainty must be non-zero and positive.");

    // Don't allow observed or predicted values of less than or equal to 0
    if (obs <= 1e2*std::numeric_limits<double>::min()) return std::numeric_limits<double>::lowest();
    if (theory <= 1e2*std::numeric_limits<double>::min()) return std::numeric_limits<double>::lowest();

    // If no theory error is given, go back to the standard log-normal likelihood without systematics.
    if (reltheoryerr == 0) profile_systematics = false;

    double log_obs_on_theory = log(obs / theory);
    double obserr_prime = log1p(relobserr);
    double theoryerr_prime = log1p(reltheoryerr);
    double errsq = obserr_prime*obserr_prime + theoryerr_prime*theoryerr_prime;
    double chi2 = -log(obs) - 0.5*pow(log_obs_on_theory,2)/errsq;
    double norm = profile_systematics ? log(2.0*HEPConst::pi*obserr_prime*theoryerr_prime) : 0.5*log(2.0*HEPConst::pi*errsq);
    double extra = profile_systematics ? theoryerr_prime*theoryerr_prime*(0.5*obserr_prime*obserr_prime - log_obs_on_theory)/errsq : 0.0;
    return chi2 + extra - norm;
  }
  /// Use a detection to compute a gaussian log-likelihood for an upper limit
  double gaussian_upper_limit(double theory, double obs, double theoryerr, double obserr, bool profile_systematics)
  {
    if (theoryerr < 0) utils_error().raise(LOCAL_INFO, "Theory uncertainty cannot be negative.");
    if (obserr <= 0) utils_error().raise(LOCAL_INFO, "Observational uncertainty must be non-zero and positive.");

    double loglike;
    double errsq = theoryerr*theoryerr + obserr*obserr;

    // Deal with the special case that no theory error has actually been given.  Revert to basic limit likelihood.
    if (theoryerr == 0)
      {
        double prefactor = -0.5*log(2.0*HEPConst::pi*errsq);
        return prefactor - (theory > obs ? 0.5*pow(theory-obs,2)/errsq : 0.0);
      }

    // Switch according to whether or not we are profiling.
    if (profile_systematics)
      {
        double prefactor = -log(2.0*HEPConst::pi*theoryerr*obserr);
        loglike = prefactor - (theory > obs? 0.5*pow(theory-obs,2)/errsq : 0.0);
      }
    else
      {
        double prefactor = -0.5*log(8.0*HEPConst::pi);
        double diff = obs - theory;
        double like = exp(-0.5*pow(diff,2)/errsq)/sqrt(errsq);
        like *= erfc(obserr * diff / (theoryerr * sqrt(2.0*errsq)));
        like += erfc(-diff/(sqrt(2.0)*theoryerr))/obserr;
        if (like < 0) utils_error().raise(LOCAL_INFO, "Marginalised Gaussian limit likelihood went negative!");
        loglike = prefactor + (like == 0 ? logmin : log(like));
      }

    return loglike;

  }

  /// Use a detection to compute a gaussian log-likelihood for a lower limit
  double gaussian_lower_limit(double theory, double obs, double theoryerr, double obserr, bool profile_systematics)
  {
    return gaussian_upper_limit(-theory, -obs, theoryerr, obserr, profile_systematics);
  }
  */
}
  
