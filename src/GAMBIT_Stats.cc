//   HEPLike: High Energy Physics Likelihoods
//
//   Module with statistic functions from GAMBIT
//
//////////////////////////////////////////////////

/* Copyright (c) 2017-2021, The GAMBIT Collaboration
   All rights reserved.

   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include <iostream>

#include "HL_Constants.h"
#include "GAMBIT_Stats.h"

namespace HL_Stats
{

  /// Minimum finite result returnable from log(double x);
  const double logmin = log(std::numeric_limits<double>::min());

  double gaussian_loglikelihood(double theory, double obs, double theoryerr, double obserr, bool profile_systematics)
  {
    if(theoryerr <0.) std::cerr<<"Negative theory error, in gaussian_loglikelihood"<<std::endl;
    if(obserr<0.) std::cerr<<"Negative observed error, in gaussian_loglikelihood"<<std::endl;

    if (theoryerr <= 0.) profile_systematics = false;

    double errsq = theoryerr*theoryerr + obserr*obserr;
    double chi2 = -0.5*pow(theory-obs,2)/errsq;
    double norm = profile_systematics ? log(2.0*HL_Const::pi*obserr*theoryerr) : 0.5*log(2.0*HL_Const::pi*errsq);

    return chi2 - norm;
  }


  /// Use a detection to compute a gaussian log-likelihood for an upper limit
  double gaussian_upper_limit(double theory, double obs, double theoryerr, double obserr, bool profile_systematics)
  {
    if (theoryerr < 0) std::cout << "Theory uncertainty cannot be negative." <<std::endl;
    if (obserr <= 0) std::cout << "Observational uncertainty must be non-zero and positive." <<std::endl;

    double loglike;
    double errsq = theoryerr*theoryerr + obserr*obserr;

    // Deal with the special case that no theory error has actually been given.  Revert to basic limit likelihood.
    if (theoryerr == 0)
    {
      double prefactor = -0.5*log(2.0*HL_Const::pi*errsq);
      return prefactor - (theory > obs ? 0.5*pow(theory-obs,2)/errsq : 0.0);
    }

    // Switch according to whether or not we are profiling.
    if (profile_systematics)
    {
      double prefactor = -log(2.0*HL_Const::pi*theoryerr*obserr);
      loglike = prefactor - (theory > obs? 0.5*pow(theory-obs,2)/errsq : 0.0);
    }
    else
    {
      double prefactor = -0.5*log(8.0*HL_Const::pi);
      double diff = obs - theory;
      double like = exp(-0.5*pow(diff,2)/errsq)/sqrt(errsq);
      like *= erfc(obserr * diff / (theoryerr * sqrt(2.0*errsq)));
      like += erfc(-diff/(sqrt(2.0)*theoryerr))/obserr;
      loglike = prefactor + (like == 0 ? logmin : log(like));
    }

    return loglike;

  }

}

