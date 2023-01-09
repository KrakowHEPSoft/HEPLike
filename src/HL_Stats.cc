//   HEPLike: High Energy Physics Likelihoods
//
//   Module with statistic functions
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//   author: Tomas Gonzalo
//////////////////////////////////////////////////


#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "HL_Stats.h"
#include "HL_Constants.h"


//external:
#include "yaml-cpp/yaml.h"

#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_exp.h"

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace HL_Stats
{
  /// Minimum finite result returnable from log(double x);
  const double logmin = log(std::numeric_limits<double>::min());

  /// Use a detection to compute a simple chi-square-like log likelihood, for the case when obs is Gaussian distributed.
  double gauss(double obs, double theory, double theoryerr)
  {
    double norm=1./(sqrt(2.*HL_Const::pi)*theoryerr);
    double rest=pow(HL_Const::e,-((obs-theory)*(obs-theory))/(2.*theoryerr*theoryerr));
    return (norm*rest);
  }

  double gaussian_loglikelihood_theory_syst(double obs, double theory, double theoryerr)
  {
    double errsq = theoryerr*theoryerr;
    double chi2 = -0.5*pow(theory-obs,2)/errsq;
    double norm = 0.5*log(2.0*HL_Const::pi*errsq);
    return chi2 - norm;
  }


  double get_sigma_from_pval(double p)
  {
    double ip=0.;
    double nsigma=0.001;
    double dsigma=0.0001;

    while (ip< p) {
      ip= gsl_sf_erf(nsigma/M_SQRT2);
      nsigma+=dsigma;
    }
    return nsigma;

  }

  bool InvertMatrix (const ublas::matrix<double>& input, ublas::matrix<double>& inverse)
  {
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;

    // create a working copy of the input
    matrix<double> A(input);
    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A,pm);
    if ( res != 0 ) return false;

    // create identity matrix of "inverse"
    inverse.assign(identity_matrix<double>(A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
  }

}

