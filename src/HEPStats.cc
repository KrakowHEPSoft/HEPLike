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

#include "HEPStats.h"
#include "HEPConstants.h"


//external:
#include "yaml-cpp/yaml.h"

#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_exp.h"

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>


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
  
