//   HEPLike: High Energy Physics Likelihoods
//
//   Header with statistic functions
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//   author: Tomas Gonzalo
//////////////////////////////////////////////////
#ifndef HL_STATS_H
#define HL_STATS_H

#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "GAMBIT_Stats.h"

using namespace std;

namespace HL_Stats
{
  namespace ublas = boost::numeric::ublas;

  double gauss(double obs, double theory, double theoryerr);
  double gaussian_loglikelihood_theory_syst(double obs, double theory, double theoryerr);
  double gauss(double x, double mean, double sigma);

  double get_sigma_from_pval(double p);
  template<class T>
  bool InvertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse)
  {
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;

    // create a working copy of the input
    matrix<T> A(input);
    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A,pm);
    if ( res != 0 )
    {
      throw std::runtime_error("The matrix inversion failed in HL_Stats::InvertMatrix");
      return false;
    }
    // create identity matrix of "inverse"
    inverse.assign(identity_matrix<T>(A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
  }

}

#endif
