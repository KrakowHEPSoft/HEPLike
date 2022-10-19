//   HEPLike: High Energy Physics Likelihoods
//
//   Module for minimizing
//
//   author: Tomas Gonzalo
//////////////////////////////////////////////////

#include "HL_Minimizer.h"

HL_Minimizer::HL_Minimizer(std::string type, size_t n)
{
  ndim = n;

  if(type == "ConjugateFR")
    s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, n);

  x = gsl_vector_alloc(ndim);
}

HL_Minimizer::~HL_Minimizer()
{
  gsl_min_fdfminimizer_free(s);
  gsl_vector_free (x);
}

void HL_Minimizer::SetMaxIterations(const int iters)
{
  maxiters = iters;
}

void HL_Minimizer::SetTolerance(const double tol)
{
  tolerance = tol;
}

void HL_Minimizer::SetFunction(HL_Function &f)
{
  my_func.n = ndim;
  my_func.f = &(f->DoEval);
}

void HL_Minimizer::SetVariable(int i, double value, double step)
{
  gsl_vector_set (x, i, value);
  if(step_size > step) step_size = step;
}

void HL_Minimizer::Minimize()
{
  size_t iter = 0;
  int status;

  gsl_multimin_fdfminimizer_set (s, &my_func, x, step_size, tolerance);

  do
  {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate (s);

    if (status)
      break;

    status = gsl_multimin_test_gradient (s->gradient, tolerance);
  }
  while (status == GSL_CONTINUE && iter < maxiters);
}

std::vector<double> HL_Minimizer::X()
{
  std::vector<double> X;
  for(size_t i = 0; i < ndim; ++i)
    X.push_back(gsl_vector_get(s->x, i));
  return X;
}
