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
  step_size = 0.;

  if(type == "ConjugateFR")
    s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, n);
  else if (type == "ConjugatePR")
    s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_pr, n);
  else
    throw std::runtime_error("Unkown minimizer algorithm");

  x = gsl_vector_alloc(ndim);
}

HL_Minimizer::HL_Minimizer(const HL_Minimizer &min)
 : ndim(min.GetNDim())
 , maxiters(min.GetMaxIters())
 , step_size(min.GetStepSize())
 , tolerance(min.GetTolerance())
 , my_func(min.my_func)
{
  s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, ndim);
  x = gsl_vector_alloc(ndim);
}

HL_Minimizer &HL_Minimizer::operator=(const HL_Minimizer &min)
{
  ndim = min.GetNDim();
  maxiters = min.GetMaxIters();
  step_size = min.GetStepSize();
  tolerance = min.GetTolerance();
  my_func = min.my_func;

  s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, ndim);
  x = gsl_vector_alloc(ndim);

  return *this;
}

HL_Minimizer::~HL_Minimizer()
{
  if(s) gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);
}

void HL_Minimizer::SetMaxIterations(const size_t iters)
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

  const void * p = &f;
  assert (p != 0);
  my_func.f  = &HL_Function::F;
  my_func.df = &HL_Function::DF;
  my_func.fdf = &HL_Function::FDF;
  my_func.params = const_cast<void *>(p);

}

void HL_Minimizer::SetVariable(int i, double value, double step)
{
  gsl_vector_set (x, i, value);
  if(!step_size or step_size > step) step_size = step;
}

size_t HL_Minimizer::GetNDim() const
{
  return ndim;
}

size_t HL_Minimizer::GetMaxIters() const
{
  return maxiters;
}

double HL_Minimizer::GetStepSize() const
{
  return step_size;
}

double HL_Minimizer::GetTolerance() const
{
  return tolerance;
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
