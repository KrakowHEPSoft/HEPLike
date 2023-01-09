//   HEPLike: High Energy Physics Likelihoods
//
//   Module for minimizing
//
//   author: Tomas Gonzalo
//////////////////////////////////////////////////
#ifndef HL_MINIMIZER_H
#define HL_MINIMIZER_H

#include <string>

#include "HL_Function.h"

#include "gsl/gsl_vector.h"
#include "gsl/gsl_multimin.h"

class HL_Minimizer
{
  private:

    size_t ndim;
    size_t maxiters;
    double step_size;
    double tolerance;

  public:

    gsl_multimin_fdfminimizer *s;
    gsl_multimin_function_fdf my_func;
    gsl_vector *x;

    HL_Minimizer(std::string, size_t);
    HL_Minimizer(const HL_Minimizer&);
    HL_Minimizer &operator=(const HL_Minimizer&);

    ~HL_Minimizer();

    void SetMaxIterations(const size_t);
    void SetTolerance(const double);

    void SetFunction(HL_Function&);

    void SetVariable(int, double, double);

    size_t GetNDim() const;
    size_t GetMaxIters() const;
    double GetStepSize() const;
    double GetTolerance() const;

    void Minimize();

    std::vector<double> X();
};

#endif
