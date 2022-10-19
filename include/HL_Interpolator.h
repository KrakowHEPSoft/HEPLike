//   HEPLike: High Energy Physics Likelihoods
//
//   Module for interpolation
//
//   author: Tomas Gonzalo
//////////////////////////////////////////////////
#ifndef HL_INTERPOLATOR_H
#define HL_INTERPOLATOR_H

//C++ headers
#include <vector>
#include <string>

//external:
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>


class HL_Interpolator1D
{
  private:

    // Member variables
    gsl_spline* spline;
    gsl_interp_accel* x_accel;

    double x_min, x_max;

  public:

    // Constructor
    HL_Interpolator1D(int, double*, double*);

    // Destructor
    ~HL_Interpolator1D();

    // Evaluate a given interpolation
    double Eval(double) const;

    // Set limits
    void SetLimits(double, double);

    // Get limits
    double GetXmin() const;
    double GetXmax() const;

    #ifdef USE_ROOT
      TGraph TG;
      HL_Interpolator1D(TGraph tgraph);
    #endif
};

class HL_Interpolator2D
{
  private:

    // Member variables
    gsl_spline* spline2d;
    gsl_interp_accel* x_accel;
    gsl_interp_accel* y_accel;

    double x_min, xmax;
    double y_min, y_max;

  public:

    // Constructor
    HL_Interpolator2D(int, int, double*, double*, double*);

    // Destructor
    ~HL_Interpolator2D();

    // Evaluate a given interpolation
    double Eval(double, double) const;

    // Get limits
    double GetXmin() const;
    double GetXmax() const;
    double GetYmin() const;
    double GetYmax() const;

};

#endif
