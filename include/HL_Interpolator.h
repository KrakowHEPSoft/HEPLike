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
#include <gsl/gsl_spline2d.h>


class HL_Interpolator1D
{
  private:

    // Member variables
    size_t nx;
    double x_min, x_max;

  public:

    double *x_data, *y_data;

    gsl_spline* spline;
    gsl_interp_accel* x_accel;

    // Constructors
    HL_Interpolator1D(size_t, double[], double[]);
    HL_Interpolator1D(const HL_Interpolator1D&);
    HL_Interpolator1D &operator=(const HL_Interpolator1D&);

    // Destructor
    ~HL_Interpolator1D();

    // Evaluate a given interpolation
    double Eval(double) const;

    // Set limits
    void SetLimits(double, double);

    // Get number of points
    size_t nX() const;

    // Get limits
    double GetXmin() const;
    double GetXmax() const;

    #ifdef USE_ROOT
      TGraph *TG;
      HL_Interpolator1D(TGraph*);
    #endif
};

class HL_Interpolator2D
{
  private:

    // Member variables
    size_t nx, ny;
    double x_min, x_max;
    double y_min, y_max;

  public:

    double *x_data, *y_data, *z_data;

    gsl_spline2d* spline2d;
    gsl_interp_accel* x_accel;
    gsl_interp_accel* y_accel;

    // Constructors
    HL_Interpolator2D(size_t, size_t, double[], double[], double[]);
    HL_Interpolator2D(const HL_Interpolator2D&);
    HL_Interpolator2D &operator=(const HL_Interpolator2D&);

    // Destructor
    ~HL_Interpolator2D();

    // Evaluate a given interpolation
    double Eval(double, double) const;

    // Set limits
    void SetLimits(double, double, double, double);

    // Get number of points
    size_t nX() const;
    size_t nY() const;

    // Get limits
    double GetXmin() const;
    double GetXmax() const;
    double GetYmin() const;
    double GetYmax() const;

    #ifdef USE_ROOT
      TH2D *TH;
      HL_Interpolator2D(TH2D*);
    #endif

};

#endif
