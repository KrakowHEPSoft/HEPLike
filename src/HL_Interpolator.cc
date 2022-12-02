//   HEPLike: High Energy Physics Likelihoods
//
//   Module for interpolation
//
//   author: Tomas Gonzalo
//////////////////////////////////////////////////

#include <stdexcept>
#include <iostream>
#include <algorithm>

#include "HL_Interpolator.h"


HL_Interpolator1D::HL_Interpolator1D(size_t npoints, double x[], double y[])
{
  nx = npoints;
  x_data = new double[nx]; std::copy(x,x+nx,x_data);
  y_data = new double[nx]; std::copy(y,y+nx,y_data);
  x_accel = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_linear, npoints);
  gsl_spline_init(spline, x, y, npoints);
}

HL_Interpolator1D::HL_Interpolator1D(const HL_Interpolator1D &interp)
 : nx(interp.nX())
 , x_min(interp.GetXmin())
 , x_max(interp.GetXmax())
{
  x_data = new double[nx]; std::copy(interp.x_data,interp.x_data+nx,x_data);
  y_data = new double[nx]; std::copy(interp.y_data,interp.y_data+nx,y_data);

  x_accel = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_linear, nx);
  gsl_spline_init(spline, x_data, y_data, nx);
}

HL_Interpolator1D &HL_Interpolator1D::operator=(const HL_Interpolator1D &interp)
{
  nx = interp.nX();
  x_min = interp.GetXmin();
  x_max = interp.GetXmax();

  x_data = new double[nx]; std::copy(interp.x_data,interp.x_data+nx,x_data);
  y_data = new double[nx]; std::copy(interp.y_data,interp.y_data+nx,y_data);

  x_accel = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_linear, nx);
  gsl_spline_init(spline, x_data, y_data, nx);

  return *this;
}


#ifdef USE_ROOT
 HL_Interpolator1D::HL_Interpolator1D(TGraph *tgraph)
 {
   TG = tgraph;
 }
#endif

HL_Interpolator1D::~HL_Interpolator1D()
{
  if(spline) gsl_spline_free(spline);
  if(x_accel) gsl_interp_accel_free(x_accel);

  delete[] x_data;
  delete[] y_data;

  #ifdef USE_ROOT
    delete TG;
  #endif
}

double HL_Interpolator1D::Eval(double x) const
{
  if(x < x_min or x > x_max)
    throw std::runtime_error("Variable outside of interpolation range.");

  #ifdef USE_ROOT
    return TG->Eval(x, 0);
  #else
    return gsl_spline_eval(spline, x, x_accel);
  #endif
}

void HL_Interpolator1D::SetLimits(double xmin, double xmax)
{
  x_min = xmin;
  x_max = xmax;

  #ifdef USE_ROOT
    TG->GetXaxis()->SetLimits(xmin,xmax);
  #endif
}

size_t HL_Interpolator1D::nX() const
{
  return nx;
}

double HL_Interpolator1D::GetXmin() const
{
  return x_min;
}

double HL_Interpolator1D::GetXmax() const
{
  return x_max;
}

HL_Interpolator2D::HL_Interpolator2D(size_t npointsx, size_t npointsy, double x[], double y[], double z[])
: nx(npointsx)
, ny(npointsy)
{
  x_data = new double[nx]; std::copy(x,x+nx,x_data);
  y_data = new double[ny]; std::copy(y,y+ny,y_data);
  z_data = new double[nx*ny]; std::copy(z,z+nx*ny,z_data);

  x_accel = gsl_interp_accel_alloc();
  y_accel = gsl_interp_accel_alloc();
  spline2d = gsl_spline2d_alloc(gsl_interp2d_bilinear, npointsx, npointsy);
  gsl_spline2d_init(spline2d, x, y, z, npointsx, npointsy);

}

HL_Interpolator2D::HL_Interpolator2D(const HL_Interpolator2D &interp)
 : nx(interp.nX())
 , ny(interp.nY())
 , x_min(interp.GetXmin())
 , x_max(interp.GetXmax())
 , y_min(interp.GetYmin())
 , y_max(interp.GetYmax())
{
  x_data = new double[nx]; std::copy(interp.x_data,interp.x_data+nx,x_data);
  y_data = new double[ny]; std::copy(interp.y_data,interp.y_data+ny,y_data);
  z_data = new double[nx*ny]; std::copy(interp.z_data,interp.z_data+nx*ny,z_data);

  x_accel = gsl_interp_accel_alloc();
  y_accel = gsl_interp_accel_alloc();
  spline2d = gsl_spline2d_alloc(gsl_interp2d_bilinear, nx, ny);
  gsl_spline2d_init(spline2d, x_data, y_data, z_data, nx, ny);
}

HL_Interpolator2D &HL_Interpolator2D::operator=(const HL_Interpolator2D &interp)
{
  nx = interp.nX();
  ny = interp.nY();
  x_min = interp.GetXmin();
  x_max = interp.GetXmax();
  y_min = interp.GetYmin();
  y_max = interp.GetYmax();

  x_data = new double[nx]; std::copy(interp.x_data,interp.x_data+nx,x_data);
  y_data = new double[ny]; std::copy(interp.y_data,interp.y_data+ny,y_data);
  z_data = new double[nx*ny]; std::copy(interp.z_data,interp.z_data+nx*ny,z_data);

  x_accel = gsl_interp_accel_alloc();
  y_accel = gsl_interp_accel_alloc();
  spline2d = gsl_spline2d_alloc(gsl_interp2d_bilinear, nx, ny);
  gsl_spline2d_init(spline2d, x_data, y_data, z_data, nx, ny);

  return *this;
}

#ifdef USE_ROOT
 HL_Interpolator2D::HL_Interpolator2D(TH2D *th2d)
 {
   TH = th2d;
 }
#endif

HL_Interpolator2D::~HL_Interpolator2D()
{
  if(spline2d) gsl_spline2d_free(spline2d);
  if(x_accel) gsl_interp_accel_free(x_accel);
  if(y_accel) gsl_interp_accel_free(y_accel);

  delete[] x_data;
  delete[] y_data;
  delete[] z_data;

  #ifdef USE_ROOT
    delete TH;
  #endif
}

double HL_Interpolator2D::Eval(double x, double y) const
{
  if(x < x_min or x > x_max or y < y_min or y > y_max)
    throw std::runtime_error("Variable outside of interpolation range.");

  #ifdef USE_ROOT
    return TH->Interpolate(x, y);
  #else
    return gsl_spline2d_eval(spline2d, x, y, x_accel, y_accel);
  #endif
}

void HL_Interpolator2D::SetLimits(double xmin, double xmax, double ymin, double ymax)
{
  x_min = xmin;
  x_max = xmax;
  y_min = ymin;
  y_max = ymax;
}

size_t HL_Interpolator2D::nX() const
{
  return nx;
}

size_t HL_Interpolator2D::nY() const
{
  return ny;
}

double HL_Interpolator2D::GetXmin() const
{
  return x_min;
}

double HL_Interpolator2D::GetXmax() const
{
  return x_max;
}

double HL_Interpolator2D::GetYmin() const
{
  return y_min;
}

double HL_Interpolator2D::GetYmax() const
{
  return y_max;
}

