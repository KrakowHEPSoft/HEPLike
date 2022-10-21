//   HEPLike: High Energy Physics Likelihoods
//
//   Module for function wrappers
//
//   author: Tomas Gonzalo
//////////////////////////////////////////////////
#ifndef HL_FUNCTION_H
#define HL_FUNCTION_H

#include "HL_Stats.h"
#include "HL_Interpolator.h"

#include "gsl/gsl_vector.h"

class HL_Function
{
  public:
    virtual ~HL_Function() {}
    virtual double operator()(const gsl_vector *, void *) { throw std::runtime_error("Not implemented"); }
};

class HL_Function1D : public HL_Function
{

  private:

    double theory_mean;
    double theory_err;
    HL_Interpolator1D *likelihood;


  public:

    HL_Function1D(HL_Interpolator1D *l)
    {
      likelihood = l;
    };

    double operator()(const gsl_vector *x, void *)
    {
      double theory_nuisance = gsl_vector_get(x, 0);

      double loglike = likelihood->Eval(theory_nuisance);
      //double like = exp(loglike);


      double gauss_systematic=HL_Stats::gauss(theory_nuisance, theory_mean, theory_err);

      return loglike-log(gauss_systematic);// here the loglike is -\Delta LL so no minus before
    }


    void SetTheory(double mean, double err)
    {
      theory_mean=mean;
      theory_err=err;
    };

};

class HL_Function2D : public HL_Function
{

  private:


    vector <double> theory_mean;
    boost::numeric::ublas::matrix<double> theory_cov;
    HL_Interpolator2D *likelihood;

  public:

    HL_Function2D(HL_Interpolator2D *l)
    {
      likelihood = l;
    };

    double operator()(const gsl_vector *x, void *)
    {
      std::vector<double> theory_nuisance = {gsl_vector_get(x, 0), gsl_vector_get(x,1)};

      double loglikelihood_penalty=1e10;

      double xmin=likelihood->GetXmin();
      double xmax=likelihood->GetXmax();
      double ymin=likelihood->GetYmin();
      double ymax=likelihood->GetYmax();
      if(theory_nuisance[0]>xmax) return loglikelihood_penalty;
      if(theory_nuisance[0]<xmin) return loglikelihood_penalty;
      if(theory_nuisance[1]>ymax) return loglikelihood_penalty;
      if(theory_nuisance[1]<ymin) return loglikelihood_penalty;

      double negative_LL = likelihood->Eval(theory_nuisance[0], theory_nuisance[1]);

      double chi2=0.;
      double diff[2];
      for(unsigned i=0; i < 2; i++)
      {
        diff[i]=theory_nuisance[i]-theory_mean[i];
      }
      boost::numeric::ublas::matrix<double>theory_cov_inv(2,2);

      HL_Stats::InvertMatrix(theory_cov, theory_cov_inv);


      for(unsigned i=0; i <2 ; i++)
      {
        for(unsigned j=0; j< 2 ;  j++)
        {
          chi2+=diff[i] * theory_cov_inv(i,j)*diff[j] ;
        }
      }
      double deltaLL= 0.5*chi2;

      return (deltaLL+negative_LL); // here both are negative likelihoods that needs to be mimalized

    }

    void SetTheory(vector<double> mean, boost::numeric::ublas::matrix<double> cov)
    {
      theory_mean=mean;
      theory_cov=cov;
    }
};

#endif
