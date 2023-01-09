//   HEPLike: High Energy Physics Likelihoods
//
//   Header for ndim Likelihood
//
//   author: Jihyun Bhom, Marcin Chrzaszcz
//   author: Tomas Gonzalo
//////////////////////////////////////////////////


#ifndef HL_NDIMLIKELIHOOOD_H
#define HL_NDIMLIKELIHOOOD_H

//C++ headers
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>


//HL_Like headers
#include "HL_Stats.h"
#include "HL_Data.h"
#include "HL_Interpolator.h"
#include "HL_Function.h"
#include "HL_Minimizer.h"

//external:
#include "yaml-cpp/yaml.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_exp.h"
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

class HL_nDimLikelihood: public HL_Data
{

  public:

    explicit HL_nDimLikelihood() : HL_Data(), hist2D(nullptr), hist_profileX(nullptr), hist_profileY(nullptr), gmin(nullptr), fun(nullptr) {};
    explicit HL_nDimLikelihood(std::string s) :  HL_Data(s), hist2D(nullptr), hist_profileX(nullptr), hist_profileY(nullptr), gmin(nullptr), fun(nullptr) {};
    HL_nDimLikelihood(const HL_nDimLikelihood &);
    HL_nDimLikelihood &operator=(const HL_nDimLikelihood &);


    ~HL_nDimLikelihood();


    void Read();
    double GetChi2( std::vector<double> theory) ;
    double GetChi2( std::vector<double> theory,  boost::numeric::ublas::matrix<double> theory_cov);
    double GetLikelihood( std::vector<double> theory) ;
    double GetLikelihood(std::vector<double> theory, boost::numeric::ublas::matrix<double> theory_cov);
    double GetLogLikelihood(  std::vector<double> theory) ;
    double GetLogLikelihood(std::vector<double> theory, boost::numeric::ublas::matrix<double> theory_cov);


    void Profile(std::string="max");
    double GetChi2_profile( double theory, std::string);
    double GetLikelihood_profile( double theory, std::string axis) ;
    double GetLogLikelihood_profile(  double theory, std::string X);

    std::vector<std::string> GetObservables(){ return Observables;};

    double loglikelihood_penalty;

    std::string HL_RootFile;
    std::string HL_PATH;
    std::vector<std::string> Observables;

    int NoOfObservables;
    int size_restricted;

    HL_Interpolator2D *hist2D;
    //HL_Interpolator3D *hist3D;

    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;

    std::vector<double> central_mes_val;
    int dim;

    HL_Interpolator1D *hist_profileX;
    HL_Interpolator1D *hist_profileY;
    //HL_Interpolator1D *hist_profileZ;

    int n_binsX;
    int n_binsY;
    int n_binsZ;


    bool profiled;

    HL_Minimizer *gmin;
    HL_Function2D *fun;

};


#endif
