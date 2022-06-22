// [[Rcpp::depends(RcppEigen, BH)]]
#pragma once

#include <limits>
#include <R.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/math/special_functions/beta.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "GSL_Integration.h"
#include "GSL_Interpolation.h"

class WendlandParam
{

protected:
  void checkParameters();

public:
  WendlandParam();
  WendlandParam(double _range, double _sill, double _kappa, double _mu, double _nugget);

  double range = 1.0;
  double sill = 1.0;
  double kappa = 0.0;
  double mu = 2.5;
  double nugget = 0.0;
};

class Wendland
{

protected:
  WendlandParam param = WendlandParam();
  Integrator* integrator = nullptr;
  Interpolator<Wendland>* interpolator = nullptr;
  double epstol = std::numeric_limits<double>::epsilon();
  double beta_constant;
  void computeBetaConstant();
  double computeIntegral(double distance);

public:
  Wendland();
  ~Wendland();

  void setParameters(double range, double sill, double kappa, double mu, double nugget);
  void setEpsTol(double _epstol);
  void setIntegrator(double abstol, double reltol, int intervals, int qag_key);
  void setInterpolator(int num_points, int interp_type);
  void deleteIntegrator();
  void deleteInterpolator();
  double compute(const double& distance);

  Rcpp::NumericVector computeVector(const Rcpp::NumericVector& dvec);
  Eigen::MatrixXd computeMatrix(const Eigen::MatrixXd& dmat);
  Eigen::SparseMatrix<double> computeMSparse(const Eigen::SparseMatrix<double>& dmat);
  Rcpp::List computeSpam(const Eigen::MatrixXi& index, const Eigen::VectorXd& values);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
