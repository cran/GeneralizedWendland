#include <Rcpp.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>

template<typename CovarianceFunction>
class Interpolator
{

private:
  const gsl_interp_type* type = nullptr;

protected:
  CovarianceFunction* covfun = nullptr;
  gsl_interp *workspace_pointer = nullptr;
  gsl_interp_accel *accel_pointer = nullptr;

  double *evaluation_points;
  double *evaluation_results;
  int method;
  int num_points;
  double upper_limit;
  bool initialized = false;

  double interpolate(const double& distance){
    return gsl_interp_eval(workspace_pointer, evaluation_points, evaluation_results, distance, accel_pointer);
  }

public:
  Interpolator(CovarianceFunction* _covfun, int _method, int _num_points, double _upper_limit)
  : covfun(_covfun), method(_method), num_points(_num_points), upper_limit(_upper_limit){
    if (method == 1){
      type = gsl_interp_linear;
    } else if (method == 2) {
      type = gsl_interp_polynomial;
    } else if (method == 3) {
      type = gsl_interp_cspline;
    } else {
      Rcpp::stop("Undefined method.");
    }

    initialize();
  };

  ~Interpolator(){
    gsl_interp_free(workspace_pointer);
    gsl_interp_accel_free(accel_pointer);
    delete[] evaluation_points;
    delete[] evaluation_results;
  }

  void initialize(){
    initialized = false;
    double stepsize = upper_limit/(num_points - 1.0);
    evaluation_points = new double [num_points];
    evaluation_results = new double [num_points];
    workspace_pointer = gsl_interp_alloc(type, num_points);
    accel_pointer = gsl_interp_accel_alloc();

    for (int i = 0; i < num_points; i++){
      evaluation_points[i] = stepsize * i;
      evaluation_results[i] = covfun->compute(stepsize * i);
    }

    gsl_interp_init(workspace_pointer, evaluation_points, evaluation_results, num_points);
    gsl_set_error_handler_off();
    initialized = true;
  }

  double compute(const double& distance){
    return (distance < upper_limit) ? interpolate(distance) : 0.0;
  }

  bool isInitialized(){
    return initialized;
  }
};
