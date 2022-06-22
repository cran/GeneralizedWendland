#include <Rcpp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

// Build gsl_function from lambda
template <typename F>
class gsl_function_pp : public gsl_function{
  const F func;

  static double invoke(double x, void *params){
    return static_cast<gsl_function_pp*>(params)->func(x);
  }

public:
  gsl_function_pp(const F &f) : func(f){
    function = &gsl_function_pp::invoke;
    params   = this;
  }
  operator gsl_function*(){return this;}
};

template <typename F>
gsl_function_pp<F> make_gsl_function(const F &func){
  return gsl_function_pp<F>(func);
};

class Integrator
{
private:
  gsl_integration_workspace* workspace = nullptr;

protected:
  template <typename Func>
  double integrateQNG(double lower, double upper, const Func &func){
    auto integrand = make_gsl_function(func);
    gsl_set_error_handler_off();
    status.error = gsl_integration_qng(integrand, lower, upper, abstol, reltol,
      &status.result, &status.abserr, &status.neval);

    if (status.error) Rcpp::stop("Error during QNG integration");
    return status.result;
  }

  template <typename Func>
  double integrateQAG(double lower, double upper, const Func &func){
    status.neval = intervals;
    auto integrand = make_gsl_function(func);
    gsl_set_error_handler_off();
    status.error = gsl_integration_qag(integrand, lower, upper, abstol, reltol,
      intervals, qag_key, workspace, &status.result, &status.abserr);

    if (status.error) Rcpp::stop("Error during QAG integration");
    return status.result;
  }

  template <typename Func>
  double integrateQAGS(double lower, double upper, const Func &func){
    status.neval = intervals;
    auto integrand = make_gsl_function(func);
    gsl_set_error_handler_off();
    status.error = gsl_integration_qags(integrand, lower, upper, abstol, reltol,
      intervals, workspace, &status.result, &status.abserr);

    if (status.error) Rcpp::stop("Error during QAGS integration");
    return status.result;
  }

public:

  struct{
    double result;
    double abserr;
    int error;
    size_t neval;
  } status;

  double abstol;
  double reltol;
  int intervals;
  int qag_key;
  int method;


  Integrator(double _atol, double _rtol, int _int, int _key) : abstol(_atol),
    reltol(_rtol), intervals(_int), qag_key(_key){
    if (intervals) workspace = gsl_integration_workspace_alloc(intervals);
  };

  ~Integrator(){
    if (workspace) {
      gsl_integration_workspace_free(workspace);
      gsl_integration_workspace* workspace = nullptr;
    }
  }

  template <typename Func>
  double integrate(double lower, double upper, const Func &func){
    if (!workspace){
      return integrateQNG(lower, upper, func);
    } else if (qag_key) {
      return integrateQAG(lower, upper, func);
    } else {
      return integrateQAGS(lower, upper, func);
    }
  }
};

