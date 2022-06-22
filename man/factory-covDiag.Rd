\name{factory-covDiag}

\alias{covDiagFactory}
\alias{accumulated_error}
\alias{point_diagnostics}

\title{Diagnostics For Arbitrary Covariance Functions}

\description{
A suite of diagnostic tools. The functions described here provide the user with quick access to diagnostics for arbitrary target covariance functions and arbitrary reference covariance functions.
}

\usage{
covDiagFactory(target_covariance,
  diagnostic_funs = c("accumulated_error", "point_diagnostics"),
  reference_covariance = cov.askey, reference_cov.args = list())

accumulated_error(target_covFun, target_cov.theta, reference_covFun,
  reference_cov.theta, ..., absolute = TRUE, lower = 0, upper = 1,
  subdivisions = 500L, abs.tol = .Machine$double.eps^0.5,
  rel.tol = .Machine$double.eps^0.25)

point_diagnostics(target_covFun, target_cov.theta, reference_covFun,
  reference_cov.theta, ..., grid_resolution = 100)
}

\value{
\code{covDiagFactory()} is a function factory which generates a function for computing diagnostics across a large variety of parameter values and additional arguments. The manufactured function has the form \code{function(target_theta_list, target_args_list = list(), reference_cov.theta = NULL, ...} and returns a list of the same length as \code{diagnostic_funs}, each entry storing a \code{data.frame} with the results returned by the associated diagnostic function.
}

\arguments{
\item{diagnostic_funs}{A character vector for specifying which diagnostics to compute and return. Currently, \code{accumulated_error} and \code{point_diagnostics} are implemented. Note that the functions themselves are not intended to be called directly by the user.}

\item{reference_covariance}{Covariance function with formals \code{function(h, theta, ..., cov.args)}}

\item{reference_covFun}{Same as reference_covariance, but implicitly assumes that the input was generated using \code{covarianceFactory()}.}

\item{reference_cov.args}{List with additional arguments to be passed to reference_covariance.}

\item{reference_cov.theta}{Numeric vector containing parameters for reference_covariance.}

\item{target_covariance}{Covariance function with formals \code{function(h, theta, ..., cov.args)}}

\item{target_cov.theta}{Numeric vector containing parameters for target_covariance.}

\item{target_covFun}{Same as target_covariance, but implicitly assumes that the input was generated using \code{covarianceFactory()}.}

\item{...}{Other arguments}

\item{grid_resolution}{Number of points \eqn{n} to evaluate the covariance function at.}

\item{absolute}{Logical value. Whether to return absolute value.}

\item{lower}{Lower boundary for accumulated error. Defaults to 0. Passed to \code{integrate()}.}

\item{upper}{Upper boundary for accumulated error. Defaults to 1. Passed to \code{integrate()}.}

\item{subdivisions}{The maximum number of subintervals. Passed to \code{integrate()}.}

\item{abs.tol}{Absolute accuracy. Passed to \code{integrate()}}

\item{rel.tol}{Relative accuracy. Passed to \code{integrate()}}
}

\details{
The function manufactured by \code{covDiagFactory} has the form \code{function(target_theta_list, target_args_list = list(), reference_cov.theta = NULL, ...)} and serves to iterate over a large variety of parameters.

\describe{
\item{target_theta_list}{List of named numeric vectors, each providing at least one or more values for a parameter.}

\item{target_args_list (default = \code{list()})}{List of named vectors, each providing at least one or more values for each argument to be provided to target_covariance.}

\item{reference_cov.theta (default = \code{NULL})}{A numeric vector with parameters for the reference covariance. This is primarily intended to use when comparing different covariance functions, in which case the target covariance is compared to the reference with constant parameters.If this argument is left at default and the target and reference are identical, the parameters in \code{target_theta_list} will be used by both. Otherwise this will raise an exception.}
}
}

\section{Available diagnostic functions}{
All diagnostic functions have in common that, apart from the required parameters kappa and mu, any remaining parameters are optional and captured using the \code{...} operator. Any unspecified arguments are left at default. The provided arguments are then turned into a grid using \code{expand.grid}, which is then iterated over to compute the requested diagnostic metrics.
Note that there are some built-in checks which ensure that invalid configurations are dropped from the grid. For example, rows with interpolator = 'none' and a non-zero number of supports are dropped to prevent redundant computations.

\subsection{Accumulated Error}{The accumulated error corresponds to the area between the approximated and exact covariance curve.}

\subsection{Point diagnostics}{This option returns error metrics across a user-defined number of points on the interval [0,1]. The metrics returned are "error", "absolute error", "maximum error", "target covariance value", and "reference covariance value".}
}

\author{Thomas Caspar Fischer}

\seealso{\code{\link{cov.wendland}}}

\examples{

interpolators <- c("linear", "cspline", "polynomial")
diagnostics <- c("accumulated_error", "point_diagnostics")

diagnosticFun <- covDiagFactory(cov.wendland, diagnostic_funs = diagnostics,
  reference_covariance = cov.wendland)
target_theta_list <- list(range = 0.5, sill = 1, kappa = c(0, 0.5, 1), mu = 0,
  nugget = 0)
target_args_list <- list(interp.method = interpolators, interp.num_support = 25)
wendland_comparison <- diagnosticFun(target_theta_list = target_theta_list,
  target_args_list = target_args_list)

diagnosticFun <- covDiagFactory(cov.wendland, diagnostic_funs = diagnostics,
  reference_covariance = cov.askey)
target_theta_list <- list(range = 0.5, sill = 1, kappa = 0, mu = 0.5, nugget = 0)
target_args_list <- list(interp.method = interpolators, interp.num_support = 25)
askey_comparison <- diagnosticFun(target_theta_list = target_theta_list,
  target_args_list = target_args_list, reference_cov.theta = c(0.5, 1, 0.5, 0))
}

\concept{wendland}