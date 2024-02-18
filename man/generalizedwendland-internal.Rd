\name{generalizedwendland-internal}
\alias{generalizedwendland-internal}
\alias{control.wendland}
\alias{control.askey}

\title{Internal Functions defined in package \pkg{GeneralizedWendland}}

\description{
The functions or variables listed here are used internally, i.e. they are not intended to be called directly by the user.
}

\usage{
control.wendland(numint.abstol = getOption("wendland.numint.abstol"),
  numint.reltol = getOption("wendland.numint.reltol"),
  numint.qag_key = getOption("wendland.numint.qag_key"),
  numint.subintervals = getOption("wendland.numint.subintervals"),
  interp.method = getOption("wendland.interp.method"),
  interp.num_support = getOption("wendland.interp.num_support"),
  cov.reparameterize = getOption("wendland.cov.reparameterize"),
  cov.eps = getOption("wendland.cov.eps"),
  cov.d_value = getOption("wendland.cov.d_value"), ...)

control.askey(cov.reparameterize = getOption("wendland.cov.reparameterize"),
  cov.eps = getOption("wendland.cov.eps"),
  cov.d_value = getOption("wendland.cov.d_value"), ...)
}

\arguments{
\item{numint.abstol}{Absolute tolerance for numerical integration (default: 1e-3).}

\item{numint.reltol}{Relative tolerance for numerical integration (default: 1e-3).}

\item{numint.qag_key}{Method to use in QAG integration, from 0-6 (default: 0).}

\item{numint.subintervals}{Number of subintervals to use in QAG/QAGS integration (default: 0).}

\item{interp.method}{Method to use for covariance interpolation. Valid methods are 'none', 'linear', 'polynomial', and 'cspline' (default: 'none').}

\item{interp.num_support}{Number of support points to use for covariance interpolation (default: 0).}

\item{cov.reparameterize}{Whether to apply the reparameterization \eqn{\mu=\frac{1+d}{2} + \kappa + \nu}, where \eqn{\nu} takes the place of \eqn{\mu} in input vector \eqn{\vec{\theta}} (default: TRUE). This allows users to use box constraints in maximum likelihood estimation, as the covariance function is valid for \eqn{\nu \in [0,\infty)} rather than \eqn{\mu \in [\frac{1+d}{2} + \kappa,\infty)}.}

\item{cov.eps}{The threshold distance \eqn{\epsilon} below which the function will return \eqn{\sigma + \theta} (default: .Machine$double.eps^0.5).}

\item{cov.d_value}{Dimensionality of space in which measurements were taken (default: 2). This only takes effect if \code{cov.reparameterize} is TRUE.}

\item{...}{Other arguments. Unused.}
}

\value{
A named list.
}

\details{
The functions \code{control.wendland} and \code{control.askey} are helper functions which create a named list containing the relevant default arguments passed via the \code{cov.args} argument of \code{cov.wendland} and \code{cov.askey} respectively. In addition, users can override the default values.
}

\author{Thomas Caspar Fischer}

\keyword{internal}