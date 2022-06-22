\name{neg2loglikDiagFactory}

\alias{neg2loglikDiagFactory}

\title{Diagnostics for arbitrarily specified, likelihood-based Gaussian process models}

\description{
A helper function for rapidly exploring the parameter space around the maximum likelihood estimate
}

\usage{
neg2loglikDiagFactory(y, X = data.frame(), distmat, covariance, ...)
}

\value{
Returns a function of the form \code{function(theta_list, cov.args_list = list(), chol.args_list = list())} which returns a \code{data.frame} containing the neg2loglikelihood at all permutations of the provided arguments.
}

\arguments{
\item{y}{Dependent variable}

\item{X}{Optional design matrix with covariates}

\item{distmat}{Distance matrix. Can be provided either as a dense matrix or spam object.}

\item{covariance}{Covariance function.}

\item{...}{Other arguments to be passed on.}
}

\details{

\describe{
\item{theta_list}{Named list of vectors with parameters to be passed to covariance.}

\item{cov.args_list (default = \code{list()})}{Named list of vectors with arguments to be passed to covariance}

\item{chol.args_list (default = \code{list()})}{Named list of vectors with arguments to be passed to \code{\link{choleskyFactory}}.}

}
}

\note{
The function manufactured by \code{neg2loglikDiagFactory} in principle also accepts covariance functions generated using \code{\link{covarianceFactory}}. However, the function is not yet compatible with the arguments \code{fixed_range_value} and \code{fixed_nugget_value}. For now, these should be left at default when using \code{\link{covarianceFactory}}.
}

\author{Thomas Caspar Fischer}

\examples{
set.seed(63)
n <- 50
range <- 0.7
theta  <- c(range, 1, 1, 0, 0)

locs <- data.frame(x = runif(n), y = runif(n))
dmat  <- as.matrix(dist(locs))
Sigma <- cov.wendland(h = dmat, theta = theta)
y <- c(spam::rmvnorm(1, Sigma = Sigma))

neg2loglikIterator <- neg2loglikDiagFactory(y = y, distmat = dmat,
  covariance = cov.wendland)
theta_list <- list(range = 0.5, sill = 1, kappa = 0, mu = c(0, 0.25, 0.5),
  nugget = 0)
cov.args_list <- list(numint.abstol = c(1e-1, 1e-3, 1e-6), numint.reltol = c(1e-3))

results <- neg2loglikIterator(theta_list, cov.args_list = cov.args_list)
}
