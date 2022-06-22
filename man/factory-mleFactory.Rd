\name{mleFactory}

\alias{mleFactory}

\title{Function Factory for Generating mle Function with unified arguments}

\description{
A factory function which returns a function of the form \code{function(y, X = data.frame(), distmat, init_parameters, theta_llim, theta_ulim)} which can be called to compute the maximum likelihood estimates for a Kriging model.
}

\usage{
mleFactory(covariance, cov.args = list(), chol.args = list(),
  optim.args = list(), hessian = FALSE, optimParallel.args = list())
}

\value{
A function of the form \code{function(y, X = data.frame(), distmat, init_parameters, theta_llim, theta_ulim)}. The manufactured function itself has the form \code{function(y, X = data.frame(), distmat, init_parameters, theta_llim, theta_ulim} and returns the output of \code{optim} or \code{optimParallel} when \code{optimParallel.args} was specified.
}

\arguments{
\item{covariance}{A function of the form \code{function(h, theta, ..., cov.args = list())}, where h is an object storing a distance matrix, theta is a numeric vector of parameters for the linear predictor and covariance function, and cov.args is a list of optional arguments for the covariance function.}

\item{cov.args}{A list of optional settings for a covariance function.}

\item{chol.args}{A list of optional settings for a cholesky function. \emph{Note}: Valid input arguments change depending on whether the distance matrix provided to the output function is sparse. This may change in a future version.}

\item{optim.args}{A list of optional settings for optim. See \code{\link{optim}} for documentation of valid arguments.}

\item{hessian}{A logical value which specifies whether the hessian matrix is to be returned in the output. Is FALSE by default.}

\item{optimParallel.args}{A list of optional settings for optimParallel. See \code{\link{optimParallel}} for documentation of valid arguments.}
}

\details{
The purpose of this function factory is to return an mle function with unified arguments. The returned function performs the same task as for example \code{spam::mle()}, but simplifies the process in two ways: The returned function detects whether the Gaussian process is a zero-mean process through the input argument X and whether methods from the \pkg{spam} package should be used based on the type of input argument distmat, and autonomously chooses appropriate methods to compute the neg2loglikelihood. Hence the user does not need to choose a specialized method themselves.
}

\references{
Hadley Wickham (2015) \emph{Advanced R}, CRC Press.
}

\author{Thomas Caspar Fischer}

\seealso{
\code{\link{optim}},
\code{\link{optimParallel}},
\code{\link{covarianceFactory}},
\code{\link{choleskyFactory}} and
\code{\link{optimFactory}}
}

\examples{
set.seed(57)
n <- 50
range <- 0.4
theta  <- c(range, 1, 1, 0, 0)

locs <- data.frame(x = runif(n), y = runif(n))
dmat  <- as.matrix(dist(locs))
Sigma <- cov.wendland(h = dmat, theta = theta)
y <- c(spam::rmvnorm(1, Sigma = Sigma))

init_parameters   <- c(0.7, 2, 0, 2, 2)
lower_constraints <- c(0.1, 0.1, 0, 0, 0)
upper_constraints <- c(sqrt(2), 2, 2, 2, 2)

mleFunction <- mleFactory(covariance = cov.wendland)
mle_result1 <- mleFunction(y = y, distmat = dmat,
                           init_parameters = init_parameters, theta_llim = lower_constraints,
                           theta_ulim = upper_constraints)

mleFunctionDM <- mleFactory(covariance = cov.wendland,
                            cov.args = list(fixed_range_value = range))
mle_result2 <- mleFunctionDM(y = y, X = data.frame(), distmat = dmat,
                             init_parameters = init_parameters[-1],
                             theta_llim = lower_constraints[-1],
                             theta_ulim = upper_constraints[-1])
}

\concept{factory}
\concept{wendland}