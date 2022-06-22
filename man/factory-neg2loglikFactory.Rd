\name{neg2loglikFactory}

\alias{neg2loglikFactory}

\title{Function Factory for \code{neg2loglikelihood} with Unified Input Arguments}

\description{
A function factory which generates a function of the form \code{function(parameters)} which returns the neg2loglikelihood.
}

\usage{
neg2loglikFactory(y, X, distmat, covariance = NULL, cov.args = list(),
  chol.args = list(), Rstruct = NULL, covarianceFunction = NULL,
  choleskyFunction = NULL)
}

\arguments{
\item{y}{Numeric vector. Dependent variable.}

\item{X}{Optional data.frame containing covariates.}

\item{distmat}{Distance matrix, either a numeric matrix or a spam object.}

\item{covariance}{A function which takes as input an object containing distances (h), a vector of parameters (theta), and a list of optional settings (cov.args).}

\item{cov.args}{A list of optional settings for a covariance function.}

\item{chol.args}{A list of optional settings for a cholesky function.}

\item{Rstruct}{A 'spam.chol.NgPeyton' object which represents the sparsity structure.}

\item{covarianceFunction}{A function returned by \code{covarianceFactory()}.}

\item{choleskyFunction}{A function returned by \code{choleskyFactory()}.}
}

\value{
Returns function of the form \code{function(parameters)}.
}

\details{
This function factory returns a function of the form \code{function(parameters)} which computes the neg2loglikelihood for given input parameters. The purpose of this is to reduce the number of arguments that need to be specified by the user in a call to \code{optim}, or \code{optimParallel}. Furthermore, the function detects whether the input distmat is a spam object, and autonomously selects the appropriate method for computing the neg2loglikelihood.

The function is intended to be called from within \code{\link{mleFactory}}, but is also exported by NAMESPACE for users wishing to make use of the function. There are two distinct strategies available for using the function.

\bold{Option 1}: the user may specify covariance, cov.args, chol.args, and Rstruct in the call. This syntax is more in line with the corresponding functions found in the \pkg{spam} package, yet still allows passing arguments for customizing the behaviour of the cholesky decomposition.

\bold{Option 2}: the user may instead specify covarianceFunction and choleskyFunction, obtained from calls to \code{covarianceFactory} and \code{choleskyFactory}, respectively.

In both cases, the arguments y, X, and distmat are required input. Note that the two options are equivalent, apart from the second option allowing for more concise code.
}

\references{
Hadley Wickham (2015) \emph{Advanced R}, CRC Press.
}

\author{Thomas Caspar Fischer}

\seealso{
\code{\link{covarianceFactory}} and
\code{\link{choleskyFactory}}
}

\examples{
set.seed(63)
n <- 50
range <- 0.7
theta  <- c(range, 1, 1, 0, 0)

locs <- data.frame(x = runif(n), y = runif(n))
dmat  <- as.matrix(dist(locs))
Sigma <- cov.wendland(h = dmat, theta = theta)
y <- c(spam::rmvnorm(1, Sigma = Sigma))
X <- data.frame()

neg2loglikFun <- neg2loglikFactory(y = y, X = X, distmat = dmat,
  covariance = cov.wendland, cov.args = list(), chol.args = list())
result1 <- neg2loglikFun(theta)

covarianceFun <- covarianceFactory(cov.wendland, cov.args = list())
choleskyFun <- choleskyFactory(chol.args = list())
neg2loglikFun <- neg2loglikFactory(y = y, X = X, distmat = dmat,
  covarianceFunction = covarianceFun, choleskyFunction = choleskyFun)
result2 <- neg2loglikFun(theta)
}

\concept{factory}