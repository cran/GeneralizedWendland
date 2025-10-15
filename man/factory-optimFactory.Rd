\name{optimFactory}

\alias{optimFactory}

\title{Function Factory for Optimization function with Unified Arguments}

\description{
  A function factory which returns a function with unified input arguments,
  and provides compatibility with the package \pkg{optimParallel}.
}

\usage{
  optimFactory(optim.args = list(), hessian = FALSE, optimParallel.args = list())
}

\value{
  A function of the form \code{function(par, fn, gr = NULL, ..., lower, upper)}
  which returns the output obtained from calls to \code{optim} or \code{optimParallel}
}

\arguments{

  \item{optim.args}{
    A list of optional settings for optim. See \code{\link{optim}} for documentation
    of valid arguments.
  }

  \item{hessian}{
    A logical which specifies whether the hessian matrix is to be returned with the output.
  }

  \item{optimParallel.args}{
    A list which is passed to \code{optimParallel} from the \pkg{optimParallel} package. See
    details.
  }
}

\details{
  The function factory \code{optimFactory()} returns a function of the form
  \code{function(par, fn, gr = NULL, ..., lower, upper)}. It is intended to replace
  calls to \code{optim} or \code{optimParallel} by wrapping both functions.
  By default, it returns a function that corresponds to \code{optim} with default
  arguments.

  To use optimParallel, users may specify the following arguments in optimParallel.args:
  \describe{
    \item{num_cores (default = \code{NULL})}{
      The number of cores to use during numerical optimization. Is NULL by default, which
      corresponds to using \code{stats::optim}. When \code{num_cores} is a numeric value,
      the actual number of cores is set to \code{min(detectCores()-1, num_cores)} to avoid
      accidentally overloading the user's system.
    }

    \item{forward (default = \code{FALSE})}{
      A logical value which controls whether \code{optimParallel} should use central difference
      approximation of the gradient (\code{FALSE}) or forward difference approximation (\code{TRUE}).
    }

    \item{loginfo (default = \code{FALSE})}{
      A logical value which controls whether \code{optimParallel} should return
      additional information about the optimization process. See \code{\link[optimParallel]{optimParallel}}.
    }
  }
}

\note{
  When using the \pkg{optimParallel} package, there are some caveats to consider
  depending on the user's operating system. On Windows systems, the function
  needs to be provided with a character vector with all required packages. Under
  Linux, this should not be required, but the functionality is currently bugged.
}

\references{
  Hadley Wickham (2015) \emph{Advanced R}, CRC Press.

  Florian Gerber and Reinhard Furrer (2019) optimParallel: An R package
  providing a parallel version of the L-BFGS-B optimization method,
  \emph{The R Journal}, \bold{11}(1), 352--358
}

\author{
  Thomas Caspar Fischer
}

\seealso{
  \code{\link{optim}} and
  \code{\link[optimParallel]{optimParallel}}
}

\examples{
library(GeneralizedWendland)
library(optimParallel)

set.seed(43)
n <- 50
range <- 0.4
dist_max <- 2
theta  <- c(range, 1, 1, 0, 0)

locs <- data.frame(x = runif(n, 0, sqrt(dist_max)),
                   y = runif(n, 0, sqrt(dist_max)))
dmat  <- spam::nearest.dist(locs, locs, delta = dist_max)
Sigma <- cov.wendland(h = dmat, theta = theta)
y <- c(spam::rmvnorm(1, Sigma = Sigma))

init_parameters   <- c(0.7, 2, 0, 2, 2)
lower_constraints <- c(0.1, 0.1, 0, 0, 0)
upper_constraints <- c(sqrt(2), 2, 2, 2, 2)

mleFunction <- mleFactory(covariance = cov.wendland)
(mle_result <- mleFunction(y = y, distmat = dmat, init_parameters = init_parameters,
                       theta_llim = lower_constraints, theta_ulim = upper_constraints))

mleFunctionPar <- mleFactory(covariance = cov.wendland, optimParallel.args = list(num_cores = 2))
(mle_result_par <- mleFunctionPar(y = y, distmat = dmat, init_parameters = init_parameters,
                       theta_llim = lower_constraints, theta_ulim = upper_constraints))

}

\concept{factory}
