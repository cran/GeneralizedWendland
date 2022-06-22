\name{covarianceFactory}

\alias{covarianceFactory}

\title{covarianceFactory}

\description{
A function factory which sets up a covariance function.
}

\usage{
covarianceFactory(covariance, cov.args = list())
}

\arguments{
\item{covariance}{A function which takes as input an object containing distances (h), a vector of parameters (theta), and a list of optional settings (cov.args)}

\item{cov.args}{A list of optional settings for a covariance function.}
}

\value{
A function of the form \code{function(h, theta)}. This function is enclosed in the execution environment of \code{choleskyFactory} and hence has access to the arguments \code{covariance} and \code{cov.args}. The manufactured function returns the result of \code{covariance(h = h, theta = theta, cov.args = cov.args)}.
}

\details{
\code{covarianceFactory()} is a function factory which takes as input an arbitrary covariance function and a list of additional arguments and returns a function with these arguments. The argument \code{cov.args} mainly serves to pass the relevant arguments to the covariance function, but also allows users to specify the following two arguments:

\describe{
\item{cov.fixed_range_value (default: \code{NA})}{A constant value for the range parameter. Note that the code assumes that the first value of \code{theta} corresponds to the range parameter, as this notation is used in the \pkg{spam} package.}

\item{cov.fixed_nugget_value (default: \code{NA})}{A constant value for the nugget parameter. Note that the code assumes that the last value of \code{theta} corresponds to the nugget parameter, as this notation is used in the \pkg{spam} package.}
}

Note that if either of these arguments are specified, the corresponding entries should also be omitted in \code{theta}.
}

\references{
Wickham, H. (2015) \emph{Advanced R}, CRC Press.
}

\author{Thomas Caspar Fischer}

\seealso{
\code{\link[spam]{covmat}} and
\code{\link{cov.wendland}}
}

\examples{
library(GeneralizedWendland)
library(spam)

hs <- seq(0, 1, 0.05)
covSph <- covarianceFactory(cov.sph)
covExp <- covarianceFactory(cov.exp)
covWend <- covarianceFactory(cov.wendland)

plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1))
lines(hs, covSph(hs, theta = c(0.5, 1, 0)))
lines(hs, covExp(hs, theta = c(0.5, 1, 0)), col = "red")
lines(hs, covWend(hs, theta = c(0.5, 1, 0.5, 0.5, 0)), col = "green")
}

\concept{factory}