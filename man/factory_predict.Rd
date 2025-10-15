\name{factory-predict}

\alias{predictionFactory}

\title{Function factory for prediction generator}

\description{
A function factory which generates a function for sampling from the predictive distribution
}

\usage{
predictionFactory(y, locs0, locs1, covariance, X0 = list2DF(nrow = nrow(locs0)),
  X1 = list2DF(nrow = nrow(locs1)), ..., cov.args = list(), chol.args = list(),
  use_spam = TRUE)
}

\arguments{
\item{y}{Observed measurements.}

\item{locs0}{Locations at which measurements were obtained.}

\item{locs1}{Locations at which to predict.}

\item{covariance}{Covariance function.}

\item{X0}{Covariates. By default corresponds to an empty data.frame() with \code{nrow(locs0)} rows.}

\item{X1}{Covariates for new locations. By default corresponds to an empty data.frame() with \code{nrow(locs1)} rows.}

\item{...}{Additional arguments. Unused.}

\item{cov.args}{Additional arguments for covariance function.}

\item{chol.args}{Additional arguments for function used to compute cholesky decomposition.}

\item{use_spam}{Logical value. If TRUE, use spam to compute predictions.}
}

\value{
Returns a function of the form \code{function(n, param)}, where n corresponds to the number of simulations to run, and param is a parameter vector as returned by \code{optim}.
}

\author{Thomas Caspar Fischer}

\examples{
set.seed(23)
n <- 100
res <- c(20, 20)

locs <- data.frame(x = runif(n), y = runif(n))
locs_new <- expand.grid(x = seq(0, 1, length.out = res[1]),
                        y = seq(0, 1, length.out = res[2]))
range <- 0.3
dmat <- as.matrix(dist(locs))
theta <- c(range, 1, 1, 0, 0)
cov.args <- list()
chol.args <- list()
Sigma <- cov.wendland(h = dmat, theta = theta, cov.args = cov.args)
y <- c(spam::rmvnorm(n = 1, Sigma = Sigma))

predictionFun <- predictionFactory(y = y, locs0 = locs, locs1 = locs_new,
  covariance = cov.wendland, cov.args = cov.args, chol.args = chol.args,
  use_spam = FALSE)

predictions <- predictionFun(n=10, param=theta)

image(x = seq(0, 1, length.out=res[1]),
      y = seq(0, 1, length.out=res[2]),
      z = matrix(apply(predictions,2,mean), res[1], res[2]),
      col = hcl.colors(9, "Blue-Red"),
      breaks = qnorm(seq(0.05, 0.95, 0.1)),
      xlab = "", ylab = "", xaxt = "n", yaxt = "n",
      useRaster = TRUE)

}

\concept{wendland}