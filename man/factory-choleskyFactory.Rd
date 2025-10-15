\name{choleskyFactory}
\alias{choleskyFactory}

\title{Function Factory for Generating chol Function with unified arguments}

\description{
A function factory which returns a function of the form \code{function(Sigma)} which performs a cholesky decomposition using an approach tailored to the type of input Sigma. Currently works for vectors, matrices, spam objects, and dgCMatrix objects from the \pkg{Matrix} package.
}

\usage{
choleskyFactory(chol.args = list(), Rstruct = NULL)
}

\arguments{
\item{chol.args}{A list of optional settings for a cholesky function.}

\item{Rstruct}{A \code{spam.chol.NgPeyton} object which represents the sparsity structure of covariance matrix Sigma.}
}

\value{
A function of the form \code{function(Sigma)}.
}

\details{
The output of \code{choleskyFactory} is intended to replace calls to \code{chol.default} or \code{chol.spam}. The object type is determined during runtime, after which the appropriate function is called to obtain the cholesky decomposition. For spam objects, the function attempts to use \code{update.spam.chol.NgPeyton} if \code{Rstruct} is specified, and upon failure defaults to \code{chol.spam}. The result is then assigned in the execution environment of \code{choleskyFactory}, so that \code{Rstruct} will be defined in the next call.
}

\references{
Hadley Wickham (2015) \emph{Advanced R}, CRC Press.
Reinhard Furrer and Roman Flury and Florian Gerber (2022) spam: SPArse Matrix, R package version 2.8-0.
}

\author{Thomas Caspar Fischer}

\seealso{
\code{\link[base]{chol}},
\code{\link[spam]{chol.spam}},
\code{\link[spam]{update.spam.chol.NgPeyton}}
}

\examples{
set.seed(1234)
locations <- data.frame(x = runif(10), y = runif(10))
theta <- c(0.5,1,1,0,0)
dmat <- as.matrix(dist(locations))
Sigma <- cov.wendland(dmat, theta)

cholFun <- choleskyFactory(chol.args = list())
cholD <- cholFun(Sigma)

cholFun <- choleskyFactory(chol.args = list(pivot = TRUE))
cholD_pivot <- cholFun(Sigma)

cholFun <- choleskyFactory(chol.args = list(pivot = "RCM"))
cholS_RCM <- cholFun(spam::as.spam(Sigma))
}

\concept{factory}