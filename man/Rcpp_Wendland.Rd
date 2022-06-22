
\name{Rcpp_Wendland-class}
\docType{class}
\alias{Rcpp_Wendland-class}

\title{Class "Rcpp_Wendland"}

\description{
Rcpp class which serves as an interface to the C++ implementation of the generalized Wendland covariance function.
}

\details{
Parts of the covariance function require C and C++ dependencies. This Rcpp class serves as an interface to these dependencies.
}

\section{Methods}{

\subsection{Constructor}{
\describe{
\item{Constructor}{\code{wend <- new("Rcpp_Wendland")}: Creates an instance of the Rcpp_Wendland class.}
}
}

\subsection{Get/Set methods}{
\describe{
\item{setParameters}{\code{wend$setParameters(range, sill, kappa, mu, nugget)}: Set Parameters}

\item{setEpsTol}{\code{wend$setEpsTol(eps)}: set numeric precision.}
}
}

\subsection{Integration options}{
\describe{
\item{setIntegratorQNG}{\code{wend$setIntegrator(abstol, reltol, intervals=0, qag_key=0)}}
}
}

\subsection{Interpolation options}{
\describe{
\item{setInterpolator}{\code{wend$setInterpolatorLinear(num_points, interp_type=0)}}
}
}

\subsection{Computation}{
\describe{
\item{compute}{\code{wend$compute(d)}: Compute for single value.}

\item{computeVector}{\code{wend$computeVector(d)}: Compute for numeric vector.}

\item{computeMatrix}{\code{wend$computeMatrix(d)}: Compute for numeric matrix.}

\item{computeMSparse}{\code{wend$computeMSparse(d)}: Compute for dgCMatrix from \pkg{Matrix}.}

\item{computeSpam}{\code{wend$computeSpam(index, values)}: Compute for spam object, turned into triplet form using \code{spam::triplet}.}
}
}
}

\references{
Dirk Eddelbuettel (2013) Seamless R and C++ Integration with Rcpp. \emph{Springer}, New York
}

\author{Thomas Caspar Fischer}

\seealso{
\code{\link{cov.wendland}}
}

\keyword{classes}
\concept{wendland}