\name{generalizedwendland-warning}
\alias{generalizedwendland-warning}
\alias{update.spam.chol.warn}
\alias{wendland.interp.redundantsupport.warn}
\alias{wendland.interp.lowsupport.warn}
\alias{wendland.interp.unimplemented.warn}

\title{Simple warnings defined by \pkg{GeneralizedWendland} package.}

\description{
An overview of the \code{simpleError} objects defined by the \pkg{GeneralizedWendland} package.
}

\section{Warning definitions}{

\describe{

\item{\code{update.spam.chol.warn}: \emph{"Updated covariance entries do not match length of original one. Deleting stored Rstruct."}}{
Warning which is raised when \code{update.spam.chol.NgPeyton} fails because the length of the entries in the updated covariance matrix is not equal to the length of the original matrix. This usually indicates that the initial sparsity structure was too sparse for the given parameters.}

\item{\code{wendland.interp.redundantsupport.warn}: \emph{"Argument interp.num_support > 0 while using exact method. Set to 0."}}{Insert description here.}

\item{\code{wendland.interp.lowsupport.warn}: \emph{"Argument interp.method != 'none' with less than 3 support points. Forced to 'none'."}}{Insert description here.}

\item{\code{wendland.interp.unimplemented.warn}: "Interpolator not implemented. Forcing exact method."}{Insert description here.}

}
}

\author{Thomas Caspar Fischer}