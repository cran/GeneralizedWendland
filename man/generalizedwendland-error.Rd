\name{generalizedwendland-error}
\alias{generalizedwendland-error}
\alias{update.spam.chol.error}
\alias{covfun.notfunction.error}


\title{Simple errors defined by \pkg{GeneralizedWendland} package}

\description{
An overview of the \code{simpleError} objects defined by the \pkg{GeneralizedWendland} package.
}

\section{Error definitions}{
\describe{

\item{\code{update.spam.chol.error}: \emph{"Updated covariance entries do not match length of original one."}}{ Warning which is raised when \code{update.spam.chol.NgPeyton} fails because the length of the entries in the updated covariance matrix is not equal to the length of the original matrix. This usually indicates that the initial sparsity structure was too sparse for the given parameters.}

\item{\code{wendland.insuffparam.error}: \emph{"Too few parameters for Wendland."}}{Insert description here}

\item{\code{wendland.excessparam.error}: \emph{"Too many parameters for Wendland. Did you supply fix range or nugget?"}}{Insert description here}


\item{\code{covfun.notfunction.error}: \emph{"Argument covariance must be a function."}}{Insert description here}

}
}

\author{Thomas Caspar Fischer}
