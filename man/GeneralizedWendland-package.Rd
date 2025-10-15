\name{GeneralizedWendland-package}

\alias{GeneralizedWendland-package}
\alias{GeneralizedWendland}
\alias{generalizedwendland-package}
\alias{generalizedwendland}

\docType{package}

\title{Fully Parametrized Generalized Covariance Function for R}

\description{
This package provides a fully parametrized generalized Wendland covariance function for use in geostatistical modeling, as well as various options for approximations and adjustments. In addition, the package comes with a modified process for parameter estimation, based on the spam-implementations, but compartmentalized using function factories.
}

\details{
The package requires the user to complete some initial steps before it can be used. Apart from the R dependencies, the user will need to manually install the GNU Scientific Library on their system and ensure that R can find the required libraries. On Windows systems, users can install devtools and then install gsl via MSYS2 > pacman.
}

\seealso{
\code{\link{cov.wendland}}.
}

\author{Thomas Caspar Fischer}

\keyword{documentation}
\keyword{package}
\concept{wendland}


