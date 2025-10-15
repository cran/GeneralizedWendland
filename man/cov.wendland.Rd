\name{cov.wendland}
\alias{cov.wendland}

\title{Generalized Wendland Covariance Function}

\description{
A fully parametrized generalized Wendland covariance function for use in geostatistical modeling, as well as multiple methods of obtaining computationally inexpensive approximations.

\deqn{\rho_{\beta,\kappa,\mu} = \begin{cases}\sigma + \theta \quad 0 \leq r < \epsilon \\ \frac{\sigma}{B(1+2\kappa,\mu)} \int_r^1 (u^2-r^2)^\kappa (1-u)^{\mu-1} du \quad \epsilon \leq r < 1 \\ 0 \quad 1 \leq r \end{cases}}
where \eqn{r=h/\beta}
}

\usage{
cov.wendland(h, theta, ..., cov.args = list())
}

\arguments{

\item{h}{A numeric vector, matrix, or spam object storing distances.}

\item{theta}{Numeric vector \eqn{\vec{\theta}=(\beta, \sigma, \kappa, \mu, \theta)} storing parameters.}

\item{\dots}{Other arguments.}

\item{cov.args}{Named list of arguments. See Details.}
}

\value{
Returns an object of the same type as input object h which stores the computed covariance values, i.e. a spam object if input h was also a spam object.
}

\details{
Using the list cov.args, users can provide the following additional arguments:

\describe{
\item{numint.abstol (default: \code{1e-3})}{Absolute tolerance for numerical integration.}

\item{numint.reltol (default: \code{1e-3})}{Relative tolerance for numerical integration.}

\item{numint.qag_key (default: \code{0})}{Method to use in QAG integration (Values 1 - 6)}

\item{numint.subintervals (default: \code{0})}{Number of subintervals to use in QAG/QAGS integration.}

\item{interp.method (default: \code{'none'})}{Method to use for covariance interpolation. Valid methods are 'none', 'linear', 'polynomial', and 'cspline'.}

\item{interp.num_support (default: \code{0})}{Number of support points to use for covariance interpolation.}

\item{cov.reparameterize (default: \code{TRUE})}{Whether to apply the reparameterization \eqn{\mu=\frac{1+d}{2} + \kappa + \nu}, where \eqn{\nu} takes the place of \eqn{\mu} in input vector \eqn{\vec{\theta}}. This allows users to use box constraints in maximum likelihood estimation, as the covariance function is valid for \eqn{\nu \in [0,\infty)} rather than \eqn{\mu \in [\frac{1+d}{2} + \kappa,\infty)}.}

\item{cov.eps (default: \code{.Machine$double.eps^0.5})}{The threshold distance \eqn{\epsilon} below which the function will return \eqn{\sigma + \theta}.}

\item{cov.d_value (default: \code{2})}{Dimensionality of space in which measurements were taken. This only takes effect if \code{cov.reparameterize} is TRUE.}

}
}

\references{
Moreno Bevilacqua and Tarik Faouzi and Reinhard Furrer and Emilio Porcu (2019) Estimation and prediction using generalized Wendland covariance functions under fixed domain asymptotics, \emph{Annals of Statistics}, \bold{47}(2), 828--856.
}

\author{Thomas Caspar Fischer}

\examples{
h <- seq(0, 1, 0.01)
plot(0, type = "n", xlab = "Distance", ylab = "Covariance",
     xlim = c(0, 1), ylim = c(0,1))

theta <- c(range=1, sill=1, kappa=1, mu=0, nugget=0)
cov.args <- list()
lines(x = h, y = cov.wendland(h, theta, cov.args = cov.args),
      lwd = 2)

theta <- c(range=1, sill=1, kappa=1, mu=0, nugget=0)
cov.args <- list(cov.reparametrize = FALSE, cov.d_value = 2)
theta[4] <- (1 + cov.args[[2]])/2 + theta[3] + theta[4]
lines(x = h, y = cov.wendland(h, theta, cov.args = cov.args),
      col = "red", lty = 3, lwd = 3.5)

theta <- c(range=0.5, sill=1, kappa=1, mu=0, nugget=0)
cov.args <- list(interp.method="cspline", interp.num_support=100)
lines(x = h, y = cov.wendland(h, theta, cov.args = cov.args),
      col = "green", lwd = 2)

legend("topright", legend = c("Default", "No reparameterization",
                              "Cubic spline interpolation"),
       col = c(1, 2, 3), lty = c(1,3,1), lwd = c(2, 3.5, 2))
}

\concept{wendland}
