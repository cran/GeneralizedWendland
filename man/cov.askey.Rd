\name{cov.askey}
\alias{cov.askey}

\title{Askey covariance function}

\description{
A covariance function of the form

\deqn{\rho_{\beta,\kappa,\mu} = \begin{cases} \sigma + \theta \quad 0 \leq r < \epsilon \\ (1-r)^{\mu} \quad  \epsilon \leq r < 1 \\ 0 \quad 1 \leq r \end{cases}}
where \eqn{r=h/\beta}. This is equivalent to the generalized Wendland covariance with \eqn{\kappa=0}, but much more computationally efficient.
}

\usage{
cov.askey(h, theta, ..., cov.args = list())
}

\arguments{
\item{h}{A numeric vector, matrix, or spam object storing distances.}

\item{theta}{Numeric vector \eqn{\vec{\theta}=(\beta, \sigma, \mu, \theta)} storing parameters.}

\item{...}{Other arguments.}

\item{cov.args}{Named list of arguments. See Details.}
}

\value{
Returns an object of the same type as input object h which stores the computed covariance values, i.e. a spam object if input h was also a spam object.
}

\details{
Using the list cov.args, users can provide the following arguments:

\describe{
\item{cov.eps (default: \code{.Machine$double.eps^0.5})}{The threshold distance \eqn{\epsilon} below which the function will return \eqn{\sigma + \theta}.}
}
}


\references{
Moreno Bevilacqua and Tarik Faouzi and Reinhard Furrer and Emilio Porcu (2019) Estimation and prediction using generalized Wendland covariance functions under fixed domain asymptotics, \emph{Annals of Statistics}, \bold{47}(2), 828--856.
}

\author{Thomas Caspar Fischer}


\seealso{
\code{\link{cov.wendland}}
}
