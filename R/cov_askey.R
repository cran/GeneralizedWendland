


cov.askey <- function(h, theta, ..., cov.args = list()){

  if (length(theta) != 4) stop("Invalid number of parameters")

  cov.args <- do.call(control.askey, cov.args)

  range  <- theta[1]
  sill   <- theta[2]
  mu     <- theta[3] + (cov.args[["cov.reparameterize"]]) * ((1+cov.args[["cov.d_value"]])/2)
  nugget <- theta[4]

  h <- h/range

  if (is.spam(h)) {
    h@entries <- ifelse(h@entries < 1,(1-h@entries)^mu,0) + nugget*(h@entries < cov.args[["cov.eps"]])
  } else {
    h <- sill*(((1-h)*(h<1))^mu) + nugget*(h<cov.args[["cov.eps"]])
  }

  return(h)
}
