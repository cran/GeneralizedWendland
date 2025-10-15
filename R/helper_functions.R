
# Function which assigns global values to cov.args in calls to cov.wendland
control.wendland <- function(
  numint.abstol          = getOption("wendland.numint.abstol"),
  numint.reltol          = getOption("wendland.numint.reltol"),
  numint.qag_key         = getOption("wendland.numint.qag_key"),
  numint.subintervals    = getOption("wendland.numint.subintervals"),
  interp.method          = getOption("wendland.interp.method"),
  interp.num_support     = getOption("wendland.interp.num_support"),
  cov.reparameterize     = getOption("wendland.cov.reparameterize"),
  cov.eps                = getOption("wendland.cov.eps"),
  cov.d_value            = getOption("wendland.cov.d_value"),
  ...
){
  cov.args <- list( numint.abstol          = numint.abstol,
                    numint.reltol          = numint.reltol,
                    numint.qag_key         = numint.qag_key,
                    numint.subintervals    = numint.subintervals,
                    interp.method          = interp.method,
                    interp.num_support     = interp.num_support,
                    cov.reparameterize    = cov.reparameterize,
                    cov.eps                = cov.eps,
                    cov.d_value            = cov.d_value,
                    ...)

  return(cov.args)
}

control.askey <- function(cov.reparameterize = getOption("wendland.cov.reparameterize"),
                          cov.eps = getOption("wendland.cov.eps"),
                          cov.d_value = getOption("wendland.cov.d_value"), ...){
  cov.args <- list(cov.reparameterize = cov.reparameterize, cov.eps = cov.eps,
                   cov.d_value = cov.d_value, ...)
  return(cov.args)
}
