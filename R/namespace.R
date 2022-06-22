

# On loading, add global default arguments
.onLoad <- function(libname, pkgname) {
  defaults <- list(
    wendland.numint.abstol = 1e-3,
    wendland.numint.reltol = 1e-3,
    wendland.numint.qag_key = 0,
    wendland.numint.subintervals = 0,
    wendland.interp.method = "none",
    wendland.interp.num_support = 0,
    wendland.cov.reparameterize = TRUE,
    wendland.cov.eps = .Machine$double.eps^0.5,
    wendland.cov.d_value = 2,
    wendland.cov.sparse = FALSE
  )

  options(defaults)
  loadModule("Wendland")
}
