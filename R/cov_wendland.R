
cov.wendland <- function(h, theta, ..., cov.args = list()){
  args <- c("numint.abstol", "numint.reltol", "numint.qag_key",
            "numint.subintervals", "interp.method", "interp.num_support",
            "cov.reparameterize", "cov.eps", "cov.d_value", "cov.sparse")
  interp.all <- c("none", "linear", "polynomial", "cspline")

  for (arg in args) {
    if (!(arg %in% names(cov.args))) {
      cov.args[[arg]] <- getOption(paste0("wendland.", arg))
    }
  }

  interp.type <- which(cov.args[["interp.method"]] == interp.all) - 1

  checks <- c(length(interp.type) == 0,
              cov.args[["interp.method"]] == "none" & cov.args[["interp.num_support"]] > 0,
              cov.args[["interp.method"]] != "none" & cov.args[["interp.num_support"]] < 3)
  errors <- c("Unrecognized interpolation method.",
              "No interpolator chosen, but num_support > 0",
              "Interpolation requires num_support >= 3")

  for (i in 1:length(checks)) {
    if (checks[i]) stop(errors[i])
  }

  range  <- theta[1]
  sill   <- theta[2]
  kappa  <- theta[3]
  mu     <- theta[4]+cov.args[["cov.reparameterize"]]*(1+cov.args[["cov.d_value"]]+2*kappa)/2
  nugget <- theta[5]

  wend <- new("Rcpp_Wendland")
  wend$setParameters(range, sill, kappa, mu, nugget)
  wend$setEpsTol(cov.args[["cov.eps"]])
  wend$setIntegrator(cov.args[["numint.abstol"]], cov.args[["numint.reltol"]],
                     cov.args[["numint.subintervals"]], cov.args[["numint.qag_key"]])
  wend$setInterpolator(cov.args[["interp.num_support"]], interp.type)


  if (is.vector(h)){
    computeFunction <- function(d){
      return(wend$computeVector(d))
    }

  } else if (is.matrix(h)){
    computeFunction <- function(d){
      result <- wend$computeMatrix(d)
      if (cov.args[["cov.sparse"]]) result <- as.spam(result, eps = cov.args[["cov.eps"]])

      return(result)
    }

  } else if (spam::is.spam(h)){
    computeFunction <- function(d){
      triplet <- spam::triplet(d)
      triplet <- wend$computeSpam(triplet[["indices"]], triplet[["values"]])

      return(spam::spam(list(i=triplet[["indices"]][,1],
                             j=triplet[["indices"]][,2],
                             values = triplet[["values"]])))
    }

  } else if (is(h, "dgCMatrix")){
    computeFunction <- function(d){
      return(wend$computeSparse(d))
    }

  } else {
    stop("Unknown datatype")
  }

  return(computeFunction(h))
}
