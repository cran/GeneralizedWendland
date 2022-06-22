
covarianceFactory <- function(covariance, cov.args = list()){

  #######################################################################
  # Run checks
  # [Note: reasonable to allow specifying function by character argument?]
  #----------------------------------------------------------------------

  covariance <- tryCatch(match.fun(covariance),
                         error = function(cond){
                           stop("Covariance function does not exist",
                                "or is not a valid function.")
  })

  fixed_range  <- cov.args[["fixed_range_value"]]
  fixed_nugget <- cov.args[["fixed_nugget_value"]]


  #######################################################################
  # Function blueprint
  #----------------------------------------------------------------------
  covarianceFunction <- function(h, theta, ...) {

    theta <- c(fixed_range, theta, fixed_nugget)
    result <- covariance(h = h, theta = theta, ..., cov.args = cov.args)

    return(result)
  }

  return(covarianceFunction)
}
