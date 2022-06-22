######################################################################
# Diagnostic functions for generalized wendland covariance function
# approximations obtained through interpolation.
# Author: Thomas Caspar Fischer
#---------------------------------------------------------------------

covDiagFactory <- function(target_covariance,
                           diagnostic_funs = c("accumulated_error",
                                               "point_diagnostics"),
                           reference_covariance = cov.askey,
                           reference_cov.args = list()
){

  reference_covFun <- covarianceFactory(reference_covariance,
                                        reference_cov.args)

  inner_logic <- function(FUN, target_covFun, target_cov.theta, target_cov.args,
                          reference_cov.theta, ...){
    result <- do.call(FUN,
                      args = c(list(...),
                               list(target_covFun = target_covFun,
                                    target_cov.theta = target_cov.theta,
                                    reference_covFun = reference_covFun,
                                    reference_cov.theta = reference_cov.theta)))

    # Output management

    target.pars <- c(as.list(target_cov.theta), target_cov.args)
    names(target.pars) <- paste("target", ifelse(sapply(names(target.pars), function(x) is.null(x) || (x=="")),
                                                 paste0("theta", 1:length(target.pars)),
                                                 names(target.pars)),
                                sep = ".")
    reference.pars <- c(as.list(reference_cov.theta), reference_cov.args)
    names(reference.pars) <- paste("reference", ifelse(sapply(names(reference.pars), function(x) is.null(x) || (x=="")),
                                                 paste0("theta", 1:length(reference.pars)),
                                                 names(reference.pars)),
                                sep = ".")
    other <- list(...)
    run_information <- c(target.pars, reference.pars, other)

    output <- data.frame(c(run_information, result))
    return(output)
  }


  outer_logic <- function(target_theta_list, target_args_list = list(),
                reference_cov.theta = NULL, ...){

    npar_theta <- length(target_theta_list)
    grid <- do.call(expand.grid, c(target_theta_list,
                                   target_args_list))

    results <- list()

    for (i in seq_len(nrow(grid))){
      target_cov.args <- as.list(grid[i, -seq_len(npar_theta)])
      target_cov.theta <- grid[i, seq_len(npar_theta)]

      target_covFun <- covarianceFactory(covariance = target_covariance,
                                         cov.args = target_cov.args)

      valid_args <- tryCatch(is.numeric(target_covFun(0, unlist(target_cov.theta))),
                             warning = function(w) FALSE,
                             error = function(e) FALSE)

      if (!valid_args)
        next

      if (is.null(reference_cov.theta)) {
        if (identical(reference_covariance, target_covariance)) {
          reference_cov.theta <- target_cov.theta
        } else if (identical(reference_covariance, cov.askey)) {
          reference_cov.theta <- target_cov.theta[names(target_cov.theta)[c(1,2,4,5)]]
        } else {
          stop("No reference_cov.theta provided for reference covariance of different type than target")
        }
      }

      for (j in seq_along(diagnostic_funs)) {
        fun <- diagnostic_funs[j]

        if (is.null(results[[fun]])){
          results[[fun]] <- data.frame()
        }

        results[[fun]] <- rbind(results[[fun]], inner_logic(FUN=fun, target_covFun = target_covFun,
                                target_cov.theta = target_cov.theta,
                                target_cov.args = target_cov.args,
                                reference_cov.theta = reference_cov.theta, ...))
      }
    }
    return(results)
  }
  return(outer_logic)
}


accumulated_error <- function(
  target_covFun,
  target_cov.theta,
  reference_covFun,
  reference_cov.theta,
  ...,
  absolute = TRUE,
  lower = 0,
  upper = 1,
  subdivisions = 500L,
  abs.tol = .Machine$double.eps^0.5,
  rel.tol = .Machine$double.eps^0.25
){

  f <- function(x) {

    target_cov <- target_covFun(x, unlist(target_cov.theta))
    reference_cov <- reference_covFun(x, unlist(reference_cov.theta))
    error <- target_cov - reference_cov

    if (absolute)
      return(abs(error))
    else
      return(error)
  }

  result <- tryCatch(integrate(f, lower = lower, upper = upper,
                               subdivisions = subdivisions,
                               rel.tol = rel.tol,
                               abs.tol = abs.tol)$value,
                     warning = function(w) return(NA),
                     error = function(e) return(NA))
  return(list(accumulated_error = result))
}



point_diagnostics <- function(
  target_covFun,
  target_cov.theta,
  reference_covFun,
  reference_cov.theta,
  ...,
  grid_resolution = 100
){

  h <- (seq_len(grid_resolution+1)-1)/grid_resolution

  target_cov    <- target_covFun(h = h, theta = unlist(target_cov.theta))
  reference_cov <- reference_covFun(h = h, theta = unlist(reference_cov.theta))
  error <- target_cov - reference_cov
  absolute_error <- abs(error)
  relative_error <- error/reference_cov
  max_absolute_error <- max(abs(error))

  return(list(h = h, target_cov = target_cov, reference_cov = reference_cov,
              error = error, absolute_error = absolute_error,
              relative_error = relative_error,
              max_absolute_error = max_absolute_error))
}