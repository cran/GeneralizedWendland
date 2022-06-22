
optimFactory <- function(optim.args = list(), hessian = FALSE,
                         optimParallel.args = list()
){

  num_cores_requested <- c(optimParallel.args[["num_cores"]], 0)[1]
  forward <- c(optimParallel.args[["forward"]], FALSE)[1]
  loginfo <- c(optimParallel.args[["loginfo"]], FALSE)[1]
  cl_type <- ifelse(.Platform$OS.type == "windows", "PSOCK", "FORK")

  if (num_cores_requested > 1){
    num_cores_installed <- parallel::detectCores()
    num_cores_used <- min(num_cores_installed - 1, num_cores_requested)

  } else {
    num_cores_used <- ifelse(forward | loginfo, 1, 0)
  }

  if (num_cores_used >= 1) {
    parallel.args <- list(cl = NULL, forward = forward, loginfo = loginfo)

    optimFunction <- function(par, fn, gr = NULL, ..., lower, upper){

      fn_envir <- environment(fn)
      cl <- parallel::makeCluster(num_cores_used, type = cl_type)
      parallel::setDefaultCluster(cl = cl)
      parallel::clusterExport(cl, varlist = ls(fn_envir), envir = fn_envir)
      on.exit(add = TRUE, after = TRUE, expr = {
        parallel::stopCluster(cl)
        parallel::setDefaultCluster(NULL)
      })
      parallel.args[["cl"]] <- cl

      return(optimParallel(par = par, fn = fn, gr = gr, ...,
               lower = lower, upper = upper, control = optim.args,
               hessian = hessian, parallel = parallel.args))
    }

  } else {

    optimFunction <- function(par, fn, gr = NULL, ..., lower, upper) {
      return(optim(par = par, fn = fn, gr = gr, ..., method = "L-BFGS-B",
               lower = lower, upper = upper, control = optim.args,
               hessian = hessian))
    }
  }

  return(optimFunction)
}
