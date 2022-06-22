
mleFactory <- function(covariance, cov.args = list(), chol.args = list(),
  optim.args = list(), hessian = FALSE,
  optimParallel.args = list()){

  covarianceFunction <- covarianceFactory(covariance = covariance,
    cov.args = cov.args)
  optimFunction <- optimFactory(optim.args = optim.args, hessian = hessian,
    optimParallel.args = optimParallel.args)


  mleFunction <- function(y, X = data.frame(), distmat, init_parameters,
    theta_llim, theta_ulim){

    cholFunction <- choleskyFactory(chol.args=chol.args, Rstruct=Rstruct)
    theta <- init_parameters
    beta <- numeric(ncol(X))
    Rstruct <- NULL


    neg2loglikFunction <- neg2loglikFactory(y = y, X = X, distmat = distmat,
      covarianceFunction = covarianceFunction, choleskyFunction = cholFunction)

    # Workaround for lazy evaluation behaviour.
    # neg2loglikFunction <- function(parameters) {
    #   return(neg2loglikFunction0(parameters))
    # }

    init  <- c(beta, theta)
    lower <- c(rep(-Inf, ncol(X)), theta_llim)
    upper <- c(rep( Inf, ncol(X)), theta_ulim)

    ######################################################################
    # Run maximum likelihood estimation
    #----------------------------------------------------------------------

    result <- optimFunction(par = init, fn = neg2loglikFunction, gr = NULL,
      lower = lower, upper = upper)

    return(result)
  }

  return(mleFunction)
}
