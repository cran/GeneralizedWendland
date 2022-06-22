predictionFactory <- function(y, locs0, locs1, covariance,
  X0 = list2DF(nrow = nrow(locs0)), X1 = list2DF(nrow = nrow(locs1)),
  ..., cov.args = list(), chol.args = list(), use_spam = TRUE){

  covFun <- covarianceFactory(covariance = covariance, cov.args = cov.args)
  cholFun <- choleskyFactory(chol.args)

  lpFun <- function(X, param){
    if (dim(X)[2] == 0){
      return(rep(0, dim(X)[1]))
    } else {
      return(X %*% param[seq_len(dim(X)[2])])
    }
  }

  distFun <- function(l1, l2, delta = NULL) {
    if (use_spam){
      return(spam::nearest.dist(l1, l2, delta = delta))
    } else {
      return(fields::rdist(l1, l2))
    }
  }

  solveFun <- function(cholS, rhs){

    if (use_spam){
      return(spam::backsolve.spam(cholS, spam::forwardsolve.spam(cholS, rhs,
              transpose = TRUE, upper.tri = TRUE)))
    } else {
      perm <- iperm <- diag(rep(1,dim(cholS)[2]))
      if (!is.null((pivots <- attr(cholS, which = "pivot")))){
        perm <- perm[pivots,]
        iperm <- iperm[order(pivots),]
      }
      fsolve <- base::forwardsolve(cholS, perm %*% rhs,
        upper.tri = TRUE, transpose = TRUE)
      return(iperm %*% backsolve(cholS, fsolve, dim(cholS)[2]))
    }
  }

  simulate <- function(n, param){
    theta <- param[seq(dim(X0)[2]+1, length(param))]
    if (use_spam){
      range <- c(cov.args[["fixed_range_value"]], theta[1])[1]
    } else {
      range <- NULL
    }

    Sigma0 <- covFun(distFun(locs0,locs0,range), theta = theta)
    Sigma1 <- covFun(distFun(locs1,locs1,range), theta = theta)
    Sigma2 <- covFun(distFun(locs0,locs1,range), theta = theta)

    cholS <- cholFun(Sigma0)
    mu_cond <- lpFun(X1,param) + t(Sigma2) %*% solveFun(cholS, y - lpFun(X0,param))
    Si_cond <- Sigma1 - t(Sigma2) %*% solveFun(cholS, Sigma2)
    return(mvtnorm::rmvnorm(n, mu_cond, Si_cond))
  }

  return(simulate)
}
