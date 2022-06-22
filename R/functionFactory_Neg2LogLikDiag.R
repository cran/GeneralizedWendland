
neg2loglikDiagFactory <- function(y, X = data.frame(), distmat, covariance, ...){

  covname <- deparse(substitute(covariance))

  function(theta_list, cov.args_list = list(), chol.args_list = list()) {

    theta_grid <- do.call(expand.grid, c(theta_list))
    args_grid <- do.call(expand.grid, c(cov.args_list, chol.args_list))


    results <- list()

    for (i in seq_len(nrow(args_grid))) {
      args_row <- data.frame(args_grid[i,])
      colnames(args_row) <- colnames(args_grid)
      cov.args_row <- args_row[names(cov.args_list)]
      chol.args_row <- args_row[names(chol.args_list)]

      n2ll <- neg2loglikFactory(y = y,
                                X = X,
                                distmat = distmat,
                                covariance = covariance,
                                cov.args = as.list(cov.args_row),
                                chol.args = as.list(chol.args_row))

      result <- numeric(nrow(theta_grid))

      for (j in seq_len(nrow(theta_grid))) {
        theta_row <- theta_grid[j,]
        result[j] <- n2ll(unlist(theta_row))
      }

      results[[i]] <- cbind(theta_grid, as.list(args_row), result = result)
    }

    return(do.call("rbind", results))
  }

}
