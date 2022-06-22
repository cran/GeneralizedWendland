
######################################################################
# Function factory which facilitates access to cholesky decompositions
# with given parameters
#---------------------------------------------------------------------

choleskyFactory <- function(chol.args = list(), Rstruct = NULL){

  cholSpam <- function(call.args){
    if (is.null(Rstruct) || !is(Rstruct, "spam.chol.NgPeyton")){
      chol_method <- "chol.spam"
    } else {
      chol_method <- "update.spam.chol.NgPeyton"
      call.args[["object"]] <- Rstruct
    }

    return(do.call(chol_method, call.args))
  }

  choleskyFunction <- function(Sigma) {

    call.args <- c(list(x = Sigma), chol.args)

    if (spam::is.spam(Sigma)){
      cholS <- tryCatch(cholSpam(call.args),
        error = function(err) {
          assign("Rstruct", NULL)
          return(cholSpam(call.args))
        })

      if (is.null(Rstruct)) {
        assign("Rstruct", cholS)
      }

    } else {
      tol <- call.args[["eps"]]
      pivot <- call.args[["pivot"]]

      if (is.null(pivot)) {
        pivot <- FALSE
      } else {
        pivot <- ifelse(is.logical(pivot), pivot, TRUE)
      }

      cholS <- chol(Sigma, pivot = pivot, tol = tol)
    }

    return(cholS)
  }

  return(choleskyFunction)
}
