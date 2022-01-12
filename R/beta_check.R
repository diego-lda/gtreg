#' Beta Check
#'
#' @description This function computes the values for Beta2(X) and returns them highlighting the minimum value.
#'
#' @param bmat A matrix with the betas.
#' @param Xs
#' @param nXs
#' @param nYS
#'
#' @return A list containing the minimum of Beta2(X) and the vector of Beta2(X).
#' @export
#'
#' @examples
beta_check <- function(bmat,Xs,nXs,nYS){

  beta  <- Xs%*%matrix(bmat,nr=nXs,nc=nYS)

  b2min  <- min(beta[,2])
  beta2  <- beta[,2]

  return(list(b2min=b2min,beta2=beta2))

}
