#' Beta2 Check
#'
#' @description This function computes the slope and its minimum value.
#'
#' @param bmat
#' @param Xs
#' @param nXs
#' @param nYS
#'
#' @return The slope and its minimum valye
#' @export
#'
#' @examples
beta2.check <- function(bmat,Xs,nXs,nYS){

  beta  <- Xs%*%matrix(bmat,nr=nXs,nc=nYS)

  b2min  <- min(beta[,2])
  beta2  <- beta[,2]

  return(list(b2min=b2min,beta2=beta2))

}
