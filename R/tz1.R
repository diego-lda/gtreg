#' TZ Function generator (v1)
#'
#' @description This function generates the TZ matrix.
#'
#' @param X
#' @param Y
#'
#' @return The TZ function.
#' @export
#'
#' @examples
TZ1 <- function(X,Y){

  if(!is.matrix(Y)){Y <- matrix(Y,nr=1)}
  nx <- ncol(X)
  ny <- ncol(Y)
  nobs <- nrow(X)
  TZ <- matrix(0,nrow=nobs,ncol=nx*ny)
  i  <- 0

  for(j in 1:ny){
    for(k in 1:nx){
      i <- i+1
      TZ[,i] <- X[,k]*Y[,j]
    }}
  return(TZ)
}
