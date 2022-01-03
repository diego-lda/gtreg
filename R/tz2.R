#' TZ Function generator (v2)
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
TZ2 <- function(X,Y){
  #X is nrx by nx
  #Y is nry by ny
  #produce TZ nry by (nx*ny)

  if(!is.matrix(Y)){Y <- matrix(Y,nr=1)}
  if(!is.matrix(X)){X <- matrix(X,nr=nrow(Y),nc=length(X),byrow=T)}
  nx <- ncol(X)
  ny <- ncol(Y)
  nobs <- nrow(X)
  TZ <- matrix(0,nrow=nobs,ncol=nx*ny)
  i <- 0

  for(j in 1:ny){
    for(k in 1:nx){
      i <- i+1
      TZ[,i] <- X[,k]*Y[,j]
    }}
  return(TZ)
}

