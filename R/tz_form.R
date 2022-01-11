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
tz_form <- function(X,Y){

  if(!is.matrix(Y)){Y <- matrix(Y,nr=1)}
  if(!is.matrix(X)){X <- matrix(X,nr=nrow(Y),nc=length(X),byrow=T)} #This line is removed in what used to be TZ1.
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

