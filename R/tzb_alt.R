#' tZB Function (v1)
#'
#' @description This function makes a giant constraint matrix, with Ygrid the same for every observation. Ygrid has already been processed into sYgrid, i.e. s(Ygrid), which is a matrix nYgrid by nSY.
#'
#' @param Xs
#' @param sYgrid
#'
#' @return A giant constraint matrix.
#' @export
#'
#' @examples
tzb_alt <- function(Xs,sYgrid){

  nobs <- nrow(Xs)
  nWx <- ncol(Xs)
  nYgrid <- nrow(sYgrid)
  nsY <- ncol(sYgrid)

  #OK, allocate a monster matrix A
  nrA <- nobs*nYgrid
  ncA <- nWx*nsY   #aka nb, the number of coefs in e=b.T(Y,X)

  A <- matrix(0,nr=nrA,nc=ncA)
  #compute A by observation
  #I can write fortran in any language
  for(i in 1:nobs){
    sdex <- (i-1)*nYgrid+1
    edex <- i*nYgrid
    #Xsnow <- NULL
    #for(jj in 1:nYgrid){
    #    Xsnow <- rbind(Xsnow,Xs[i,])
    #   }
    Xsnow <- matrix(rep(Xs[i,],times=nYgrid),nr=nYgrid,byrow=T)
    #Xsnow.look <<- Xsnow

    #look1 <<- TZ1(Xsnow,sYgrid)
    A[sdex:edex,] <- (TZ1(Xsnow,sYgrid))
  }
  cat("\nA is set with dim(A):=",dim(A),"\n")
  return(A)
}
