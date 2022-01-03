#' tZB Function (v2)
#'
#' @description This function makes a giant constraint matrix, with Ygrid the same for every observation. Ygrid has already been processed into sYgrid, i.e. s(Ygrid), which is a matrix nYgrid by nSY.
#'
#' @param Xs
#' @param SYfunc
#' @param sYfunc
#' @param Ygrid
#' @param nsim
#'
#' @return
#' @export
#'
#' @examples
tZB2 <- function(Xs,SYfunc,sYfunc,Ygrid,nsim=29){

  ylo <- min(Ygrid)
  yhi <- max(Ygrid)
  Ygrid.look <<- Ygrid
  trysY <- EvalBasis(object=sYfunc,x=range(Ygrid))
  nsY <- ncol(trysY)+1  #sYfunc produces an object with nsY columns
  #add 1 because sY basis starts with [0,1,...]
  nobs <- nrow(Xs)
  nWx <- ncol(Xs)
  #nYgrid <- nrow(sYgrid)
  nYgrid <- nsim+1
  print(nYgrid)
  #nsY <- ncol(sYgrid)

  #OK, allocate a monster matrix A
  nrA <- nobs*(nYgrid)
  ncA <- nWx*nsY   #aka nb, the number of coefs in e=b.T(Y,X)

  A <- matrix(0,nr=nrA,nc=ncA)
  #compute A by observation
  #I can write fortran in any language
  #ysim.mat <- matrix(0,nr=nobs,nc=nsim+1)
  ntotysim <- (nobs*nsim)
  ydex <- sample(1:ntotysim,ntotysim)
  BYgrid <- seq(ylo,yhi,len=ntotysim)
  BYgrid.s <- BYgrid[ydex]
  midval <- 0
  Ytake.look <<- NULL
  BYgrid.look <<- BYgrid
  BYgrid.s.look <<- BYgrid.s

  for(i in 1:nobs){
    sdex <- (i-1)*nYgrid+1
    edex <- i*nYgrid
    sdex2 <- (i-1)*nsim+1
    edex2 <- i*nsim

    ytake <- sort(c(BYgrid.s[sdex2:(edex2)],midval))
    #Ytake.look <<- rbind(Ytake.look,ytake) #tmp just to check
    if(i <4) print(ytake)
    #ysim <- sort(runif(nsim,ylo,yhi))
    sYnow <- EvalBasis(object=sYfunc,x=ytake)
    #do the replacement
    #SYgrid <- cbind(1,Ygrid,SYy[,-1])
    sYgrid <- cbind(0,1,sYnow[,-1])
    Xsnow <- matrix(rep(Xs[i,],times=nYgrid),nr=nYgrid,byrow=T)
    #Xsnow.look <<- Xsnow

    #look1 <<- TZ1(Xsnow,sYgrid)
    A[sdex:edex,] <- (TZ1(Xsnow,sYgrid))
  }
  cat("\nA is set with dim(A):=",dim(A),"\n")
  return(A)
}
