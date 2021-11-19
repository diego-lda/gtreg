
#' @title califunc
#'
#' @description This function is for simulations. It calibrates the coefficients.
#'
#' @param Y
#' @param X
#' @param qY
#' @param dgp
#' @param link
#'
#' @return
#' @export
#'
#' @examples
califunc <- function(Y, X, qY = NULL, dgp = "dr", link = "probit"){;

  if(dgp == "loc"){

    coef       <- lm(Y ~ X)$coef
    intercoef  <- c(coef[1],1);      # Calibrated coefs for intercept
    Xcoef	     <- c(coef[2],0);      # Calibrated coefs for Z
    eps        <- Y - (coef[1]+coef[2]*X)
    sigma      <- sd(eps)
  }

  if(dgp == "dr"){

    gridY       <- quantile(Y,qY)     # Grid of Y values for distribution regression estimator
    coefs       <- distrfun(Y = Y, X = X, gridY = gridY, link = link)$coefs0

    intercoef   <- coef( lm(coefs[1,] ~ gridY) );      # Calibrated coefs for intercept
    Xcoef	   	  <- coef( lm(coefs[2,] ~ gridY) );      # Calibrated coefs for X

    sde         <- NULL
  }


  # Store calibrated coefs
  bvec <- c( intercoef[1], Xcoef[1], intercoef[2], Xcoef[2] )

  ans <- list(bvec = bvec, sde = sde)

  return(ans)

}


#' @title simdata
#'
#' @description This function is for simulations. It generates fake data that simulates the one provided.
#'
#'
#' @param Ydata
#' @param Xdata
#' @param bvec
#' @param n
#' @param sde
#' @param dgp
#' @param ugrid
#' @param xgrid
#'
#' @return
#' @export
#'
#' @examples
simdata <- function(Ydata, Xdata, bvec, n = 100, sde = 1, dgp = "dr", ugrid = NULL, xgrid = NULL){

  Dat     <- NULL;
  Dat$eps <- rnorm(n);
  Dat$X   <- runif(n,min=min(Xdata),max=max(Xdata))

  if(dgp == "dr"){
    Dat$Y   <- as.numeric(-((bvec[1:2] %*% rbind(1,Dat$X)) / (bvec[3:4] %*%rbind(1,Dat$X))) + (1/(bvec[3:4] %*%rbind(1,Dat$X)))*Dat$eps)
    if(length(xgrid)>0){
      Dat$Qyx <- matrix(0,length(xgrid),length(ugrid))
      for(i in 1:length(ugrid)){
        Dat$Qyx[,i] <- -(bvec[1:2] %*% rbind(1,xgrid)) / (bvec[3:4] %*%rbind(1,xgrid)) + (1/(bvec[3:4] %*% rbind(1,xgrid)))*qnorm(ugrid[i])
      }
      Dat$Fyx <- pnorm((bvec[c(1,3)] %*% t(cbind(1,Dat$Y))) + (bvec[c(2,4)] %*%t(cbind(1,Dat$Y))) * Dat$X)
    }
  }

  ifelse(length(xgrid)>0, ans <- list(Y = Dat$Y, X = Dat$X, Qyx = Dat$Qyx, Fyx = Dat$Fyx ),
         ans <- list(Y = Dat$Y, X = Dat$X ))

  return(ans)
}


#' @title TZ1
#'
#' @description This function is internal. It combines X and Y to generate TZ.
#'
#' @param X
#' @param Y
#'
#' @return
#' @export
#'
#' @examples
TZ1 <- function(X,Y){
  #X is nobs by nx
  #Y is nobs by ny
  #produce TZ nobs by (nx*ny)

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
