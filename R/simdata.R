#' Simulate Data
#'
#' @description This function generates some fake data.
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
#' @return Some fake data.
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
