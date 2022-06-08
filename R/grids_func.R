#' Grids Function
#'
#' @description This generates the grids of values for the explanatory variable.
#'
#' @param x The explanatory data.
#' @param X.type
#' @param nxgrid
#' @param grid.user
#' @param gridx.cont
#'
#' @return
#' @export
#'
#' @examples
grids_func <- function(x, X.type, nxgrid=101, grid.user = NULL, gridx.cont=F){
  # For multivariate continuous X, grid.user is a list of evaluation points for Xs - not including the first continuous X.
  # x=x; X.type=X.type; nxgrid=ng.plot; grid.user=grid.user; gridx.cont=T

  x         <- data.frame(x)
  ncov      <- NCOL(x)
  X.c.dex   <- which(X.type=="continuous")
  ncov.fac  <- ncov-length(X.c.dex)   # Number of discrete covariates
  ncov.cont <- ncov-ncov.fac          # Number of continuous covariates
  xgrid     <- NULL  # IS THIS NEEDED?
  #xvals     <- NULL

  if(gridx.cont==T){
    xgrid <- list()
    # Grid for univariate X and QGM
    if(length(nxgrid)==1 && length(grid.user)==0){
      for(i in 1:ncov.cont){
        xgrid[[i]] <- seq(min(x[X.c.dex[i]]),max(x[X.c.dex[i]]),len=nxgrid)
      }
    }

    # Grid for plots with multivariate X
    if(length(nxgrid)>1 && length(grid.user)==0){
      xgrid[[1]] <- seq(min(x[X.c.dex[1]]),max(x[X.c.dex[1]]),len=nxgrid[1])
      for(i in 2:ncov.cont){
        qx         <- seq(1/(nxgrid[i]+1),nxgrid[i]/(nxgrid[i]+1),len=nxgrid[i])
        xgrid[[i]] <- quantile(as.matrix(x[X.c.dex[i]]),probs=qx)
      }
    }

    if(length(nxgrid)==1 && length(grid.user)!= length(X.type) && length(grid.user)>0){
      xgrid[[1]] <- seq(min(x[X.c.dex[1]]),max(x[X.c.dex[1]]),len=nxgrid[1])
      for(i in 2:ncov.cont){
        xgrid[[i]] <- grid.user[[i-1]]
      }
    }

    if(length(nxgrid)==1 && length(grid.user) == length(X.type)){
      for(i in 1:ncov.cont){
        xgrid[[i]] <- grid.user[[i]]
      }
    }

    xgrid <- expand.grid(xgrid)
  }

  if(gridx.cont==F){ xgrid <- x[X.c.dex] }

  if(ncov.fac==0){ xgrid.list <- list(xgrid); xvals.list <- NULL }

  if(ncov.fac>0){
    # Define X grids for discrete variables
    X.fac.dex    <- which(X.type=="discrete")
    xvals        <- list()
    for(i in 1:length(X.fac.dex)){
      if(length(grid.user)==0){
        xvals[[i]] <- sort(unique(unlist(x[X.fac.dex[i]])))
      }
      if(length(grid.user)>0){
        if(length(grid.user) != length(X.type)){
          xvals[[i]] <- sort(unique(grid.user[[(ncov.cont-1)+i]]))
        }
        if(length(grid.user)>0 && length(grid.user) != length(X.type) && length(grid.user)>0){
          xvals[[i]] <- sort(unique(grid.user[[(ncov.cont-1)+i]]))
        }
        if(length(grid.user)>0 && length(grid.user) == length(X.type)){
          xvals[[i]] <- sort(unique(grid.user[[ncov.cont+i]]))
        }
      }
    }
    xvals.list <- as.matrix(expand.grid(xvals))

    if(ncov.cont>0){
      xgrid.list <- list()
      for(i in 1:nrow(xvals.list)){
        xgrid.list[[i]] <- xgrid
        for(j in 1:ncol(xvals.list)){
          xgrid.list[[i]] <- cbind(xgrid.list[[i]],xvals.list[i,j])
        }
      }
    }

    if(ncov.cont==0){
      xgrid.list <- list()
      for(i in 1:nrow(xvals.list)){
        xgrid.list[[i]] <- xgrid
        for(j in 1:ncol(xvals.list)){
          xgrid.list[[i]] <- cbind(xgrid.list[[i]],xvals.list[i,j])
        }
      }
    }

  }

  ans <- list(xgrid=xgrid.list, xvals=xvals.list)

  return(ans)
}
