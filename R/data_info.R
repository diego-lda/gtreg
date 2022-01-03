#' Data Info Function
#'
#' @description This function creates the necessary info for the data.
#'
#' @param x
#' @param X.type
#' @param xdf
#' @param xknots
#' @param coord.bare
#' @param coord.spline
#' @param coord.tensor
#' @param nxgrid
#' @param gridx.cont
#' @param delta.ok
#'
#' @return The information for this data.
#' @export
#'
#' @examples
data.info <- function(x, X.type, xdf, xknots=NULL, coord.bare=NULL, coord.spline=NULL, coord.tensor=NULL, nxgrid=101, gridx.cont=F,delta.ok=F){

  xknots.missing <- F
  if(length(xknots)==0){ xknots.missing <- T }
  # x should have continuous variables first.
  ncov      <- NCOL(x)
  if(length(coord.bare) != 0 || length(coord.spline) != 0 ){
    # Check that each variable's type is defined
    stopifnot(length(X.type)==ncov)
    # Check that each variable's form is defined
    stopifnot(mean(1:ncov %in% c(coord.bare,coord.spline))==1)
    # Check that each variable's form is defined only once
    stopifnot(length(intersect(coord.bare,coord.spline))==0)
  }

  #x         <- data.frame(x)
  x         <- as.matrix(x)
  X.c.dex   <- which(X.type=="continuous")
  ncov.fac  <- length(X.c.dex)   # Number of discrete covariates
  ncov.cont <- ncov-ncov.fac     # Number of continuous covariates
  info      <- list()
  ntensors  <- 0
  if(!is.null(coord.tensor)){ ntensors <- length(coord.tensor) }

  # Define info list for W(X) construction
  for(i in coord.bare){
    ifelse(X.type[i]=="continuous",
           info[[i]] <- list(type = "bare", coords = i, center = mean(x[,i])),
           info[[i]] <- list(type = "bare", coords = i) )
  }

  for(i in coord.spline){
    # ifelse(delta.ok, delta <- (max(x[,i]) - min(x[,i]))/4, delta <- 0)
    # if(delta.ok && min(x[,i])>=0 && abs(delta)>=min(x[,i])){ delta <-  min(x[,i]) }
    # if(xknots.missing){ xknots <- seq(min(x[,i])-abs(delta),max(x[,i])+abs(delta),length.out=xdf[i]) }
    # if(xknots.missing){ xknots <- seq(min(x[,i])-abs(delta),max(x[,i])+abs(delta),length.out=xdf[i]) }
    #xabsmax <- max(abs(x[,i]))
    #if(xknots.missing){ xknots <- seq(min(x[,i])-xabsmax,max(x[,i]),length.out=xdf[i]) }
    if(!delta.ok)if(xknots.missing){ xknots <- seq(min(x[,i]),max(x[,i]),length.out=xdf[i]) }
    if(delta.ok){
      delta <- (max(x[,i]) - min(x[,i]))/4
      #if(min(x[,i])>=0 && abs(delta)>=min(x[,i])){ delta <-  min(x[,i]) }
      if(min(x[,i])>=0 && abs(delta)>=min(x[,i])){ delta <-  0 }
      if(xknots.missing){ xknots <- seq(min(x[,i])-abs(delta),max(x[,i])+abs(delta),length.out=xdf[i]) }
    }
    info[[i]] <- list(type = "spline", knots = xknots, coords = i, center = mean(x[,i]))
  }

  if(ntensors>0){

    for(j in 1:ntensors){

      linfo         <- length(info)+1
      info[[linfo]] <- list()

      for(i in 1:length(coord.tensor[[j]])){

        tdex <- coord.tensor[[j]][i]

        if(tdex %in% coord.bare){
          ifelse(X.type[tdex]=="continuous",
                 info[[linfo]][[i]] <- list(type = "bare", coords = tdex, center = mean(x[,tdex])),
                 info[[linfo]][[i]] <- list(type = "bare", coords = tdex) )
        }
        else{
          # ifelse(delta.ok, delta <- (max(x[,tdex]) - min(x[,tdex]))/4, delta <- 0)
          # if(delta.ok && min(x[,tdex])>=0 && abs(delta)>=min(x[,tdex])){ delta <-  min(x[,tdex]) }
          # if(xknots.missing){ xknots <- seq(min(x[,tdex])-abs(delta),max(x[,tdex])+abs(delta),length.out=xdf[tdex]) }
          #if(xknots.missing){ xknots <- seq(min(x[,tdex]),max(x[,tdex]),length.out=xdf[tdex]) }
          if(!delta.ok){ if(xknots.missing){ xknots <- seq(min(x[,tdex]),max(x[,tdex]),length.out=xdf[tdex]) } }
          if(delta.ok){
            delta <- (max(x[,tdex]) - min(x[,tdex]))/4
            #if(min(x[,tdex])>=0 && abs(delta)>=min(x[,tdex])){ delta <-  min(x[,tdex]) }
            if(min(x[,tdex])>=0 && abs(delta)>=min(x[,tdex])){ delta <-  0 }
            if(xknots.missing){ xknots <- seq(min(x[,tdex])-abs(delta),max(x[,tdex])+abs(delta),length.out=xdf[tdex]) }
          }
          #xabsmax <- max(abs(x[,tdex]))
          #if(xknots.missing){ xknots <- seq(min(x[,tdex])-xabsmax,max(x[,tdex]),length.out=xdf[tdex]) }
          info[[linfo]][[i]] <- list(type = "spline", knots = xknots, coords = tdex, center = mean(x[,tdex]))
        }

      }

      info[[linfo]]$type <- "tensor"

    }

  }

  return(info)
}
