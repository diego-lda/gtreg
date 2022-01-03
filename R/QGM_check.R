#' Quasi-Gaussian Monotonicity Check
#'
#' @description This function checks whether there is quasi-Gaussian monotonicity in the gtr result.
#'
#' @param y
#' @param x
#' @param res
#' @param xgrid.qgm
#' @param ugrid
#' @param nxgrid
#' @param nygrid
#' @param ng.qgm
#' @param nXs
#' @param nYS
#' @param ydf
#' @param yorder
#' @param iyknots
#' @param info
#' @param Ysing
#' @param delta
#' @param easy
#' @param zeros
#' @param e0mode
#'
#' @return A design matrix with the result of the check.
#' @export
#'
#' @examples
QGM.check <- function(y,x,res,xgrid.qgm=seq(min(x),max(x),len=101),ugrid,nxgrid=101,nygrid,ng.qgm=201,nXs,nYS,ydf,yorder,iyknots=NULL,
                      info=info,Ysing=FALSE,delta=1,easy=F,zeros=NULL,e0mode=T){

  nobs        <- length(y)
  ygrid.qgm   <- seq(min(y),max(y),len=ng.qgm)
  M <- diag(1,nXs*nYS,nXs*nYS)
  if(length(zeros)>0){
    M            <- matrix(0,nXs*nYS,nXs*nYS-length(zeros))
    M[-zeros,]   <- diag(1,nXs*nYS-length(zeros),nXs*nYS-length(zeros))
  }
  Bmat        <- matrix(M%*%res$bmat,nr=nXs,nc=nYS)

  if(nYS==2){
    datmat.mono <- dataprep(y=as.numeric(y),x=x,xgrid=xgrid.qgm,ygrid=ygrid.qgm,bmat=Bmat,
                            info=info,iyknots=iyknots,ydf=ydf,addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,
                            Ysing=Ysing,e0mode=e0mode,returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,nxgrid=nxgrid,nygrid=nygrid)
    #    lmin        <- beta2.check(bmat=Bmat,Xs=datmat.mono$Xsgrid,nXs=nXs,nYS=nYS)$b2min
    if(!is.list(datmat.mono$Xsgrid)){
      lmin <- beta2.check(bmat=Bmat,Xs=datmat.mono$Xsgrid,nXs=nXs,nYS=nYS)$b2min
    }

    if(is.list(datmat.mono$Xsgrid)){
      lmin <- NULL
      for(i in 1:length(datmat.mono$Xsgrid)){
        for(j in 1:dim(datmat.mono$Xsgrid[[i]])[3]){
          # b2min is the smallest value of beta2(X). It should be strictly positive.
          lmin <- c(b2min,beta2.check(bmat=Bmat,Xs=datmat.mono$Xsgrid[[i]][,,j],nXs=nXs,nYS=nYS)$b2min)
        }
      }
    }

    datmat.mono$dedygrid <- min(lmin)
  }

  if(nYS!=2 && nXs==2){
    datmat.mono <- dataprep(y=as.numeric(y),x=x,xgrid=xgrid.qgm,ygrid=ygrid.qgm,bmat=Bmat,
                            info=info,iyknots=iyknots,ydf=ydf,addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,
                            Ysing=Ysing,e0mode=e0mode,returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,nxgrid=nxgrid,nygrid=nygrid)
    dBetaY <- datmat.mono$sYgrid%*%t(Bmat)
    #lmin   <- cbind(1,matrix(c(min(x)-.1*(max(x)-min(x)),max(x)+.1*(max(x)-min(x))),nr=2,nc=1))%*%t(dBetaY)
    lmin   <- cbind(1,matrix(c(min(x),max(x)),nr=2,nc=1))%*%t(dBetaY)
    datmat.mono$dedygrid <- min(lmin)
  }

  if(nYS!=2 && nXs>2){
    if(easy){
      datmat.mono <- dataprep(y=as.numeric(y),x=x,ygrid=ygrid.qgm,xgrid=xgrid.qgm,bmat=Bmat,
                              info=info,iyknots=iyknots,ydf=ydf,addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,
                              Ysing=Ysing,e0mode=e0mode,returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,nxgrid=nxgrid,nygrid=nygrid)
      mdedy <- NULL
      if(length(xgrid.qgm[[2]])>0){
        for(i in 1:nrow(xgrid.qgm[[2]])){
          mdedy       <- c(mdedy, min(datmat.mono$dedygrid[[i]]))
        }
      }
      #if(length(xgrid.qgm[[2]])==0){
      else{
        mdedy       <- min(unlist(datmat.mono$dedygrid))
      }
      datmat.mono$dedygrid <- min(mdedy)
    }

    if(!easy){

      ngrids <- length(xgrid.qgm$xgrid[[1]][,1])/nxgrid  # This is the number of replications of the first X grid - one for each grid value of other continuous X
      lmin   <- NULL
      i      <- 0

      #for(i in 1:length(xgrid.qgm$xgrid)){   # This is a loop over grid values of discrete x components.
      while(i < length(xgrid.qgm$xgrid)  && (lmin > 0 || length(lmin) == 0)){   # This is a loop over grid values of discrete x components.

        i <- i+1
        j <- 0

        #   for(j in 1:ngrids){             # This is a loop over grid values of continuous x components
        while(j < ngrids && (min(lmin) > 0 || length(lmin) == 0)){

          j           <- j+1

          nlev        <- 0
          kdex        <- 0
          maa.down    <- 0
          maa.up      <- 0

          ygrid.qgm   <- seq(min(y),max(y),len=ng.qgm)
          xgridnow    <- as.matrix(xgrid.qgm$xgrid[[i]][(1+nxgrid*(j-1)):(nxgrid*j),])

          while(nlev!=1 && (maa.down==0 || maa.up==0) && kdex <= 9){

            kdex        <- kdex+1
            datmat.mono <- dataprep(y=as.numeric(y),x=x,ygrid=ygrid.qgm,xgrid=xgridnow,bmat=Bmat,
                                    info=info,iyknots=iyknots,ydf=ydf,addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,
                                    Ysing=Ysing,e0mode=e0mode,returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,nxgrid=nxgrid,nygrid=nygrid)

            line.list   <- contourLines(xgridnow[,1],ygrid.qgm,pnorm(datmat.mono$egrid[[1]][,,1]),levels=ugrid)   # MODIFIED
            nlevel      <- length(line.list)

            umin.lev    <- NULL
            umax.lev    <- NULL

            for(l in 1:nlevel){
              if(line.list[[l]]$level==ugrid[1]){ umin.lev <- c(umin.lev, l) }
              if(line.list[[l]]$level==ugrid[length(ugrid)]){ umax.lev <- c(umax.lev, l) }
            }
            nlev.min    <- length(umin.lev)
            nlev.max    <- length(umax.lev)
            ranyg       <- range(ygrid.qgm)

            if(nlev.min!=1 || maa.down==0){ ygrid.qgm <- seq(min(ygrid.qgm)-5e-2*(ranyg[2]-ranyg[1]),max(ygrid.qgm),len=ng.qgm) }
            if(nlev.min!=1 || maa.up==0){   ygrid.qgm <- seq(min(ygrid.qgm),max(ygrid.qgm)+5e-2*(ranyg[2]-ranyg[1]),len=ng.qgm) }
            #print(c(kdex,ranyg))
            nlev        <- mean(c(nlev.min,nlev.max))

            if(maa.up==0){
              aa.up <- NULL
              for(ll in 1:nxgrid){
                aa.up <- c(aa.up,length((which(pnorm(datmat.mono$egrid[[1]][,,1])[ll,]>.99))))
              }
              maa.up   <- min(aa.up)
            }
            if(maa.down==0){
              aa.down <- NULL
              for(ll in 1:nxgrid){
                aa.down <- c(aa.down,length((which(pnorm(datmat.mono$egrid[[1]][,,1])[ll,]<.01))))
              }
              maa.down <- min(aa.down)
            }
          }

          for(l in 1:nlevel){

            if(ngrids == 1){
              xgrid.cqf <- line.list[[l]]$x
            }
            if(ngrids > 1){
              xgrid.cqf     <- matrix(0,length(line.list[[l]]$x),NCOL(x))
              xgrid.cqf[,1] <- line.list[[l]]$x
              for(m in 2:NCOL(x)){
                xgrid.cqf[,m] <- xgridnow[1,m]
              }
            }

            #for(lnow in 1:nlevel){print(range(as.numeric(line.list[[lnow]]$y)))}

            datmat.mono <- dataprep(y=y,x=x,ygrid=as.numeric(line.list[[l]]$y),xgrid=xgrid.cqf,
                                    bmat=Bmat,info=info,iyknots=iyknots,ydf=ydf,addxint=TRUE,yorder=yorder,
                                    yorth=FALSE,xorth=FALSE,Ysing=Ysing,e0mode=e0mode,returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,
                                    nxgrid=length(line.list[[l]]$x),nygrid=length(line.list[[l]]$y))

            lmin        <- c(lmin,min(unlist(datmat.mono$dedygrid)))
          }

        }

      }

      datmat.mono$dedygrid <- min(lmin)

    }

  }

  return(datmat.mono)

}

