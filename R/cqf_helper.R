#' Conditional Quantile Function Helper
#'
#' @description This function runs some stuff necessary to get the conditional quantile functions.
#'
#' @param y
#' @param x
#' @param xgrid
#' @param info
#' @param Ysing
#' @param bmat
#' @param ng.plot
#' @param ydf
#' @param yorder
#' @param ugrid
#' @param nr
#' @param nc
#' @param delta
#' @param e0mode
#'
#' @return
#' @export
#'
#' @examples
cqf_helper <- function(y,x,xgrid,info,Ysing=FALSE,bmat,ng.plot,ydf,yorder,ugrid,nr,nc,delta=1,e0mode=F){

  nxgrid   <- NROW(xgrid)
  exit <- F
  try <- F
  nlevel <- NULL

  while(!exit){
    exit <- T

    if(yorder == 0 && (nr>2  || (nr ==2 && try))){ # CQFs in closed-form for location-scale

      ygrid.p   <- seq(min(y),max(y),len=ng.plot)
      Bmat      <- matrix(bmat,nr=nr,nc=nc)
      datmatnow <- data_prep(y=ygrid.p,x=xgrid,ygrid=ygrid.p,xgrid=xgrid,bmat=Bmat,
                            info=info,y_knots=NULL,ydf=ydf,
                            addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,Ysing=Ysing,
                            returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,#)
                            nxgrid=nxgrid,nygrid=ng.plot,e0mode=e0mode)

      beta      <- datmatnow$Xsgrid[[1]][,,1]%*%matrix(bmat,nr=nr,nc=nc)

      line.list <- list()
      nlevel    <- length(ugrid)
      lmin      <- NULL
      for(l in 1:nlevel){
        line.list[[l]] <- list(x=as.matrix(xgrid)[,1], y=-as.vector(beta[,1])/as.vector(beta[,2]) + (1/as.vector(beta[,2]))*qnorm(ugrid[l]))
        lmin           <- c(lmin,min(beta[,2]))
      }

    }

    if(length(nlevel)==0 && (yorder > 0 || nr==2)){

      Bmat     <- matrix(bmat,nr=nr,nc=nc)
      ygrid    <- seq(min(y),max(y),len=ng.plot)

      if(nr==2){
        datmatnow <- data_prep(y=as.numeric(y),x=xgrid,ygrid=ygrid,xgrid=xgrid,bmat=Bmat,
                              info=info,y_knots=NULL,ydf=ydf,
                              addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,Ysing=Ysing,
                              returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,#)
                              nxgrid=nxgrid,nygrid=ng.plot,e0mode=e0mode)
      }

      if(nr>2){

        del.up   <- 0
        del.down <- 0
        maa.up   <- 0
        maa.down <- 0
        Ymeas    <- max(y) - min(y)

        while(maa.down==0 || maa.up==0){

          datmatnow <- data_prep(y=as.numeric(y),x=xgrid,ygrid=ygrid,xgrid=xgrid,bmat=Bmat,
                                info=info,y_knots=NULL,ydf=ydf,
                                addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,Ysing=Ysing,
                                returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,nxgrid=nxgrid,nygrid=ng.plot,e0mode=e0mode)

          eps <- 1e-2

          if(maa.up==0){
            aa.up <- NULL
            #for(ll in 1:length(xgrid)){
            for(ll in 1:NROW(xgrid)){
              aa.up <- c(aa.up,length((which(pnorm(datmatnow$egrid[[1]][,,1])[ll,]>1-eps))))
            }
            maa.up   <- min(aa.up)
          }
          if(maa.down==0){
            aa.down <- NULL
            #for(ll in 1:length(xgrid)){
            for(ll in 1:NROW(xgrid)){
              aa.down <- c(aa.down,length((which(pnorm(datmatnow$egrid[[1]][,,1])[ll,]<eps))))
            }
            maa.down <- min(aa.down)
          }
          if(maa.up==0){
            del.up    <- del.up+.5
            if(maa.down!=0){ ygrid <- seq(min(ygrid),max(ygrid)+del.up,len=ng.plot) }
            if(maa.down==0){
              del.down <- del.down+.5
              ygrid    <- seq(min(ygrid)-del.down,max(ygrid)+del.up,len=ng.plot)
            }
          }
          if(maa.down==0 && maa.up!=0){
            del.down <- del.down+.5
            ygrid    <- seq(min(ygrid)-del.down,max(ygrid),len=ng.plot)
          }

          line.list.now <- contourLines(xgrid[,1],ygrid,pnorm(datmatnow$egrid[[1]][,,1]),levels=ugrid)
          nlevel        <- length(line.list.now)
          umin.lev      <- NULL
          umax.lev      <- NULL

          for(l in 1:nlevel){
            if(line.list.now[[l]]$level==ugrid[1]){ umin.lev <- c(umin.lev, l) }
            if(line.list.now[[l]]$level==ugrid[length(ugrid)]){ umax.lev <- c(umax.lev, l) }
          }
          nlev.min <- length(umin.lev)
          nlev.max <- length(umax.lev)

          Ymeasnow <- max(ygrid) - min(ygrid)
          if(nlev.min==1 || Ymeasnow > 2*Ymeas){ maa.down <- 1 }
          if(nlev.max==1 || Ymeasnow > 2*Ymeas){ maa.up <- 1 }
          print(c(maa.up,maa.down))
          print(range(ygrid))

        }

      }

      line.list <- contourLines(as.matrix(xgrid)[,1],ygrid,pnorm(datmatnow$egrid[[1]][,,1]),levels=ugrid)
      #line.list <- contourLines(as.matrix(xgrid)[,1],ygrid,datmatnow$egrid[[1]][,,1],levels=qnorm(ugrid))
      nlevel    <- length(line.list)
      lmin      <- NULL

      ### ADDED (TO RETURN CQF AT SPECIFIED GRID POINTS)
      lev.dex <- rep(NA, length(ugrid))
      for(l in 1:length(ugrid)){
        count<-0
        for(ll in 1:nlevel){
          if(ugrid[l]==line.list[[ll]]$level && count!=1){
            count<-count+1
            lev.dex[l] <- ll
          }
        }
      }
      udex<-0
      line.list.now <- list()
      for(l in lev.dex){
        udex <- udex + 1
        qfunc <- splinefun(line.list[[l]]$x, line.list[[l]]$y, method="natural")
        line.list.now[[udex]] <- list(x=as.matrix(xgrid)[,1], y=qfunc(as.matrix(xgrid)[,1]))
      }

      line.list <- line.list.now
      nlevel <- length(ugrid)
      ### END MODIF

      lmin <- foreach(l = 1:nlevel, .combine="c",.packages=c("orthogonalsplinebasis", "splines","doParallel", "gtreg")) %dopar% {
        xgridnow <- line.list[[l]]$x
        if(NCOL(xgrid)>1){
          for(iii in 2:(NCOL(xgrid))){
            xgridnow <- cbind(xgridnow,xgrid[1,iii])
          }
        }

        datmatnow <- data_prep(y=y,x=xgrid,xgrid=xgridnow,ygrid=line.list[[l]]$y,
                              bmat=Bmat,info=info,y_knots=NULL,ydf=ydf,yorder=yorder,
                              addxint=TRUE,yorth=FALSE,xorth=FALSE,Ysing=Ysing,e0mode=e0mode,
                              returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,
                              nxgrid=length(line.list[[l]]$x),nygrid=length(line.list[[l]]$y))
        ans <- min(unlist(datmatnow$dedygrid))
        return(ans)
      }


    }
    if(nlevel != length(ugrid) && yorder == 0 && nr ==2){
      exit <- F
      try <- T
    }
    print(nlevel)
  }

  ans    <- list(line.list=line.list,lmin=lmin,nlevel=nlevel)

  return(ans)

}
