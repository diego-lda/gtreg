
#' @title dataprep
#'
#' @description This is a function that takes Y and X and returns design matrices. It is useful to get recursive e/dedy basis funcitons.
#'
#' @param y This is the dependent variable, must be continous.
#' @param x These are that covariates used.
#' @param info
#' @param ygrid
#' @param xgrid
#' @param bmat
#' @param iyknots
#' @param addxint
#' @param ydf
#' @param yorder
#' @param yorth
#' @param xorth
#' @param Ysing
#' @param e0mode
#' @param mask
#' @param returnTZ
#' @param returnTZg
#' @param de0dygrid
#' @param extrapolate
#' @param plot.mode
#' @param delta
#' @param nxgrid
#' @param nygrid
#' @param Ys
#' @param YS
#' @param Xs
#' @param Xsgrid
#' @param sYgrid
#' @param SYgrid
#' @param TZ
#' @param tZ
#' @param TZ.pl
#' @param tZ.pl
#' @param TZgrid
#' @param tZgrid
#' @param TZgrid.pl
#' @param tZgrid.pl
#'
#' @return
#' @export
#'
#' @examples
dataprep <- function(y,x,info=NULL,ygrid=NULL,xgrid=NULL,
                     bmat=NULL,iyknots=NULL,addxint=T,ydf,yorder,yorth=FALSE,
                     xorth=FALSE,Ysing=FALSE,e0mode=F,mask=NULL,returnTZ=T,returnTZg=F,
                     de0dygrid=NULL,  # For now (used in density prediction) - should be done internally
                     extrapolate=T,plot.mode=F,delta=1,
                     nxgrid,nygrid,
                     Ys=NULL,YS=NULL,
                     Xs=NULL,Xsgrid=NULL,
                     sYgrid=NULL,SYgrid=NULL,
                     TZ=NULL,tZ=NULL,
                     TZ.pl=NULL,tZ.pl=NULL,
                     TZgrid=NULL,tZgrid=NULL,
                     TZgrid.pl=NULL,tZgrid.pl=NULL
){

  process.y <- T
  if(!is.null(YS) && !is.null(Ys)){
    process.y <- F
    nYS <- ncol(YS)
  }

  process.x <- T
  if(!is.null(Xs)){
    process.x <- F
    nXs <- ncol(Xs)
    if(!is.list(xgrid)){
      xgrid.now <- list()
      xgrid.now[[1]] <- xgrid
      xgrid <- list(xgrid=xgrid.now)
    }
    ifelse(is.null(xgrid),ngrids <- 1,ifelse(length(xgrid$xgrid[[1]])!=nxgrid, ngrids <- length(xgrid$xgrid[[1]][,1])/nxgrid, ngrids <- 1))
  }


  ifelse(length(xgrid)>0, gridx.mode <- T, gridx.mode <- F)
  ifelse(length(ygrid)>0, gridy.mode <- T, gridy.mode <- F)
  ifelse(!gridx.mode && !gridy.mode, grid.mode<-F, grid.mode <- T)


  # ======================================
  # y processing begins.
  # ======================================

  if(process.y){

    if(ydf > 0){

      if(e0mode){if(is.null(iyknots)){iyknots <- quantile(y,probs=seq(0,1,length=ydf))}}
      if(!e0mode){
        yabsmax <- max(abs(y))
        knoteps <- 1.e-5
        if(is.null(iyknots)){
          iyknots    <- qnorm(seq(knoteps,1-knoteps,length=ydf))
          iyknots[1] <- -delta*yabsmax
          iyknots[length(iyknots)] <- delta*yabsmax
        }
      }
      iyknots[[1]] <- plyr::round_any(iyknots[[1]],0.001,floor)
      iyknots[[length(iyknots)]] <- plyr::round_any(iyknots[[length(iyknots)]],0.001,ceiling)
      iyknots <- round(iyknots,3)
      names(iyknots) <- NULL

      y0<-y;
      if(plot.mode){ y<-ygrid }

      yknots <- expand.knots(iyknots,order=yorder)
      ys     <- SplineBasis(yknots,order=yorder)
      if(yorth) ys <- orthogonalize(ys)
      yS     <- integrate(ys)
      Ys     <- evaluate(object=ys,x=y)
      YS1    <- evaluate(object=yS,x=y)

      if(!Ysing){
        YS <- YS1
        YS[,1] <- y
        YS <- cbind(1,YS)
        Ys <- Ys[,-1]
        Ys <- cbind(1,Ys)
        nYS <- ncol(YS)
      }
      if(Ysing){
        YS <- YS1
        YS <- cbind(1,y,YS)
        Ys <- cbind(0,1,Ys)
        nYS <- ncol(YS)
      }

      nobs    <- length(y)

      if(gridy.mode){
        Ygrid <- ygrid
        #Prepare for loc.-scale extrapolation outside empirical Y support
        if(extrapolate == T){
          ext.dex.l <- which(ygrid<iyknots[1])
          ext.dex.u <- which(ygrid>iyknots[length(iyknots)])
          ext.dex   <- c(ext.dex.l,ext.dex.u)
        }

        #Ygrid   <- ygrid  #tmp patch rhs
        SYfunc  <- yS
        sYfunc  <- ys
        SYy     <- EvalBasis(object=SYfunc,x=c(ygrid,max(y0)))
        sYy     <- EvalBasis(object=sYfunc,x=ygrid)
        if(extrapolate==T){
          if(length(ext.dex.l)>0){
            SYy[ext.dex.l,] <- matrix(0,length(ext.dex.l),ncol(SYy))
            sYy[ext.dex.l,] <- matrix(0,length(ext.dex.l),ncol(SYy))
          }
          if(length(ext.dex.u)>0){
            SYy[ext.dex.u,] <- matrix(SYy[nrow(SYy),],length(ext.dex.u),ncol(SYy),byrow=TRUE)
            sYy[ext.dex.u,] <- matrix(0,length(ext.dex.u),ncol(sYy))
          }
        }
        SYy     <- SYy[1:length(ygrid),]

        #Ygrid <- ygrid   # This is bad

        #do the replacement
        if(!Ysing){
          SYgrid <- cbind(1,Ygrid,SYy[,-1])
          sYgrid <- cbind(0,1,sYy[,-1])
        }
        if(Ysing){
          SYgrid <- cbind(1,Ygrid,SYy)
          sYgrid <- cbind(0,1,sYy)
        }

      }

    }

    if(ydf==0){
      YS     <- cbind(1,y)
      Ys     <- cbind(0,1)
      nYS    <- ncol(YS)
      nobs   <- length(y)
      Ygrid  <- ygrid  #tmp patch rhs
      SYgrid <- cbind(1,Ygrid)
      sYgrid <- cbind(0,rep(1,length(ygrid)))
    }

    #in old style, Ys needed a 'zero' col as in:
    if(!Ysing && ydf>0){ Ys <-cbind(0,Ys) }  #so now it has it

  }

  # ======================================
  # x processing begins.
  # ======================================

  if(process.x){

    if(!is.data.frame(x)){x <- data.frame(x)}  #cast x to data.frame always

    #cat("xnss2 beginning x processing, checking for factors\n")
    #the length of a data.frame is the number of columns
    isfactor <- logical(length(x))
    for(i in seq_along(x)){
      isfactor[[i]] <- is.factor(x[[i]])
    }
    nofactors <- !any(isfactor) #make life easier

    #=======info is present, process it==========================
    #process through gX6 BEFORE model.matrix();(xm not yet defined)
    if(!is.null(info)){
      #cat("xnss2 has found info and is processing it \n")
      #info.look <<- info
      #x.look <<- x
      xnow <- gX6(X=x,info=info,orth=xorth)
      #xnow includes splined variables, no intercept, factor variables unexpanded
      xxnow <- data.frame(xnow[[1]],xnow)
      #xxnow.look <<- xxnow
    }
    else xnow <- x

    if(nofactors){
      #xxnow <- xxnow[,-1]
      Xs <- xnow
      xm <- x
    }

    #cat("xnss2 processing factors, if any\n")
    if(any(isfactor)){
      Xs <- model.matrix(formula(xxnow),data=xxnow)
      #Xsbd.look <<- Xs
      Xs <- Xs[,-1] #drop the intercept that is added by model.matrix()
      # Xs <- xm
      xorigformm <- data.frame(x[[1]],x)
      xm <- model.matrix(formula(xorigformm), data=xorigformm)
      xm <- xm[,-1]
      infom <- expandinfo(info, x)
      #xm <- x
    }
    #xm <- if( ncol(xm)>1) xm[, !apply(xm==0,2,all)] else xm[!apply(xm==0,2,all)]
    #xm <- if(!isdf | ncol(xm)>1) xm[, !apply(xm==0,2,all)] else xm[!apply(xm==0,2,all)]


    #Xs <- xnow  #basically density estimation
    #cat("xnss2 processing at 106\n")
    if(!is.matrix(Xs)) Xs <- as.matrix(Xs)
    if(addxint){
      Xs <- cbind(1,Xs)
      xm <- cbind(1,xm)
    } #add one intercept for the whole thing
    #Xs <- xnow
    nXs <- ncol(Xs)


    if(gridx.mode){

      if(is.null(info)){
        Xsgrid <- cbind(1,xgrid)
      }

      if(!is.null(info)){

        Xsgrid <- list()

        if(!is.list(xgrid)){
          nvars          <- 1
          xgrid.now      <- list()
          xgrid.now[[1]] <- xgrid
          xgrid          <- list(xgrid=xgrid.now)
        }
        else{ nvars <- length(xgrid$xgrid) }

        for(i in 1:nvars){
          ifelse(length(as.matrix(xgrid$xgrid[[i]]))!=nxgrid, ngrids <- length(xgrid$xgrid[[i]][,1])/nxgrid, ngrids <- 1)
          xgrid.subset <- split(as.data.frame(xgrid$xgrid[[i]]), rep(1:ngrids, each = nxgrid))
          Xsgrid[[i]]  <- array(0,c(nxgrid,nXs,ngrids))
          for(j in 1:ngrids){
            ifelse(addxint,
                   Xsgrid[[i]][,,j] <- cbind(1,as.matrix(gX6(X=as.data.frame(xgrid.subset[[j]]),info=info,orth=xorth))),
                   Xsgrid[[i]][,,j] <- as.matrix(gX6(X=as.data.frame(xgrid.subset[[j]]),info=info,orth=xorth)))
          }
        }

      }
    }

    if(!gridx.mode){ Xsgrid <- NULL }

    #so Xs and YS are relevant for TZ
    #and Xs and Ys for tZ
  }

  # ======================================
  # TZ forming.
  # ======================================

  #Stu Feldman
  #cat("line 101 calls to TZ1 in progress,\n")
  #cat("xnss2 processing at 119\n")
  #Xs.look <<- Xs

  if(returnTZ){
    TZ <- TZ1(X=Xs,Y=YS)
    tZ <- TZ1(X=Xs,Y=Ys)
    if(is.null(mask)){ TZ.pl <- NULL; tZ.pl <- NULL }
    if(!is.null(mask)){
      TZ.pl <- TZ[,which(mask==1)]
      tZ.pl <- tZ[,which(mask==1)]
    }
  }

  if(returnTZg){

    if(nxgrid==length(ygrid)){
      if(!gridx.mode){Xsgrid <- Xs}
      if(!gridy.mode){SYgrid <- YS; sYgrid <- Ys}
      TZgrid    <- list()
      tZgrid    <- list()
      TZgrid.pl <- list()
      tZgrid.pl <- list()
      for(i in 1:length(Xsgrid)){
        TZgrid[[i]]    <- array(0,c(nxgrid,nXs*nYS,ngrids))
        tZgrid[[i]]    <- array(0,c(nxgrid,nXs*nYS,ngrids))
        TZgrid.pl[[i]] <- array(0,c(nxgrid,nXs*nYS,ngrids))
        tZgrid.pl[[i]] <- array(0,c(nxgrid,nXs*nYS,ngrids))
        for(j in 1:ngrids){
          TZgrid[[i]][,,j] <- TZ1(X=Xsgrid[[i]][,,j],Y=SYgrid)
          tZgrid[[i]][,,j] <- TZ1(X=Xsgrid[[i]][,,j],Y=sYgrid)
          if(is.null(mask)){ TZgrid.pl <- NULL; tZgrid.pl <- NULL }
          if(!is.null(mask)){
            TZgrid.pl[[i]][,,j] <- TZgrid[[i]][,,j][,which(mask==1)]
            tZgrid.pl[[i]][,,j] <- tZgrid[[i]][,,j][,which(mask==1)]
          }
        }
      }
    }

    if(nxgrid!=length(ygrid)){
      if(!gridx.mode){Xsgrid <- Xs}
      if(!gridy.mode){SYgrid <- YS; sYgrid <- Ys}
      TZgrid <- list()
      tZgrid <- list()
      TZgrid.pl <- list()
      tZgrid.pl <- list()
      for(i in 1:length(Xsgrid)){
        TZgrid[[i]] <- array(0,c(nxgrid,nXs*nYS,nrow(Xsgrid),ngrids))
        tZgrid[[i]] <- array(0,c(nxgrid,nXs*nYS,nrow(Xsgrid),ngrids))
        TZgrid.pl[[i]] <- array(0,c(nxgrid,nXs*nYS,nrow(Xsgrid),ngrids))
        tZgrid.pl[[i]] <- array(0,c(nxgrid,nXs*nYS,nrow(Xsgrid),ngrids))
        #      TZgrid <- array(0,c(nrow(SYgrid),ncol(SYgrid)*ncol(Xsgrid),nrow(Xsgrid)))
        #      tZgrid <- array(0,c(nrow(SYgrid),ncol(SYgrid)*ncol(Xsgrid),nrow(Xsgrid)))
        if(!is.null(mask)){
          TZgrid.pl <- list()
          tZgrid.pl <- list()
        }
        for(i in 1:nrow(Xsgrid)){
          TZgrid[,,i] <- TZ2(X=Xsgrid[i,],Y=SYgrid)
          tZgrid[,,i] <- TZ2(X=Xsgrid[i,],Y=sYgrid)
          if(!is.null(mask)){
            TZgrid.pl[[i]] <- TZgrid[,which(mask==1),i]
            tZgrid.pl[[i]] <- tZgrid[,which(mask==1),i]
          }
        }
        if(is.null(mask)){ TZgrid.pl <- NULL; tZgrid.pl <- NULL }
      }
    }

  }

  # Get btarg
  muy          <- mean(YS[,2])
  sdy          <- sd(YS[,2])
  btarg        <- rep(0,len=nXs*nYS)
  btarg[1]     <- -muy/sdy
  btarg[nXs+1] <- 1/sdy

  # Storage
  ans          <- list()

  if(!is.null(bmat)){
    beta         <- Xs%*%bmat
    dedy         <- beta%*%t(Ys)
    e            <- beta%*%t(YS)
    ans$beta     <- beta
    ans$dedy     <- dedy
    ans$e        <- e
    if(grid.mode){
      ifelse(is.null(de0dygrid), de0dygrid <- 1, de0dygrid <- de0dygrid)
      if(!gridx.mode){Xsgrid <- Xs}
      betagrid <- list()
      egrid    <- list()
      dedygrid <- list()
      phiegrid <- list()
      fgrid    <- list()
      for(i in 1:length(Xsgrid)){
        betagrid[[i]] <- array(0,c(nxgrid,nYS,ngrids))
        egrid[[i]]    <- array(0,c(nxgrid,nygrid,ngrids))
        dedygrid[[i]] <- array(0,c(nxgrid,nygrid,ngrids))
        phiegrid[[i]] <- array(0,c(nxgrid,nygrid,ngrids))
        fgrid[[i]]    <- array(0,c(nxgrid,nygrid,ngrids))
        for(j in 1:ngrids){
          betagrid[[i]][,,j] <- Xsgrid[[i]][,,j]%*%bmat
          egrid[[i]][,,j]    <- betagrid[[i]][,,j]%*%t(SYgrid)
          dedygrid[[i]][,,j] <- betagrid[[i]][,,j]%*%t(sYgrid)
          phiegrid[[i]][,,j] <- dnorm2(egrid[[i]][,,j])
          fgrid[[i]][,,j]    <- phiegrid[[i]][,,j]*dedygrid[[i]][,,j]*de0dygrid
        }
      }
      ans$betagrid <- betagrid
      ans$egrid    <- egrid
      ans$dedygrid <- dedygrid
      ans$fgrid    <- fgrid
    }
  }

  if(returnTZ){
    ans$TZ       <- TZ
    ans$tZ       <- tZ
    ans$TZ.pl    <- TZ.pl
    ans$tZ.pl    <- tZ.pl
  }
  if(returnTZg){
    ans$TZgrid       <- TZgrid
    ans$tZgrid       <- tZgrid
    ans$TZgrid.pl    <- TZgrid.pl
    ans$tZgrid.pl    <- tZgrid.pl
  }
  if(grid.mode){
    ans$sYgrid <- sYgrid
    ans$SYgrid <- SYgrid
    ans$Xsgrid <- Xsgrid
  }
  ans$Xs       <- Xs
  ans$YS       <- YS
  ans$Ys       <- Ys
  ans$iyknots  <- iyknots
  ans$btarg    <- btarg

  return(ans)
}


#' @title TIC.function
#'
#' @description This function does some internal work
#'
#' @param res
#' @param TZ
#' @param tZ
#' @param Ginv
#' @param Adex
#' @param method
#'
#' @return
#' @export
#'
#' @examples
TIC.func <- function(res,TZ,tZ,Ginv=T,Adex=NULL,method="TIC"){

  if(length(Adex)==0){ Adex <- 1:ncol(TZ) }
  TZ       <- TZ[,Adex]
  tZ       <- tZ[,Adex]
  bmat     <- res$bmat[Adex]
  ifelse(length(Adex)==1, e <- TZ*bmat, e <- TZ%*%bmat)
  ifelse(length(Adex)==1, eta <- tZ*bmat, eta <- tZ%*%bmat)

  if(method=="TIC"){

    if(min(eta)>0){

      llf      <- sum(log(dnorm(e)*eta))
      gradmat  <- -TZ*as.vector(e)+tZ*as.vector(1/eta)
      grad2    <- t(gradmat)%*%gradmat
      hessmat  <- -(t(TZ)%*%TZ+t(tZ/as.vector(eta))%*%(tZ/as.vector(eta)))
      Hinv     <- NULL
      ifelse(length(Ginv)==0, try(Hinv <- solve(-hessmat,tol=1e-40), silent=T), Hinv <- ginv(-hessmat))

      bias <- sum(diag(Hinv%*%grad2))
      ans  <- -2*llf + 2*bias
    }

    if(min(eta)<=0){ ans <- 1e22 }
  }

  if(method=="BIC"){

    if(min(eta)>0){

      llf <- sum(log(dnorm(e)*eta))
      df  <- ncol(TZ)
      ans <- -2*llf+df*log(length(y))

    }

    if(min(eta)<=0){ ans <- 1e22 }

  }

  return(ans)

}



#' @title plot3d.gtr
#'
#' @description This Makes a 3D plot.
#'
#' @param y
#' @param x
#' @param ugrid
#' @param res
#' @param bmat
#' @param type
#' @param ydf
#' @param yorder
#' @param info
#' @param Ysing
#' @param ng.plot
#' @param ygrid
#' @param xgrid
#' @param dx
#' @param dy
#' @param FOV
#' @param alpha
#' @param ncolors
#' @param coltype
#' @param xlab
#' @param ylab
#' @param zlab
#' @param theta
#' @param phi
#' @param d
#' @param expand
#' @param r
#' @param liness
#' @param shade
#' @param blackgrid
#' @param projrq
#' @param BW
#' @param cex
#' @param ticktype
#' @param sizept
#' @param colpt
#' @param colines
#' @param lwd.liness
#' @param lwd.projrq
#' @param col.liness
#' @param col.projrq
#' @param col.surface
#' @param zlim
#' @param main
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
plot3d.gtr <- function(y,x,ugrid=seq(.1,.9,by=.1),res,bmat,type = "CDF",ydf=NULL,yorder=NULL,info = NULL, Ysing = FALSE, ng.plot = 30,
                       ygrid = NULL, xgrid = NULL, dx = 0, dy = 0, FOV = 1, alpha = 1, ncolors = 500, coltype = 1, xlab = "X",
                       ylab = "Y", zlab = "", theta = 70, phi = 30, d = 1, expand = 1, r = sqrt(3), liness = F, shade=NA,
                       blackgrid = NULL, projrq = F, BW = NULL, cex = 1, ticktype = "simple", sizept = 3, colpt = "blue", colines = "red",
                       lwd.liness=1,lwd.projrq=1,col.liness=1,col.projrq=1,col.surface="lightblue",zlim=NULL,main=NULL,delta=1){

  ifelse(length(xgrid)==0, gx <- seq(min(x)-dx, max(x)+dx, len = ng.plot), gx <- xgrid)
  ifelse(length(ygrid)==0, gy <- seq(min(y)-dy, max(y)+dy, len = ng.plot), gy <- ygrid)

  datmat    <- dataprep(y=as.numeric(y),x=x,xgrid=gx,ygrid=gy,bmat=bmat,info=info,iyknots=NULL,ydf=ydf,
                        addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,Ysing=Ysing,e0mod=F,delta=delta)
  z   <- datmat$egrid
  range(pnorm(z))

  if(type=="CDF"){

    z         <- datmat$egrid
    e         <- res$e
    line.list <- contourLines(gx,gy,pnorm(z),levels=ugrid)
    z[which(z<=range(e)[1])] <- range(e)[1]
    z[which(z>=range(e)[2])] <- range(e)[2]

    Fn        <- approxfun(e, (rank(e)-.5)/length(e));
    Fninv     <- approxfun((rank(e)-.5)/length(e), e);
    Fnz       <- matrix(Fn(z), ng.plot, ng.plot)

    res.p     <- persp(gx,gy,Fnz,col=col.surface,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,d=d,r=r,expand=expand,ticktype=ticktype,border=NA,shade=shade,zlim=c(0,1),main=main)

    if(liness){

      for(i in 1:length(ugrid)){
        lines(trans3d(x=line.list[[i]]$x, y=line.list[[i]]$y, z=matrix(ugrid[i],length(line.list[[i]]$y),1), pmat = res.p), col = col.liness,lwd=lwd.liness)
      }

    }

    if(projrq){

      for(i in 1:length(ugrid)){
        lines(trans3d(x=line.list[[i]]$x, y=line.list[[i]]$y, z= matrix(0,length(line.list[[i]]$y),1), pmat = res.p), col = col.projrq,lwd=lwd.projrq)
      }

    }

  }


  if(type=="PDF"){
    eta <- res.al$eta

    if(length(zlim)==0){ zlim <- range(datmat$fgrid, na.rm = TRUE) }
    res.p  <- persp(gx,gy,datmat$fgrid,col=col.surface,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,d=d,r = r,expand=expand,ticktype=ticktype,border=NA,shade=.5,zlim=zlim,main=main)

    for(i in round(seq(1,ng.plot,len=10))){
      lines(trans3d(x=rep(gx[i],ng.plot),y=gy,z=datmat$fgrid[i,],pmat=res.p), col = 1,lwd=min(2,.5+1.5*i/ng.plot))
    }
  }

}



#' @title QGM.check
#'
#' @description This checks whether the solution that GTR provides satisfies Monotonicity. It creates a grid for the continous variables and checks for each value of discrete values.
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
#'
#' @return
#' @export
#'
#' @examples
QGM.check <- function(y,x,res,xgrid.qgm=seq(min(x),max(x),len=101),ugrid,nxgrid=101,nygrid,ng.qgm=201,nXs,nYS,ydf,yorder,iyknots=NULL,
                      info=info,Ysing=FALSE,delta=1,easy=F,zeros=NULL){

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
                            Ysing=Ysing,e0mode=F,returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,nxgrid=nxgrid,nygrid=nygrid)
    #    lmin        <- beta2.check(bmat=Bmat,Xs=datmat.mono$Xsgrid,nXs=nXs,nYS=nYS)$b2min
    if(!is.list(datmat.mono$Xsgrid)){
      lmin <- beta2.check(bmat=Bmat,Xs=datmat.mono$Xsgrid,nXs=nXs,nYS=nYS)$b2min
    }

    if(is.list(datmat.mono$Xsgrid)){
      lmin <- NULL
      for(i in 1:length(datmat.mono$Xsgrid)){
        for(j in 1:dim(datmat.mono$Xsgrid[[i]])[3]){

          lmin <- c(b2min,beta2.check(bmat=Bmat,Xs=datmat.mono$Xsgrid[[i]][,,j],nXs=nXs,nYS=nYS)$b2min)
        }
      }
    }

    datmat.mono$dedygrid <- min(lmin)
  }

  if(nYS!=2 && nXs==2){
    datmat.mono <- dataprep(y=as.numeric(y),x=x,xgrid=xgrid.qgm,ygrid=ygrid.qgm,bmat=Bmat,
                            info=info,iyknots=iyknots,ydf=ydf,addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,
                            Ysing=Ysing,e0mode=F,returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,nxgrid=nxgrid,nygrid=nygrid)
    dBetaY <- datmat.mono$sYgrid%*%t(Bmat)
    lmin   <- cbind(1,matrix(c(min(x),max(x)),nr=2,nc=1))%*%t(dBetaY)
    datmat.mono$dedygrid <- min(lmin)
  }

  if(nYS!=2 && nXs>2){
    if(easy){
      datmat.mono <- dataprep(y=as.numeric(y),x=x,ygrid=ygrid.qgm,xgrid=xgrid.qgm,bmat=Bmat,
                              info=info,iyknots=iyknots,ydf=ydf,addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,
                              Ysing=Ysing,e0mode=F,returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,nxgrid=nxgrid,nygrid=nygrid)
      mdedy <- NULL
      if(length(xgrid.qgm[[2]])>0){
        for(i in 1:nrow(xgrid.qgm[[2]])){
          mdedy       <- c(mdedy, min(datmat.mono$dedygrid[[i]]))
        }
      }
      if(length(xgrid.qgm[[2]])==0){
        mdedy       <- min(unlist(datmat.mono$dedygrid))
      }
      datmat.mono$dedygrid <- min(mdedy)
    }

    if(!easy){

      ngrids <- length(xgrid.qgm$xgrid[[1]][,1])/nxgrid  # This is the number of replications of the first X grid - one for each grid value of other continuous X
      lmin   <- NULL
      i      <- 0


      while(i < length(xgrid.qgm$xgrid)  && (lmin > 0 || length(lmin) == 0)){   # This is a loop over grid values of discrete x components.

        i <- i+1
        j <- 0

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
                                    Ysing=Ysing,e0mode=F,returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,nxgrid=nxgrid,nygrid=nygrid)

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


            datmat.mono <- dataprep(y=y,x=x,ygrid=as.numeric(line.list[[l]]$y),xgrid=xgrid.cqf,
                                    bmat=Bmat,info=info,iyknots=iyknots,ydf=ydf,addxint=TRUE,yorder=yorder,
                                    yorth=FALSE,xorth=FALSE,Ysing=Ysing,e0mode=F,returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,
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


#' @title beta2.check
#'
#' @description This function checks what the smallest beta is as well as returning them.
#'
#' @param bmat
#' @param Xs
#' @param nXs
#' @param nYS
#'
#' @return
#' @export
#'
#' @examples
beta2.check <- function(bmat,Xs,nXs,nYS){

  beta  <- Xs%*%matrix(bmat,nr=nXs,nc=nYS)

  b2min  <- min(beta[,2])
  beta2  <- beta[,2]

  return(list(b2min=b2min,beta2=beta2))

}


#' @title cqf.func4
#'
#' @description This function is internal.
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
#'
#' @return
#' @export
#'
#' @examples
cqf.func4 <- function(y,x,xgrid,info,Ysing=FALSE,bmat,ng.plot,ydf,yorder,ugrid,nr,nc,delta=1){

  #y=y;x=x;xgrid=xx;info=info;Ysing=Ysing;bmat=res$bmat;ng.plot=ng.plot;ydf=ydf.dex[j];yorder=yorder.dex[i];ugrid=ugrid;nr=nr;nc=nc;delta=1


  if(yorder == 0 && nr>2){ # CQFs in closed-form for location-scale

    ygrid.p   <- seq(min(y),max(y),len=ng.plot)
    Bmat      <- matrix(bmat,nr=nr,nc=nc)
    datmatnow <- dataprep(y=ygrid.p,x=xgrid,ygrid=ygrid.p,xgrid=xgrid,bmat=Bmat,
                          info=info,iyknots=NULL,ydf=ydf,
                          addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,Ysing=Ysing,e0mod=F,
                          returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,#)
                          nxgrid=ng.plot,nygrid=ng.plot)

    beta      <- datmatnow$Xsgrid[[1]][,,1]%*%matrix(bmat,nr=nr,nc=nc)

    line.list <- list()
    nlevel    <- length(ugrid)
    lmin      <- NULL
    for(l in 1:nlevel){
      line.list[[l]] <- list(x=xgrid, y=-beta[,1]/as.vector(beta[,2]) + (1/as.vector(beta[,2]))*qnorm(ugrid[l]))
      lmin           <- c(lmin,min(beta[,2]))
    }

  }

  if(yorder > 0 || nr==2){

    #ng.plot   <- 201
    Bmat     <- matrix(bmat,nr=nr,nc=nc)
    ygrid    <- seq(min(y),max(y),len=ng.plot)

    if(nr==2){
      datmatnow <- dataprep(y=as.numeric(y),x=xgrid,ygrid=ygrid,xgrid=xgrid,bmat=Bmat,
                            info=info,iyknots=NULL,ydf=ydf,
                            addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,Ysing=Ysing,e0mod=F,
                            returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,#)
                            nxgrid=ng.plot,nygrid=ng.plot)
    }

    if(nr>2){

      del.up   <- 0
      del.down <- 0
      maa.up   <- 0
      maa.down <- 0
      Ymeas    <- max(y) - min(y)

      while(maa.down==0 || maa.up==0){

        datmatnow <- dataprep(y=as.numeric(y),x=xgrid,ygrid=ygrid,xgrid=xgrid,bmat=Bmat,
                              info=info,iyknots=NULL,ydf=ydf,
                              addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,Ysing=Ysing,e0mod=F,
                              returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,nxgrid=ng.plot,nygrid=ng.plot)
        if(maa.up==0){
          aa.up <- NULL
          for(ll in 1:length(xgrid)){
            aa.up <- c(aa.up,length((which(pnorm(datmatnow$egrid[[1]][,,1])[ll,]>.99))))
          }
          maa.up   <- min(aa.up)
        }
        if(maa.down==0){
          aa.down <- NULL
          for(ll in 1:length(xgrid)){
            aa.down <- c(aa.down,length((which(pnorm(datmatnow$egrid[[1]][,,1])[ll,]<.01))))
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

        line.list.now <- contourLines(xgrid,ygrid,pnorm(datmatnow$egrid[[1]][,,1]),levels=ugrid)
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

    line.list <- contourLines(xgrid,ygrid,pnorm(datmatnow$egrid[[1]][,,1]),levels=ugrid)
    nlevel    <- length(line.list)
    lmin      <- NULL
    for(l in 1:nlevel){
      datmatnow <- dataprep(y=y,x=xgrid,xgrid=line.list[[l]]$x,ygrid=line.list[[l]]$y,
                            bmat=Bmat,info=info,iyknots=NULL,ydf=ydf,yorder=yorder,
                            addxint=TRUE,yorth=FALSE,xorth=FALSE,Ysing=Ysing,e0mode=F,
                            returnTZ=F,returnTZg=F,plot.mode=T,delta=delta,#)
                            nxgrid=length(line.list[[l]]$x),nygrid=length(line.list[[l]]$y))
      lmin   <- c(lmin, min(unlist(datmatnow$dedygrid)))
    }

    dedy.min <- min(lmin)

  }


  ans    <- list(line.list=line.list,lmin=lmin,nlevel=nlevel)


  return(ans)

}



#' @title TZ2
#'
#' @description This function is internal.
#'
#' @param X
#' @param Y
#'
#' @return
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



#' @title logspace
#'
#' @description This function is internal.
#'
#' @param a
#' @param b
#' @param n
#'
#' @return
#' @export
#'
#' @examples
logspace   <- function( a, b, n){exp(log(10)*seq(a, b, length.out=n))}




#' @title data.info
#'
#' @description This function is internal. It generates the data info for a variable.
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
#'
#' @return
#' @export
#'
#' @examples
data.info <- function(x, X.type, xdf, xknots=NULL, coord.bare=NULL, coord.spline=NULL, coord.tensor=NULL, nxgrid=101, gridx.cont=F){

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
    if(xknots.missing){ xknots <- seq(min(x[,i]),max(x[,i]),length.out=xdf[i]) }
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
          if(xknots.missing){ xknots <- seq(min(x[,tdex]),max(x[,tdex]),length.out=xdf[tdex]) }
          info[[linfo]][[i]] <- list(type = "spline", knots = xknots, coords = tdex, center = mean(x[,tdex]))
        }

      }

      info[[linfo]]$type <- "tensor"

    }

  }

  return(info)
}



#' @title grids.func
#'
#' @description This function is internal. It generates the grids for a variable.
#'
#' @param x
#' @param X.type
#' @param nxgrid
#' @param gridx.cont
#'
#' @return
#' @export
#'
#' @examples
grids.func <- function(x, X.type, nxgrid=101, gridx.cont=F){

  x         <- data.frame(x)
  ncov      <- NCOL(x)
  X.c.dex   <- which(X.type=="continuous")
  ncov.fac  <- ncov-length(X.c.dex)   # Number of discrete covariates
  ncov.cont <- ncov-ncov.fac          # Number of continuous covariates
  xgrid     <- NULL  # IS THIS NEEDED?
  #xvals     <- NULL

  if(gridx.cont==T){
    xgrid <- list()
    for(i in 1:ncov.cont){
      xgrid[[i]] <- seq(min(x[X.c.dex[i]]),max(x[X.c.dex[i]]),len=nxgrid)
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
      xvals[[i]] <- sort(unique(unlist(x[X.fac.dex[i]])))
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


#' @title iyknots.func
#'
#' @description This function is internal. It generates the knots for the splines.
#'
#' @param y
#' @param e0mode
#' @param ydf
#' @param yorder
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
iyknots.func <- function(y,e0mode=F,ydf,yorder,delta=1){
  if(e0mode){ iyknots <- quantile(y,probs=seq(0,1,length=ydf)) }
  if(!e0mode){
    yabsmax <- max(abs(y))
    knoteps <- 1.e-5
    iyknots    <- qnorm(seq(knoteps,1-knoteps,length=ydf))
    iyknots[1] <- -delta*yabsmax
    iyknots[length(iyknots)] <- delta*yabsmax
  }
  # iyknots[[1]] <- plyr::round_any(iyknots[[1]],0.001,floor)
  # iyknots[[length(iyknots)]] <- plyr::round_any(iyknots[[length(iyknots)]],0.001,ceiling)
  # iyknots <- round(iyknots,3)
  # names(iyknots) <- NULL
  #
  # yknots <- expand.knots(iyknots,order=yorder)

  return(iyknots)
}
