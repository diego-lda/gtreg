#' GTR Adaptive Lasso
#'
#' @description This is a wrap-up function for adaptive lasso GTR/spline-spline. It allows for a Lasso first-step.
#'
#' @param y
#' @param x
#' @param X.type
#' @param gam
#' @param lam.vec
#' @param xdf
#' @param ydf
#' @param yorder
#' @param nyg.max
#' @param Ysing
#' @param iyknots
#' @param ng.qgm
#' @param addxint
#' @param info
#' @param yorth
#' @param xorth
#' @param reltol
#' @param feastol
#' @param abstol
#' @param tol.res
#' @param cval
#' @param lb
#' @param threshold
#' @param e0mode
#' @param unconstrained
#' @param constrained
#' @param AL
#' @param doprimal
#' @param fac
#' @param maxit
#' @param model
#' @param bounded
#' @param beta2
#' @param Cbound
#' @param ugrid
#' @param gam.grid
#' @param algor
#' @param coord.bare
#' @param coord.spline
#' @param coord.tensor
#' @param easy
#' @param method
#' @param delta.ok
#' @param parallel
#'
#' @return
#' @export
#'
#' @examples
gtr.al <- function(y,x,X.type=NULL,gam=0,lam.vec=rep(1,6),xdf=0,ydf=4,yorder=4,nyg.max=101,Ysing=FALSE,iyknots=NULL,ng.qgm=101,
                   addxint=T,info=NULL,yorth=FALSE,xorth=FALSE,reltol=1.e-03,feastol=1.e-03,abstol=1.e-03,tol.res=1e-1,cval=1e-1,lb=1,threshold=1e-5,
                   e0mode=F,unconstrained=T,constrained=T,AL=F,doprimal=F,fac=1,maxit=90,model=NULL,bounded=F,beta2=F,Cbound=Inf,ugrid=seq(.1,.9,by=.1),
                   gam.grid=logspace(-3,3,20),algor="ECOS",coord.bare=NULL,coord.spline=NULL,coord.tensor=NULL,easy=T,method="BIC",delta.ok=F,parallel=F

){
  objval.bic  <- NULL
  res.al.bic  <- NULL
  mdedy       <- NULL
  status      <- NULL

  if(ydf==0 || xdf==0){ cval <- 0 }

  if(Cbound==Inf){ bounded = F }

  if(length(model)==0){

    nobs       <- length(y)
    res0       <- NULL
    res.sol    <- NULL
    res.al     <- NULL
    BIC.vec    <- NULL
    BIC.vec2   <- NULL
    fac.now2   <- NULL
    gam.al     <- NULL
    xknots     <- NULL
    cvg        <- NULL
    b2min.now  <- -1
    nyg.now    <- 0
    if(!beta2){ b2min <- 1 }

    if(ydf==0){yorder<-2} # tmp patch ss
    ifelse(Ysing || ydf==0, nYS <- ydf+yorder, nYS <- ydf+yorder-1)
    # Generate X data types if needed
    if(is.null(X.type)){ X.type <- rep("continuous",NCOL(x)) }
    # Check that each spline component as xdf defined:
    if(length(coord.spline)>0){stopifnot(length(xdf)==length(coord.spline))}

    xgrid.qgm  <- grids.func(x=x, X.type=X.type, nxgrid=ng.qgm, gridx.cont=T)
    ygrid.qgm  <- seq(min(y),max(y),len=ng.qgm)

    if(length(info)==0){ info <- data.info(x=x, X.type=X.type, xdf=xdf, coord.bare=coord.bare, coord.spline=coord.spline, coord.tensor=coord.tensor, nxgrid=ng.qgm, delta.ok=delta.ok) }
    if(dim(x)[2] == 1 && all(x[,1] == 1)){
      datmat     <- dataprep(y=as.numeric(y),x=x,ygrid=ygrid.qgm,
                             Xs=x,nxgrid=0,
                             iyknots=iyknots,ydf=ydf,addxint=addxint,yorder=yorder,
                             yorth=yorth,xorth=xorth,Ysing=Ysing,e0mode=e0mode,nygrid=ng.qgm)
    }else{
      datmat     <- dataprep(y=as.numeric(y),x=x,ygrid=ygrid.qgm,xgrid=xgrid.qgm,
                             info=info,iyknots=iyknots,ydf=ydf,addxint=addxint,yorder=yorder,
                             yorth=yorth,xorth=xorth,Ysing=Ysing,e0mode=e0mode,nxgrid=ng.qgm,nygrid=ng.qgm)
    }

    TYX        <- datmat$TZ;
    tYX        <- datmat$tZ
    Xs         <- datmat$Xs
    nXs        <- ncol(datmat$Xs);
    nYS        <- ncol(datmat$YS)
    Xsgrid.qgm <- datmat$Xsgrid
    sYgrid.qgm <- datmat$sYgrid
    rm(datmat)


    # ===============
    # Estimation
    # ===============

    dedy.min    <- -1
    cvg.mono    <- NULL
    nyg         <- -1/2
    #easy        <- F
    if(mean(xdf) == 0){ easy <- T }

    while( dedy.min <= .Machine$double.eps && nyg < nyg.max ){

      if(nyg>201 && easy==F && xdf !=0){ easy <- T; nyg <- 1}
      nyg     <- 2*nyg + 1
      if(nyg==1){ nyg <- 3 }; print(paste("nyg now =",nyg))
      if(nyg==0){
        Xsgrid <- Xs; #rm(Xs);
        sYgrid <- NULL
        if(dim(x)[2] == 1 && all(x[,1] == 1)){sYgrid <- sYgrid.qgm}
      }

      # Generate grid of Ys for QGM constraints
      if(nyg>0){
        nygnow <- 0+nyg
        xgrid  <- grids.func(x=x, X.type=X.type, nxgrid=nygnow, gridx.cont=T)
        ygrid  <- seq(min(y),max(y),len=nygnow)
        if(dim(x)[2] == 1 && all(x[,1] == 1) ){
          datmat <- dataprep(y=as.numeric(y),x=x,ygrid=ygrid,
                             Xs=x,nxgrid=0,
                             iyknots=iyknots,ydf=ydf,addxint=addxint,yorder=yorder,
                             yorth=yorth,xorth=xorth,Ysing=Ysing,e0mode=e0mode)
        }else{
          datmat <- dataprep(y=as.numeric(y),x=x,ygrid=ygrid,xgrid=xgrid,
                             info=info,iyknots=iyknots,ydf=ydf,addxint=addxint,yorder=yorder,
                             yorth=yorth,xorth=xorth,Ysing=Ysing,e0mode=e0mode,nxgrid=nygnow)
        }

        Xsgrid <- datmat$Xsgrid
        sYgrid <- datmat$sYgrid
        rm(datmat)
      }

      ans <- inner.gtr.c( nyg=nyg,TYX=TYX,tYX=tYX,gam=gam,reltol=reltol,feastol=feastol,abstol=abstol,
                          doprimal=doprimal,lam.vec=lam.vec,
                          y=y,x=x,info=info,yorder=yorder,ydf=ydf,
                          Ysing=Ysing,nYS=nYS,nXs=nXs,
                          xgrid.qgm=xgrid.qgm,
                          Xs=Xsgrid,sYgrid=sYgrid,ugrid=ugrid,
                          Xs.qgm=Xsgrid.qgm,sYgrid.qgm=sYgrid.qgm,
                          maxit=maxit,ng.qgm=ng.qgm,dedy.min=dedy.min,
                          tol.res=reltol,algor=algor,easy=easy,e0mode=e0mode )

      dedy.min    <- ifelse(length(ans$res.sol$eta)>0,min(ans$res.sol$eta,ans$dedy.min),ans$dedy.min)
      tol.res     <- ans$tol.res
      cvg.mono    <- ans$cvg.mono
      algor       <- ans$algor

      if(constrained == F){ nyg <- nyg.max+1 }   # This is to exit the loop with increasing number of QGM constraints

    }

    if(parallel){
      res.sol <- list(bmat=ans$res.sol$bmat, eta=ans$res.sol$eta)
      rm(ans)
    }
    if(!parallel){
      res.sol  <- ans$res.sol;
      nyg.now  <- ans$nyg.now;
      cvg      <- ans$cvg;
      fac.now2 <- ans$fac.now2;
    }

    if(dim(x)[2] == 1 && all(x[,1] == 1)){
      datmat.now <- dataprep(y=as.numeric(y),x=x,#xgrid=xgrid,ygrid=ygrid,
                             iyknots=iyknots,ydf=ydf,
                             Xs=x,nxgrid=0,
                             addxint=addxint,yorder=yorder,yorth=FALSE,xorth=FALSE,
                             Ysing=Ysing,e0mode=e0mode,returnTZ=T,returnTZg=F)
    }
    else{
      datmat.now <- dataprep(y=as.numeric(y),x=x,#xgrid=xgrid,ygrid=ygrid,
                             info=info,iyknots=iyknots,ydf=ydf,
                             addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,
                             Ysing=Ysing,e0mode=e0mode,returnTZ=T,returnTZg=F)
    }

    ################ STORE FIRST STEP RESULTS ###########################

    if( !parallel && (length(res.sol)>0 && min(res.sol$eta)>.Machine$double.eps )
        || parallel && (length(res.sol$bmat)>0 && min(res.sol$eta)>.Machine$double.eps)){
      nXs       <- ncol(datmat.now$Xs)
      if(gam==0){ BIC.vec   <- TIC.func(res=res.sol,TZ=datmat.now$TZ,tZ=datmat.now$tZ,method=method,Ginv=T) }
      if(gam>.Machine$double.eps){
        Adex <- which(abs(res.sol$bmat)>=threshold)
        BIC.vec <- TIC.func(res=res.sol,TZ=datmat.now$TZ,tZ=datmat.now$tZ,Adex=Adex,method=method,Ginv=T)
      }
      # BIC.vec   <- TIC.func(res=res.sol,TZ=datmat.now$TZ,tZ=datmat.now$tZ,method=method,Ginv=T)
      # ifelse(constrained && length(cvg)>0, nyg <- nyg.now, nyg <- nyg.now <- 0 )
    }
    rm(datmat.now)
  }


  ################### BEGIN ADAPTIVE LASSO ##########################

  # ===============================================
  # Retrieve info if gtr.al() was previously run separately without adaptive lasso

  if(length(model)>0){

    res.sol   <- model$res
    res.al    <- NULL
    gam.al    <- NULL
    BIC.vec   <- model$BIC
    BIC.vec2  <- NULL
    nXs       <- model$nr
    nYS       <- model$nc
    fac.now2  <- model$fac
    xknots    <- model$xknots
    nyg.now   <- model$nyg
    # nyg      <- model$nyg
    # nyg.now  <- nyg
    #dedymin.u <- model$dedymin.u
    cvg       <- model$cvg
    Xs        <- model$Xs
    #algor     <- mod$algor

  }

  # ===============================================

  if( AL && min(res.sol$eta)>0 ){

    doprimal <- T
    nglam    <- length(gam.grid)

    # Adaptive Lasso
    ans0      <- list()

    # Generate X data types if needed
    if(is.null(X.type)){ X.type <- rep("continuous",NCOL(x)) }
    xgrid.qgm  <- grids.func(x=x, X.type=X.type, nxgrid=ng.qgm, gridx.cont=T)
    ygrid.qgm  <- seq(min(y),max(y),len=ng.qgm)
    if(length(info)==0){ info <- data.info(x=x, X.type=X.type, xdf=xdf, coord.bare=coord.bare, coord.spline=coord.spline, coord.tensor=coord.tensor, nxgrid=ng.qgm, delta.ok=delta.ok) }
    datmat     <- dataprep(y=as.numeric(y),x=x,ygrid=ygrid.qgm,xgrid=xgrid.qgm,bmat=matrix(res.sol$bmat,nr=nXs,nc=nYS),
                           info=info,iyknots=iyknots,ydf=ydf,addxint=addxint,yorder=yorder,
                           yorth=yorth,xorth=xorth,Ysing=Ysing,e0mode=e0mode,nxgrid=ng.qgm,nygrid=ng.qgm)
    TYX        <- datmat$TZ;
    tYX        <- datmat$tZ
    Xs         <- datmat$Xs
    Xsgrid.qgm <- datmat$Xsgrid
    sYgrid.qgm <- datmat$sYgrid
    #rm(datmat)


    # ===============
    # Estimation
    # ===============
    ans0  <- list()
    ggnow <- 0

    for(gg in 1:nglam){
      print( paste("gg =",gg) )

      dedy.min    <- -1
      nyg         <- -1/2
      if(mean(xdf) == 0){ easy <- T }

      while( dedy.min <= .Machine$double.eps && nyg < nyg.max ){

        if(nyg>201 && easy==F && xdf !=0){ easy <- T; nyg <- 1}
        nyg     <- 2*nyg + 1
        if(nyg==1){ nyg <- 3 }; print(paste("nyg now =",nyg))
        if(nyg==0){
          Xsgrid <- Xs; #if(gg==5){rm(Xs)};
          sYgrid <- NULL
        }

        # Generate grid of Ys for QGM constraints
        if(nyg>0){
          xgrid  <- grids.func(x=x, X.type=X.type, nxgrid=nyg, gridx.cont=T)
          ygrid  <- seq(min(y),max(y),len=nyg)
          datmat <- dataprep(y=as.numeric(y),x=x,ygrid=ygrid,xgrid=xgrid,
                             info=info,iyknots=iyknots,ydf=ydf,addxint=addxint,yorder=yorder,
                             yorth=yorth,xorth=xorth,Ysing=Ysing,e0mode=e0mode,nxgrid=nyg)
          Xsgrid <- datmat$Xsgrid
          sYgrid <- datmat$sYgrid
          rm(datmat)
        }

        weights <- as.vector(abs(res.sol$bmat))^1
        zeros   <- which(weights<threshold)
        weights <- 1/as.numeric(weights)
        TYXnow  <- TYX
        tYXnow  <- tYX
        if(length(zeros>0)){
          weights <- weights[-zeros]
          TYXnow  <- TYX[,-zeros]
          tYXnow  <- tYX[,-zeros]
        }

        ans     <- inner.gtr.c( nyg=nyg,TYX=TYXnow,tYX=tYXnow,gam=gam.grid[gg],weights=weights,zeros=zeros,
                                reltol=reltol,feastol=feastol,abstol=abstol,
                                doprimal=doprimal,lam.vec=lam.vec,
                                y=y,x=x,info=info,yorder=yorder,ydf=ydf,
                                Ysing=Ysing,nYS=nYS,nXs=nXs,
                                xgrid.qgm=xgrid.qgm,
                                Xs=Xsgrid,sYgrid=sYgrid,ugrid=ugrid,
                                Xs.qgm=Xsgrid.qgm,sYgrid.qgm=sYgrid.qgm,
                                maxit=maxit,ng.qgm=ng.qgm,dedy.min=dedy.min,
                                tol.res=reltol,algor=algor,easy=easy,threshold=threshold,e0mode=e0mode )

        dedy.min <- ans$dedy.min

        if(constrained == F){ nyg <- nyg.max+1 }   # This is to exit the loop with increasing number of QGM constraints

      }

      if(length(ans$res.sol)>0){
        ggnow         <- ggnow + 1
        ans0[[ggnow]] <- ans$res.sol;
        #Adex          <- which(abs(ans0[[ggnow]]$bmat)>threshold)
        #BIC.now       <- TIC.func(res=ans0[[ggnow]],TZ=TYX,tZ=tYX,Adex=Adex,method=method)
        if(length(zeros)>0){
          TYXnow <- TYX[,-zeros]
          tYXnow <- tYX[,-zeros]
        }
        BIC.now       <- TIC.func(res=ans0[[ggnow]],TZ=TYXnow,tZ=tYXnow,method=method,Ginv=T)
        objval.bic    <- c(objval.bic, BIC.now)
      }
    }

    select.bic <- which.min(objval.bic)
    res.al     <- ans0[[select.bic]]
    BIC.vec2   <- objval.bic[select.bic]
    gam.al     <- gam.grid[select.bic]

    if(length(zeros)>0){
      M            <- matrix(0,nXs*nYS,nXs*nYS-length(zeros))
      M[-zeros,]   <- diag(1,nXs*nYS-length(zeros),nXs*nYS-length(zeros))
      res.al$bmat  <- M%*%res.al$bmat
    }

  }

  if(!parallel){
    res    <- list(bmat=res.sol$bmat,e=res.sol$e,eta=res.sol$eta)
    res.al <- list(bmat=res.al$bmat,e=res.al$e,eta=res.al$eta)
    if(length(res.sol)==0){ nXs <- NULL }
  }

  ifelse(parallel, ans    <- list(BIC=BIC.vec),
         ans    <- list(res=res,res.al=res.al,gam.al=gam.al,BIC=BIC.vec,BIC.al=BIC.vec2,
                        nr=nXs,nc=nYS,fac=fac.now2,xknots=xknots,nyg=nyg.now,cvg=cvg,
                        mdedy=mdedy,objval.bic=objval.bic,algor=algor,status=res.sol$result$status,
                        X.type=X.type,coord.bare=coord.bare,coord.spline=coord.spline,coord.tensor=coord.tensor))


  return(ans)

}
