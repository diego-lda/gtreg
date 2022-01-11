#' Dataprep
#'
#' @description A function that takes Y and X and returns design matrices
#'
#' @param y
#' @param x
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
#' @return The design matrices for Y and X.
#' @export
#'
#' @examples
data_prep <- function(y,x,info=NULL,ygrid=NULL,xgrid=NULL,
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
){  # if extrapolate=T location-scale extrapolation outside empirical support of Y

  nobs    <- length(y)

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
    ngrids <- 1
    if(!is.null(xgrid) && nxgrid!=0){
      if(length(xgrid$xgrid[[1]])!=nxgrid){ ngrids <- length(xgrid$xgrid[[1]][,1])/nxgrid }
    }
  }


  ifelse(length(xgrid)>0, gridx.mode <- T, gridx.mode <- F)
  ifelse(length(ygrid)>0, gridy.mode <- T, gridy.mode <- F)
  ifelse(!gridx.mode && !gridy.mode, grid.mode<-F, grid.mode <- T)


  # ======================================
  # y processing begins.
  # ======================================

  if(process.y){

    if(ydf > 0){

      if(is.null(iyknots)){iyknots <- quantile(y,probs=seq(0,1,length=ydf))}

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

      if(gridy.mode){
        Ygrid <- ygrid

        if(extrapolate == T){
          ext.dex.l <- which(ygrid<iyknots[1])
          ext.dex.u <- which(ygrid>iyknots[length(iyknots)])
          ext.dex   <- c(ext.dex.l,ext.dex.u)
        }

        SYfunc  <- yS
        sYfunc  <- ys
        SYy     <- eval_basis(object=SYfunc,x=c(ygrid,max(y0)))
        sYy     <- eval_basis(object=sYfunc,x=ygrid)
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
      Ys     <- matrix(c(0,1),nr=nobs,nc=2,byrow=T)
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

        # if(!("factor" %in% X.type)){
        #     ngrids <- length(xgrid$xgrid[,1])/nxgrid
        #     xgrid.subset <- split(as.data.frame(xgrid$xgrid), rep(1:ngrids, each = nxgrid))
        #     nXs <- length(coord.bare)+length(coord.spline)*(xdf+1)+length(coord.tensor)*length(coord.spline%in%coord.tensor)*(xdf+1) #Fix that and only for one tensor for now
        #     if(addxint){ nXs <- nXs + 1 }
        #     Xsgrid <- array(0,c(nxgrid,nXs,ngrids))
        #     for(j in 1:ngrids){
        #       ifelse(addxint,
        #              Xsgrid[,,j] <- as.matrix(gX6(X=data.frame(xgrid.subset[[j]]),info=info,orth=xorth)),
        #              Xsgrid[,,j] <- cbind(1,as.matrix(gX6(X=data.frame(xgrid.subset[[j]]),info=info,orth=xorth))))
        #     }
        #   }
        #
        # if("factor" %in% X.type){

        # ==============================================================
        #   Xsgrid <- gX6(X=data.frame(xgrid),info=info,orth=xorth)
        #   if(!is.matrix(Xsgrid)) Xsgrid <- as.matrix(Xsgrid)
        #   if(addxint){
        #     Xsgrid <- cbind(1,Xsgrid)
        #   }
        # }
        # ==============================================================

        Xsgrid <- list()
        Xsgrid.tmp <- list()

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
          Xsgrid.tmp[[i]] <- foreach(j = 1:ngrids, .packages=c("orthogonalsplinebasis", "splines","gtreg")) %dopar% {
            if(addxint) Xsgrid.now <- cbind(1,as.matrix(gX6(X=as.data.frame(xgrid.subset[[j]]),info=info,orth=xorth)))#,
            if(!addxint) Xsgrid.now <- as.matrix(gX6(X=as.data.frame(xgrid.subset[[j]]),info=info,orth=xorth))
            return(Xsgrid.now)
          }

          Xsgrid[[i]]  <- array(0,c(nxgrid,nXs,ngrids))
          for(j in 1:ngrids){
            Xsgrid[[i]][,,j]  <- Xsgrid.tmp[[i]][[j]]
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

    # if(nxgrid==length(ygrid)){
    #   if(!gridx.mode){Xsgrid <- Xs}
    #   if(!gridy.mode){SYgrid <- YS; sYgrid <- Ys}
    #   TZgrid    <- list()
    #   tZgrid    <- list()
    #   TZgrid.pl <- list()
    #   tZgrid.pl <- list()
    #   for(i in 1:length(Xsgrid)){
    #      TZgrid[[i]]    <- array(0,c(nxgrid,nXs*nYS,ngrids))
    #      tZgrid[[i]]    <- array(0,c(nxgrid,nXs*nYS,ngrids))
    #      TZgrid.pl[[i]] <- array(0,c(nxgrid,nXs*nYS,ngrids))
    #      tZgrid.pl[[i]] <- array(0,c(nxgrid,nXs*nYS,ngrids))
    #      for(j in 1:ngrids){
    #         TZgrid[[i]][,,j] <- TZ1(X=Xsgrid[[i]][,,j],Y=SYgrid)
    #         tZgrid[[i]][,,j] <- TZ1(X=Xsgrid[[i]][,,j],Y=sYgrid)
    #         if(is.null(mask)){ TZgrid.pl <- NULL; tZgrid.pl <- NULL }
    #         if(!is.null(mask)){
    #           TZgrid.pl[[i]][,,j] <- TZgrid[[i]][,,j][,which(mask==1)]
    #           tZgrid.pl[[i]][,,j] <- tZgrid[[i]][,,j][,which(mask==1)]
    #           }
    #       }
    #    }
    # }
    # if(length(xgrid)==length(ygrid)){
    #   if(!gridx.mode){Xsgrid <- Xs}
    #   if(!gridy.mode){SYgrid <- YS; sYgrid <- Ys}
    #   TZgrid <- TZ1(X=Xsgrid,Y=SYgrid)
    #   tZgrid <- TZ1(X=Xsgrid,Y=sYgrid)
    #   if(is.null(mask)){ TZgrid.pl <- NULL; tZgrid.pl <- NULL }
    #   if(!is.null(mask)){
    #     TZgrid.pl <- TZgrid[,which(mask==1)]
    #     tZgrid.pl <- tZgrid[,which(mask==1)]
    #   }
    # }
    # if(length(xgrid)!=length(ygrid)){
    #   if(!gridx.mode){Xsgrid <- Xs}
    #   if(!gridy.mode){SYgrid <- YS; sYgrid <- Ys}
    #   TZgrid <- array(0,c(nrow(SYgrid),ncol(SYgrid)*ncol(Xsgrid),nrow(Xsgrid)))
    #   tZgrid <- array(0,c(nrow(SYgrid),ncol(SYgrid)*ncol(Xsgrid),nrow(Xsgrid)))
    #   if(!is.null(mask)){
    #     TZgrid.pl <- list()
    #     tZgrid.pl <- list()
    #   }
    #   for(i in 1:nrow(Xsgrid)){
    #     TZgrid[,,i] <- TZ2(X=Xsgrid[i,],Y=SYgrid)
    #     tZgrid[,,i] <- TZ2(X=Xsgrid[i,],Y=sYgrid)
    #     if(!is.null(mask)){
    #       TZgrid.pl[[i]] <- TZgrid[,which(mask==1),i]
    #       tZgrid.pl[[i]] <- tZgrid[,which(mask==1),i]
    #     }
    #   }
    #   if(is.null(mask)){ TZgrid.pl <- NULL; tZgrid.pl <- NULL }
    # }

    #  if(nxgrid!=length(ygrid)){
    if(!gridx.mode){Xsgrid <- Xs}
    if(!gridy.mode){SYgrid <- YS; sYgrid <- Ys}
    TZgrid <- list()
    tZgrid <- list()
    TZgrid.pl <- list()
    tZgrid.pl <- list()
    if(is.list(Xsgrid)){
      for(i in 1:length(Xsgrid)){
        TZgrid[[i]] <- array(0,c(nygrid,nXs*nYS,nrow(Xsgrid[[i]]),ngrids))
        tZgrid[[i]] <- array(0,c(nygrid,nXs*nYS,nrow(Xsgrid[[i]]),ngrids))
        TZgrid.pl[[i]] <- array(0,c(nygrid,nXs*nYS,nrow(Xsgrid[[i]]),ngrids))
        tZgrid.pl[[i]] <- array(0,c(nygrid,nXs*nYS,nrow(Xsgrid[[i]]),ngrids))
        #      TZgrid <- array(0,c(nrow(SYgrid),ncol(SYgrid)*ncol(Xsgrid),nrow(Xsgrid)))
        #      tZgrid <- array(0,c(nrow(SYgrid),ncol(SYgrid)*ncol(Xsgrid),nrow(Xsgrid)))
        if(!is.null(mask)){
          TZgrid.pl <- list()
          tZgrid.pl <- list()
        }
        for(j in 1:ngrids){
          TZgrid[[i]][,,,j] <- TZ1(X=matrix(Xsgrid[[i]][1,,j],nr=nygrid,nc=ncol(Xsgrid[[i]]),byrow=T),Y=SYgrid)
          tZgrid[[i]][,,,j] <- TZ1(X=matrix(Xsgrid[[i]][1,,j],nr=nygrid,nc=ncol(Xsgrid[[i]]),byrow=T),Y=sYgrid)
          if(is.null(mask)){ TZgrid.pl <- NULL; tZgrid.pl <- NULL }
          if(!is.null(mask)){
            TZgrid.pl[[i]][,,,j] <- TZgrid[[i]][,,j][,which(mask==1)]
            tZgrid.pl[[i]][,,,j] <- tZgrid[[i]][,,j][,which(mask==1)]
          }
        }
      }
      if(is.matrix(Xsgrid)){
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
    #}

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
    dedy         <- diag(beta%*%t(Ys))
    e            <- diag(beta%*%t(YS))
    ans$beta     <- beta
    ans$dedy     <- dedy
    ans$e        <- e
    if(grid.mode){
      ifelse(is.null(de0dygrid), de0dygrid <- 1, de0dygrid <- de0dygrid)
      if(!gridx.mode){Xsgrid <- Xs}
      # betagrid     <- Xsgrid%*%bmat
      # egrid        <- betagrid%*%t(SYgrid)
      # dedygrid     <- betagrid%*%t(sYgrid)
      # phiegrid     <- data_norm(egrid)
      # fgrid        <- phiegrid*dedygrid*de0dygrid
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
          for(k in 1:nxgrid){
            egrid[[i]][k,,j]    <- betagrid[[i]][,,j][k,]%*%t(SYgrid)
            dedygrid[[i]][k,,j] <- betagrid[[i]][,,j][k,]%*%t(sYgrid)
            phiegrid[[i]][k,,j] <- data_norm(egrid[[i]][k,,j])
            fgrid[[i]][k,,j]    <- phiegrid[[i]][k,,j]*dedygrid[[i]][k,,j]*de0dygrid
          }
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
  #if(ydf>0){ ans$yknots <- yknots }
  ans$iyknots  <- iyknots
  ans$btarg    <- btarg

  return(ans)
}
