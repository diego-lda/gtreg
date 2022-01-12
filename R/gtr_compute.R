#' Inner GTR Computation Function
#'
#' @description This function does the inner computation for GTR.
#'
#' @param TYX
#' @param tYX
#' @param y
#' @param x
#' @param nyg
#' @param ng.qgm
#' @param weights
#' @param zeros
#' @param reltol
#' @param feastol
#' @param abstol
#' @param pen
#' @param cvg.mono
#' @param res.sol
#' @param fac
#' @param fac.now2
#' @param Xs
#' @param sYgrid
#' @param xgrid.qgm
#' @param Xs.qgm
#' @param sYgrid.qgm
#' @param info
#' @param yorder
#' @param ydf
#' @param iyknots
#' @param Ysing
#' @param maxit
#' @param nXs
#' @param nYS
#' @param lam.vec
#' @param gam
#' @param ugrid
#' @param dedy.min
#' @param doprimal
#' @param tol.res
#' @param bounded
#' @param Cbound
#' @param cval
#' @param algor
#' @param easy
#' @param threshold
#' @param e0mode
#'
#' @return
#' @export
#'
#' @examples
inner.gtr.c <- function(TYX,tYX,y,x,
                        nyg=0,ng.qgm,weights=1,zeros=NULL,
                        reltol=1e-3,feastol=1e-3,abstol=1e-3,pen=NULL,
                        cvg.mono=NULL,res.sol=NULL,fac=1,fac.now2=NULL,
                        Xs=NULL,sYgrid=NULL,
                        xgrid.qgm=NULL,Xs.qgm=NULL,sYgrid.qgm=NULL,
                        info,yorder,ydf,iyknots=NULL,Ysing,maxit,nXs,nYS,
                        lam.vec=NULL,gam=0,ugrid=NULL,dedy.min=NULL,doprimal=F,
                        tol.res,bounded=F,Cbound=Inf,cval=1e-1,algor="ECOS",easy=T, threshold=1e-5,e0mode=F){

  M <- diag(1,nXs*nYS,nXs*nYS)
  nzeros <- length(zeros)
  if(nzeros>0){
    M            <- matrix(0,nXs*nYS,nXs*nYS-nzeros)
    M[-zeros,]   <- diag(1,nXs*nYS-nzeros,nXs*nYS-nzeros)
  }

  cvg         <- NULL
  fac.now     <- fac
  res0        <- 1
  nyg.now     <- 0

  while(length(cvg)==0 && tol.res*fac.now < 10 && length(res0)>0){ # until cvg converges (to 1)

    print(paste("fac =",fac.now))
    res0 <- NULL
    if(nyg != 0){ pen <- T }
    if(length(pen)>0 && doprimal == F){
      print("Switching to primal")
      doprimal <- T
    }

    try(res0 <- gtr_solve(TYX=TYX,tYX=tYX,algor=algor,maxit=maxit,doprimal=doprimal,nXs=nXs,nYS=nYS,Xs=Xs,threshold=threshold,
                           pen=pen,gam=gam,weights=weights,zeros=zeros,lam.vec=lam.vec,sYgrid=sYgrid,cval=cval,reltol=tol.res*fac.now,feastol=feastol,abstol=abstol,bounded=bounded,Cbound=Cbound,beta2=F))

    if(length(res0)<=1 || max(res0$eta) > 1e6 || ((!(dim(x)[2] == 1 && all(x[,1] == 1)) && length(unique(round(res0$eta,digits=3)))==1))){# || res0$result$status=="optimal_inaccurate"){ # means it doesn't work
      # Try a different algorithm
      ifelse(algor=="ECOS", algor<-"SCS", algor<-"ECOS"); ifelse(algor=="SCS", tol.res <- 1.e-01, tol.res <- 1e-4);
      try(res0 <- gtr_solve(TYX=TYX,tYX=tYX,algor=algor,maxit=maxit,doprimal=doprimal,nXs=nXs,nYS=nYS,Xs=Xs,threshold=threshold,
                             pen=pen,gam=gam,weights=weights,zeros=zeros,lam.vec=lam.vec,sYgrid=sYgrid,cval=cval,reltol=tol.res*fac.now,feastol=feastol,abstol=abstol,bounded=bounded,Cbound=Cbound,beta2=F))
    }

    if(length(res0)>1 && ((dim(x)[2] == 1 && all(x[,1] == 1)) || length(unique(round(res0$eta,digits=3)))!=1)){

      if(res0$result$num_iters==maxit){ fac <- fac*2 } # not converged
      # Check integrates to 1
      if(nYS>2 && nXs >1){

        if(!is.list(Xs)){
          b2min <- beta_check(bmat=matrix(M%*%res0$bmat,nr=nXs,nc=nYS),Xs=Xs,nXs=nXs,nYS=nYS)$b2min
        }
        if(is.list(Xs)){
          b2min <- NULL
          for(i in 1:length(Xs)){
            for(j in 1:dim(Xs[[i]])[3]){
              # b2min is the smallest value of beta2(X). It should be strictly positive.
              b2min <- c(b2min,beta_check(bmat=matrix(M%*%res0$bmat,nr=nXs,nc=nYS),Xs=Xs[[i]][,,j],nXs=nXs,nYS=nYS)$b2min)
            }
          }
          b2min <- min(b2min)
        }

        if(b2min<=.Machine$double.eps){
          try( res0 <- gtr_solve(TYX=TYX,tYX=tYX,algor=algor,lam.vec=lam.vec,doprimal=T,nXs=nXs,nYS=nYS,Xs=Xs,threshold=threshold,
                                  pen=pen,gam=gam,weights=weights,zeros=zeros,sYgrid=sYgrid,cval=cval,reltol=tol.res*fac.now,feastol=feastol,abstol=abstol,maxit=maxit,bounded=bounded,Cbound=Cbound,beta2=T) )
        }
        if(length(res0)==1
           || length(res0$result$status)==0
           || res0$result$status=="optimal_inaccurate"
           || length(which(round(res0$eta,digits=3)==cval))>0
           || (!(dim(x)[2] == 1 && all(x[,1] == 1)) && length(unique(round(res0$eta,digits=3)))==1)){

          # Try a different algorithm
          ifelse(algor=="ECOS",algor<-"SCS",algor<-"ECOS"); ifelse(algor=="SCS", tol.res <- 1.e-01, tol.res <- 1e-4);
          try( res0 <- gtr_solve(TYX=TYX,tYX=tYX,algor=algor,lam.vec=lam.vec,doprimal=T,nXs=nXs,nYS=nYS,Xs=Xs,threshold=threshold,
                                  pen=pen,gam=gam,weights=weights,zeros=zeros,sYgrid=sYgrid,cval=cval,reltol=tol.res*fac.now,feastol=feastol,abstol=abstol,maxit=maxit,bounded=bounded,Cbound=Cbound,beta2=T) )
        }

      }

      if(length(res0)>1 && ((dim(x)[2] == 1 && all(x[,1] == 1)) || length(unique(round(res0$eta,digits=3)))!=1)){

        # Check QGM
        if(nXs > 1){
          if(nYS>2){
            datmat.mono <- qgm_check(y=y,x=x,xgrid=xgrid.qgm,res=res0,Ysing=Ysing,iyknots=iyknots,ugrid=ugrid,ng.qgm=ng.qgm,nxgrid=ng.qgm,nygrid=ng.qgm,nXs=nXs,nYS=nYS,ydf=ydf,yorder=yorder,info=info,easy=T,zeros=zeros,e0mode=e0mode)
            mdedy       <- min(datmat.mono$dedygrid)
            if( mdedy>.Machine$double.eps && !easy ){
              print("Calling qgm_check")
              t_start <- Sys.time()
              datmat.mono <- qgm_check(y=y,x=x,xgrid.qgm=xgrid.qgm,res=res0,Ysing=Ysing,iyknots=iyknots,ugrid=ugrid,ng.qgm=ng.qgm,nxgrid=ng.qgm,nygrid=ng.qgm,nXs=nXs,nYS=nYS,ydf=ydf,yorder=yorder,info=info,easy=F,zeros=zeros,e0mode=e0mode)
              t_end <- Sys.time()
              t_check <- t_end - t_start
              print(paste("It took =",round(t_check,digits=3)))
              mdedy       <- min(datmat.mono$dedygrid)
            }
            rm(datmat.mono)
            print(paste("min dedy now is =",round(mdedy,digits=6)))
          }

          if(nYS==2){
            if(!is.list(Xs.qgm)){
              mdedy <- beta_check(bmat=matrix(M%*%res0$bmat,nr=nXs,nc=nYS),Xs=Xs.qgm,nXs=nXs,nYS=nYS)$b2min
            }
            if(is.list(Xs.qgm)){
              b2min <- NULL
              for(i in 1:length(Xs.qgm)){
                for(j in 1:dim(Xs.qgm[[i]])[3]){
                  # b2min is the smallest value of beta2(X). It should be strictly positive.
                  b2min <- c(b2min,beta_check(bmat=matrix(M%*%res0$bmat,nr=nXs,nc=nYS),Xs=Xs.qgm[[i]][,,j],nXs=nXs,nYS=nYS)$b2min)
                }
              }
              mdedy <- min(b2min)
            }
            print(paste("min dedy now is =",round(mdedy,digits=6)))
          }
        }

        if(nXs == 1){

          b2min <- sYgrid%*%matrix(M%*%res0$bmat,nr=nYS,nc=1)
          mdedy <- min(b2min)

          print(paste("min dedy now is =",round(mdedy,digits=6)))
        }

        if(mdedy>.Machine$double.eps){

          if(((length(cvg.mono)==0 || (length(cvg.mono)>0 && mdedy>.Machine$double.eps) ) ) ){

            if(length(cvg.mono)==0){

              cvg      <- 1
              nyg.now  <- nyg
              print(paste("min dedy now =",round(mdedy,digits=6)))

            }

            if(mdedy>.Machine$double.eps){

              cvg.mono <- 1
              res.sol  <- res0

              fac.now2 <- fac.now
              print(paste("fac.now =",fac.now2))
              dedy.min <- mdedy
              nyg.now  <- nyg

            }

          }

        }

      }

    }

    ifelse(algor=="ECOS" || nyg == 0, fac.now <- 11/tol.res, fac.now <- 2*fac*fac.now)

  }

  ans <- list(res.sol=res.sol,res0=res0,cvg=cvg,nyg.now=nyg.now,fac.now2=fac.now2,dedy.min=dedy.min,algor=algor,tol.res=tol.res)

  return(ans)

}
