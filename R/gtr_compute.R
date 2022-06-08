#' GTR Compute Function
#'
#' @description This function does the inner computation for GTR. Used to be called "inner gtr".
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
#' @param y_knots
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
gtr_compute <- function(TYX,tYX,y,x,
                        nyg=0,ng.qgm,weights=1,zeros=NULL,
                        reltol=1e-3,feastol=1e-3,abstol=1e-3,pen=NULL,
                        cvg.mono=NULL,res.sol=NULL,fac=1,fac.now2=NULL,
                        Xs=NULL,sYgrid=NULL,
                        xgrid.qgm=NULL,Xs.qgm=NULL,sYgrid.qgm=NULL,
                        info,yorder,ydf,y_knots=NULL,Ysing,maxit,nXs,nYS,
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
  exit        <- F

  while(length(cvg)==0 && tol.res*fac.now < 10 && length(res0)>0 && !exit){ # until cvg converges (to 1)

    cat(paste("fac =",fac.now,"\n"))
    res0 <- NULL
    if(nyg != 0){ pen <- T }
    if(length(pen)>0 && doprimal == F){
      cat("Switching to primal\n")
      doprimal <- T
    }

    try(res0 <- gtr_solve(TYX=TYX,tYX=tYX,algor=algor,maxit=maxit,doprimal=doprimal,nXs=nXs,nYS=nYS,Xs=Xs,threshold=threshold,
                           pen=pen,gam=gam,weights=weights,zeros=zeros,lam.vec=lam.vec,sYgrid=sYgrid,cval=cval,reltol=tol.res*fac.now,feastol=feastol,abstol=abstol,bounded=bounded,Cbound=Cbound,beta2=F))

    if(length(res0)<=1 || is.na(res0$llf) || (!(dim(x)[2] == 1 && all(x[,1] == 1)) && length(unique(round(res0$eta,digits=3)))==1)){
      # Try a different algorithm
      cat("Trying a different algorithm\n")
      ifelse(algor=="ECOS", algor<-"SCS", algor<-"ECOS"); ifelse(algor=="SCS", tol.res <- 1.e-01, tol.res <- 1e-4);
      try(res0 <- gtr_solve(TYX=TYX,tYX=tYX,algor=algor,maxit=maxit,doprimal=doprimal,nXs=nXs,nYS=nYS,Xs=Xs,threshold=threshold,
                             pen=pen,gam=gam,weights=weights,zeros=zeros,lam.vec=lam.vec,sYgrid=sYgrid,cval=cval,reltol=tol.res*fac.now,feastol=feastol,abstol=abstol,bounded=bounded,Cbound=Cbound,beta2=F))
      if(eta.check && length(pen)==0
         && (is.na(res0$llf) || length(res0)<=1 || ((!(dim(x)[2] == 1 && all(x[,1] == 1)) && length(unique(round(res0$eta,digits=3)))==1)))){
        exit <- T
        cat("No algorithm worked\n")
      }

    }


    # Check if high eta is not a computational issue
    # ----------------------------------------------
    if(length(res0)>1 && eta.check && !exit && gam == 0){
      #if(mean((res0$eta^2-mean(res0$eta^2))^2)>mean(res0$eta^2)){ exit <- T }
      gtest<-grubbs.test(log(res0$eta))
      print(paste("grubbs = ",gtest$p.value))
      if(gtest$p.value<.05){
        if(length(pen)==0){
          #print(sd(res0$eta))
          if(!is.na(res0$llf)){
            cat("dgap check in progress\n")
            doprimal.now <- !doprimal
            res.dual <- NULL
            try(res.dual <- gtr_solve(TYX=TYX,tYX=tYX,algor=algor,maxit=maxit,doprimal=doprimal.now,nXs=nXs,nYS=nYS,Xs=Xs,threshold=threshold,
                                       pen=pen,gam=gam,weights=weights,zeros=zeros,lam.vec=lam.vec,sYgrid=sYgrid,cval=cval,reltol=tol.res*fac.now,feastol=feastol,abstol=abstol,bounded=bounded,Cbound=Cbound,beta2=F))

            if(length(res.dual)>1 && !is.na(res.dual$llf)){
              print(abs(res.dual$llf - res0$llf))
              if(abs(res.dual$llf - res0$llf) > 1e-1 || (sd(res0$eta)>1 && abs(res.dual$llf - res0$llf) > 1e-6)){
                print(paste("dgap now =",abs(res.dual$llf - res0$llf)))
                exit<-T
              }#else{
              #if(res.dual$llf > res0$llf){
              #   res0 <- res.dual
              #   }
              #}
            }
            if(length(res.dual)<=1){
              exit <- T
              res0<-NULL
              cat("dgap check not available\n")
            }
          }
        }else{
          # Try a different algorithm
          res.alg2 <- NULL
          ifelse(algor=="ECOS", algor<-"SCS", algor<-"ECOS"); ifelse(algor=="SCS", tol.res <- 1.e-01, tol.res <- 1e-4);
          try(res.alg2 <- gtr_solve(TYX=TYX,tYX=tYX,algor=algor,maxit=maxit,doprimal=doprimal,nXs=nXs,nYS=nYS,Xs=Xs,threshold=threshold,
                                     pen=pen,gam=gam,weights=weights,zeros=zeros,lam.vec=lam.vec,sYgrid=sYgrid,cval=cval,reltol=tol.res*fac.now,feastol=feastol,abstol=abstol,bounded=bounded,Cbound=Cbound,beta2=F))
          if(length(res.alg2)>1){
            #if(length(res.alg2$eta)>0 && mean((res.alg2$eta^2-mean(res.alg2$eta^2))^2)<mean(res.alg2$eta^2)){
            if(!is.na(res.alg2$llf)){
              print(abs(res.alg2$llf - res0$llf))
              if(abs(res.alg2$llf - res0$llf) < 1e-1){
                if(res.alg2$llf > res0$llf){
                  cat("switching algo.\n")
                  res0 <- res.alg2
                }
              }else{
                exit<-T
                res0<-NULL
              }
              #}
            }
          }else{
            exit<-T
            res0<-NULL
          }
        }
      }
    }

    #print(exit)
    if(!exit){

      if(length(res0)>1 && ((dim(x)[2] == 1 && all(x[,1] == 1)) || length(unique(round(res0$eta,digits=3)))!=1)){
        #print(sd(res0$eta))
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
                cat("Calling qgm_check\n")
                t_start <- Sys.time()
                datmat.mono <- qgm_check(y=y,x=x,xgrid.qgm=xgrid.qgm,res=res0,Ysing=Ysing,iyknots=iyknots,ugrid=ugrid,ng.qgm=ng.qgm,nxgrid=ng.qgm,nygrid=ng.qgm,nXs=nXs,nYS=nYS,ydf=ydf,yorder=yorder,info=info,easy=F,zeros=zeros,e0mode=e0mode)
                t_end <- Sys.time()
                t_check <- t_end - t_start
                cat(paste("It took =",round(t_check,digits=3),"\n"))
                mdedy       <- min(datmat.mono$dedygrid)
              }
              rm(datmat.mono)
              cat(paste("min dedy now is =",round(mdedy,digits=6),"\n"))
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
              cat(paste("min dedy now is =",round(mdedy,digits=6),"\n"))
            }
          }

          if(nXs == 1){

            b2min <- sYgrid%*%matrix(M%*%res0$bmat,nr=nYS,nc=1)
            mdedy <- min(b2min)

            cat(paste("min dedy now is =",round(mdedy,digits=6),"\n"))
          }

          #      if(mdedy>.Machine$double.eps){

          if(((length(cvg.mono)==0 || (length(cvg.mono)>0 && mdedy>.Machine$double.eps) ) ) ){

            if(length(cvg.mono)==0){

              cvg      <- 1
              nyg.now  <- nyg
              cat(paste("min dedy now =",round(mdedy,digits=6),"\n"))

            }

            if(mdedy>.Machine$double.eps){

              cvg.mono <- 1
              res.sol  <- res0

              fac.now2 <- fac.now
              cat(paste("fac.now =",fac.now2,"\n"))
              dedy.min <- mdedy
              nyg.now  <- nyg

            }

            #        }

          }

        }

      }

      ifelse(algor=="ECOS" || nyg == 0, fac.now <- 11/tol.res, fac.now <- 2*fac*fac.now)

    }

  }

  ans <- list(res.sol=res.sol,res0=res0,cvg=cvg,nyg.now=nyg.now,fac.now2=fac.now2,dedy.min=dedy.min,algor=algor,tol.res=tol.res)

  return(ans)
}
