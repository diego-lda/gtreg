#' GTR Normalise
#'
#' @description This function performs a transformation to normality, or at least to the real line.
#'
#' @param y
#' @param nyg.max
#' @param plotting
#' @param eps
#' @param ygrid
#' @param log.t
#' @param doprimal
#' @param gam0
#'
#' @return
#' @export
#'
#' @examples
gtr_normalise <- function(y,nyg.max=0,plotting=F,eps=1e-5,ygrid=NULL,log.t=F,doprimal=F,gam0=0.1){

  ydf0     <- NULL
  b0 <- NULL
  e0mode <- F
  y.orig   <- y
  x0 <- matrix(1,nr=length(y))
  eta0g <- NULL
  wins <- F


  if(!log.t){

    bic0.vec <- NULL
    e0       <- NULL
    eta0     <- NULL

    exit     <- F
    kdex     <- 0

    while(length(bic0.vec)==0 && !exit){

      kdex <- kdex + 1

      # Nonlinear spec.
      for(yorder in 2:3){

        for(ydf in 2:3){

          mod0 <- gtr_al(y=y,x=x0,gam=gam0,ydf=ydf,yorder=yorder,maxit=500,doprimal=doprimal,nyg.max=nyg.max)

          if(length(mod0$BIC)>0){

            print(paste("Int = ", pnorm(mod0$res$e[which.max(y)[1]])-pnorm(mod0$res$e[which.min(y)[1]])))
            print(length(bic0.vec))

            if((length(bic0.vec)==0 || mod0$BIC<min(bic0.vec)) && (pnorm(mod0$res$e[which.max(y)[1]])-pnorm(mod0$res$e[which.min(y)[1]]) > .99#{
                                                                   || (pnorm(mod0$res$e[which.max(y)[1]])-pnorm(mod0$res$e[which.min(y)[1]]) > .95 && wins))
            ){

              print(c(yorder,ydf))
              e0 <- mod0$res$e
              eta0 <-mod0$res$eta
              b0 <- mod0$res$bmat
              ydf0 <- ydf
              e0mode <- T

              if(length(ygrid)>0){
                datmat  <- data_prep(y=as.numeric(y),x=x0,Xs=x0,ygrid=ygrid,iyknots=NULL,ydf=ydf0,
                                    addxint=TRUE,yorder=yorder,nxgrid=0,nygrid=length(ygrid),
                                    yorth=FALSE,xorth=FALSE,Ysing=Ysing,e0mode=F)
                eta0g <- datmat$sYgrid%*%b0
              }

              exit <- T
              if(plotting){
                qqnorm(mod0$res$e,main=paste("yorder=",yorder,"ydf=",ydf))
                plot(y,mod0$res$e,ylab="e0",main=paste("yorder=",yorder,"ydf=",ydf))
              }

              bic0.vec <- c(bic0.vec,mod0$BIC)
            }

          }


        }

      }

      if(!exit && kdex == 1 && min(y.orig)>0){
        print("============================")
        print("Applying log transformation.")
        print("============================")
        if(min(y)>0){
          y <- log(y)
        }
        wins <- T
      }

      if(!exit && (min(y.orig)<0 || kdex == 2)){
        exit <- T
        if(length(e0)==0 || min(e0)> -2 || max(e0)<2 || max(abs(e0)) > qnorm(1-eps) || (min(y.orig)>0 && sd(e0)>1.01)){

          if(kdex==1){
            e0 <- y
            eta0 <- 1
            if(length(ygrid)>0){
              eta0g<-1
            }
            print("No transformation applied!")
          };
          if(kdex==2){
            e0 <- y
            eta0 <- 1/exp(y)
            if(length(ygrid)>0){
              eta0g<-1/exp(ygrid)
            }
            e0mode <- T
          };
        }
      }
    }
  }

  if(log.t){
    e0 <- log(y)
    eta0 <- 1/y
    if(length(ygrid)>0){
      eta0g<-1/ygrid
    }
    e0mode <- T
  }

  stopifnot(length(e0)>0)

  return(list(e0=e0,eta0=eta0,eta0g=eta0g,b0=b0,ydf0=ydf0,e0mode=e0mode))
}
