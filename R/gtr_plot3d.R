#' GTR Plot in 3D
#'
#' @description This function is only for plotting simple gtr output with a scalar X for now.
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
#' @return A plot of the gtr output with a scalar X.
#' @export
#'
#' @examples
gtr_plot3d <- function(y,x,ugrid=seq(.1,.9,by=.1),res,bmat,type = "CDF",ydf=NULL,yorder=NULL,info = NULL, Ysing = FALSE, ng.plot = 30,
                       ygrid = NULL, xgrid = NULL, dx = 0, dy = 0, FOV = 1, alpha = 1, ncolors = 500, coltype = 1, xlab = "X",
                       ylab = "Y", zlab = "", theta = 70, phi = 30, d = 1, expand = 1, r = sqrt(3), liness = F, shade=NA,
                       blackgrid = NULL, projrq = F, BW = NULL, cex = 1, ticktype = "simple", sizept = 3, colpt = "blue", colines = "red",
                       lwd.liness=1,lwd.projrq=1,col.liness=1,col.projrq=1,col.surface="lightblue",zlim=NULL,main=NULL,delta=1){

  ifelse(length(xgrid)==0, gx <- seq(min(x)-dx, max(x)+dx, len = ng.plot), gx <- xgrid)
  ifelse(length(ygrid)==0, gy <- seq(min(y)-dy, max(y)+dy, len = ng.plot), gy <- ygrid)

  datmat    <- data_prep(y=as.numeric(y),x=x,xgrid=gx,ygrid=gy,bmat=bmat,info=info,y_knots=NULL,ydf=ydf,
                        addxint=TRUE,yorder=yorder,yorth=FALSE,xorth=FALSE,Ysing=Ysing,e0mode=e0mode,delta=delta)
  z   <- datmat$egrid
  range(pnorm(z))

  # CDF
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

  # PDF
  if(type=="PDF"){
    eta <- res.al$eta

    if(length(zlim)==0){ zlim <- range(datmat$fgrid, na.rm = TRUE) }
    res.p  <- persp(gx,gy,datmat$fgrid,col=col.surface,xlab=xlab,ylab=ylab,zlab=zlab,theta=theta,phi=phi,d=d,r = r,expand=expand,ticktype=ticktype,border=NA,shade=.5,zlim=zlim,main=main)
    for(i in round(seq(1,ng.plot,len=10))){
      lines(trans3d(x=rep(gx[i],ng.plot),y=gy,z=datmat$fgrid[i,],pmat=res.p), col = 1,lwd=min(2,.5+1.5*i/ng.plot))
    }
  }

}
