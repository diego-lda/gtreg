#' XNSS2 Function
#'
#' @description This is a workhorse solver for the problem.
#'
#' @param y
#' @param x
#' @param info
#' @param ydf
#' @param iyknots
#' @param addxint
#' @param bstart
#' @param lam
#' @param lam.vec
#' @param gam
#' @param delta
#' @param pentype
#' @param yorth
#' @param xorth
#' @param Ysing
#' @param yorder
#' @param maxit
#' @param checkg
#' @param checkH
#' @param tval
#' @param algor
#' @param itprint
#' @param mask
#' @param secalc
#' @param HJreturn
#' @param bthresh
#' @param silent
#' @param fastreturn
#' @param reltol
#' @param diary
#' @param cval
#' @param theta
#' @param maxing
#' @param ygrid
#' @param btarg
#' @param egam
#' @param etarg
#' @param etgam
#' @param R
#'
#' @return
#' @export
#'
#' @examples
xnss2 <- function(y,x,info,ydf=4,iyknots=NULL,addxint=T,
                  bstart=NULL,lam=0,lam.vec=NULL,gam=0,delta=0,pentype=0,
                  yorth=FALSE,xorth=FALSE,Ysing=FALSE,yorder=4,maxit=300,
                  checkg=FALSE,checkH=FALSE,tval=0,
                  algor="optim",itprint=1,mask=NULL,secalc=F,HJreturn=F,
                  bthresh=.001,silent=F,fastreturn=F,
                  reltol=1.e-8,diary=NULL,cval=.01,theta=.9,maxing=FALSE,
                  ygrid=seq(min(y),max(y),len=101),btarg=NULL,egam=0,
                  etarg=NULL,etgam=0,
                  R=NULL){

  #cat("At entry to xnss2 Sys.time() is:",Sys.time(),"\n")
  systime1 <- Sys.time()
  if(is.null(iyknots)){iyknots <- quantile(y,probs=seq(0,1,length=ydf))}
  iyknots[[1]] <- plyr::round_any(iyknots[[1]],0.001,floor)
  iyknots[[length(iyknots)]] <- plyr::round_any(iyknots[[length(iyknots)]],0.001,ceiling)
  iyknots <- round(iyknots,3)
  names(iyknots) <- NULL
  yknots <- expand.knots(iyknots,order=yorder)
  ys <- SplineBasis(yknots,order=yorder)
  if(yorth)ys <- orthogonalize(ys)
  yS <- integrate(ys)
  Ys <- evaluate(object=ys,x=y)
  YS1 <- evaluate(object=yS,x=y)
  #if(ycenter){yshift <- apply(YS1,2,mean)}
  #if(!ycenter){yshift <-rep(0,length=ncol(YS1))}
  #cat("xnss2 processing at 60\n")

  if(!Ysing){
    YS <- YS1
    YS[,1] <- y
    YS <- cbind(1,YS)
    Ys <- Ys[,-1]
    Ys <- cbind(1,Ys)
    nYS <- ncol(YS)
    #cat("at line 69, nYS is ",nYS,"\n")
  }
  if(Ysing){
    YS <- YS1
    #YS[,1] <- y
    YS <- cbind(1,y,YS)
    #Ys <- Ys[,-1]
    Ys <- cbind(0,1,Ys)
    nYS <- ncol(YS)
    #cat("Ysing is TRUE\n")
  }


  #new in Feb 2018--define Ygrid and SYgrid sYgrid
  #print("wow starting SYgrid")
  #sorty <- sort(y,index.return=TRUE)
  #Ygrid <- sorty$x
  #ydex <- sorty$ix
  #cat("xnss2 processing at 72\n")
  nobs <- length(y)
  Ygrid <- ygrid  #tmp patch rhs
  SYfunc <- yS
  sYfunc <- ys
  SYy <- eval_basis(object=SYfunc,x=Ygrid)
  sYy <- eval_basis(object=sYfunc,x=Ygrid)

  #do the replacement
  if(!Ysing){
    SYgrid <- cbind(1,Ygrid,SYy[,-1])
    sYgrid <- cbind(0,1,sYy[,-1])
  }
  if(Ysing){
    SYgrid <- cbind(1,Ygrid,SYy)
    sYgrid <- cbind(0,1,sYy)
  }
  #Stu Feldman
  #TZ <- tz_form(X=Xs,Y=YS)
  #tZ.Ygrid <- tz_form(X=Xs,Y=sYy)

  #x processing begins.
  if(!is.data.frame(x)) x <- data.frame(x)  #cast x to data.frame always

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
    info.look <<- info
    x.look <<- x
    xnow <- gX6(X=x,info=info,orth=xorth)
    #xnow includes splined variables, no intercept, factor variables unexpanded
    xxnow <- data.frame(xnow[[1]],xnow)
    xxnow.look <<- xxnow
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
    Xsbd.look <<- Xs
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

  #so Xs and YS are relevant for TZ
  #and Xs and Ys for tZ
  #in old style, Ys needed a 'zero' col as in:
  if(!Ysing){Ys <-cbind(0,Ys)}  #so now it has it


  #Stu Feldman
  #cat("line 101 calls to tz_form in progress,\n")
  #cat("xnss2 processing at 119\n")
  Xs.look <<- Xs
  TZ <- tz_form(X=Xs,Y=YS)
  tZ <- tz_form(X=Xs,Y=Ys)
  #make synonym for YS,Ys to follow new notation/grid notation
  SY <- YS   #SY is numbers for S(Y)
  sY <- Ys   #sY is numbers for s(Y)
  #TZgrid <- tz_form(X=Xs,Y=SYgrid) #delayed from above because Xs not set
  #tZgrid <- tz_form(X=Xs,Y=sYgrid) #delayed from above because Xs not set
  #TZgrid.look <<- TZgrid
  #tZgrid.look <<- tZgrid
  #stop(180)
  #cat("line 178 calls to tz_form completed,\n")

  #while we are TZ'ing, make the 'giant constraint matrix' via tzb_form if
  #pentype==1
  #tzb_form <- function(Xs,sYgrid){
  if(pentype==1){A <- tzb_form(Xs,sYgrid)
  cat("giant A matrix dim's are:",dim(A),"\n")


  }
  #so this is like tZ for pseudo-ygrid data

  #b.compressed logic; use in place of !is.null(mask) because
  #ECOS and SCS use coef restriction to have the effect of mask and
  #thus they do not compress b
  b.compressed <- FALSE
  if(!is.null(mask)){
    b.compressed <- TRUE
    if(algor=="SCS") b.compressed <- FALSE
    if(algor=="ECOS") b.compressed <- FALSE
  }

  #prepare for Hessian contribution (etc.) computations
  #this needs to be fixed so wt can be excised rhs Feb 2018
  wt <- rep(1,length(y))
  sqwt <- sqrt(wt)
  TZ.wt <- TZ*sqwt
  tZ.wt <- tZ*sqwt
  #if you need to contract TZ.wt for Hessian use, do it here
  #if(!is.null(mask)){
  if(b.compressed){
    TZ.wt <- cm_contract(TZ.wt,buse=mask)
    tZ.wt <- cm_contract(tZ.wt,buse=mask)
  }
  #cat("xnss2 processing at 147\n")

  grx <- function(b){
    #ans <- grad(objfn1,x=b,method="complex")
    #ans <- numDeriv::grad(objfn1,x=b,method="complex")
    ans <- numDeriv::grad(objfn1,x=b,method="simple")
    return(ans)
  }

  gr2 <- function(b,Hmode=FALSE){
    if(b.compressed){bmat <- vexpand(b,buse=mask)}
    if(!b.compressed){bmat <- matrix(b,nrow=nXs,ncol=nYS)}

    #new style for xss5
    b <- matrix(bmat,nc=1)
    e <- TZ%*%b
    dedy <- tZ%*%b
    part1 <- -TZ*as.vector(e)
    part2 <- tZ/as.vector(dedy)
    gradmat <- part1+part2

    #previously contraction of the gradient for masked cases
    #took place here
    #now with modular penalties this will be delayed until
    #after penalty code

    gradmat.look <<- gradmat
    grad <- apply(gradmat,2,sum)
    totgrad <- grad #so this is what happens for pentype==0
    #sign changed at return for !maxing

    #put the penalty here for now
    if(pentype==2){
      llfvec <- log(data_norm(e)*dedy)
      hpen <- penfunc3(bmat,Xs,Ygrid,SYgrid,sYgrid,llfvec,
                       pickmat,gradmode=TRUE)
      llf <- sum(llfvec)
      ellf <- hpen$ellf
      dedyrange <- hpen$dedyrange
      cat("llf=",llf," ellf=",ellf," dedyrange=",dedyrange,"\n")
      ellfgrad <- hpen$ellfgrad
      #grad is not in 'bmat' form so make ellfgrad similar
      ellfgrad <- as.vector(ellfgrad)
      totgrad <- tval*grad+ellfgrad
      #totgrad <- (1-gam)*grad + gam*ellfgrad
    }
    if(pentype==3){
      llfvec <- log(data_norm(e)*dedy)
      hpen <- penfunc3(bmat,Xs,Ygrid,SYgrid,sYgrid,llfvec,
                       pickmat,gradmode=TRUE)
      llf <- sum(llfvec)
      cat("llf=",llf,"\n")
      ellf <- hpen$ellf
      ellfgrad <- hpen$ellfgrad
      #grad is not in 'bmat' form so make ellfgrad similar
      ellfgrad <- as.vector(ellfgrad)
      totgrad <- grad-delta*(llf-ellf)*(grad-ellfgrad)
      #totgrad <- (1-gam)*grad + gam*ellfgrad
    }
    if(pentype==4){
      #bmat or b in this call? rhs Feb 2018
      hpen <- penfunc4(bmat,lam,gradmode=TRUE)
      #call penfunc3 to see what is going on....
      llfvec <- log(data_norm(e)*dedy)
      hpen3 <- penfunc3(bmat,Xs,Ygrid,SYgrid,sYgrid,llfvec,
                        pickmat,gradmode=FALSE)

      gradterm <- hpen$gradterm
      #cat("reporting gam as ",gam,"\n")
      #print(gradterm[1:20])
      totgrad <- grad+gam*gradterm  #is this correct????
    }
    if(pentype==5){
      #doing pentype==3 piece first
      llfvec <- log(data_norm(e)*dedy)
      llf <- sum(llfvec)
      hpen3 <- penfunc3(bmat,Xs,Ygrid,SYgrid,sYgrid,llfvec,
                        pickmat,gradmode=TRUE)
      ellf <- hpen3$ellf
      ellfgrad <- hpen3$ellfgrad
      #grad is not in 'bmat' form so make ellfgrad similar
      ellfgrad <- as.vector(ellfgrad)
      totgrad1 <- grad-delta*(llf-ellf)*(grad-ellfgrad)
      #totgrad <- (1-gam)*grad + gam*ellfgrad
      hpen4 <- penfunc4(bmat,lam,gradmode=TRUE)
      penval4 <- hpen4$penval
      penval3 <- (llf-ellf)^2
      llf <- sum(llfvec)
      overfit <- llf-ellf
      cat("llf, ellf, bnorm penalty, overfit=",
          llf,ellf,penval4,overfit,"\n")
      gradterm <- hpen4$gradterm
      totgrad <- totgrad1+gam*gradterm  #is this correct????
    }

    if(pentype==6){
      llfvec <- log(data_norm(e)*dedy)
      hpen <- penfunc6(bmat,Xs,Ygrid,SYgrid,sYgrid,
                       gradmode=TRUE)
      llf <- sum(llfvec)
      penval <- hpen$penval
      dedyrange <- hpen$dedyrange
      cat("llf=",llf," penval=",penval,
          " dedyrange=",dedyrange,"\n")
      pengrad <- as.vector(hpen$pengrad)
      totgrad <- tval*grad+pengrad
      #totgrad <- (1-gam)*grad + gam*ellfgrad
    }

    if(b.compressed){
      gradmat <- cm_contract(gradmat,buse=mask)
      totgrad <- vcontract(as.vector(totgrad),buse=mask)
    }
    if(maxing){return(totgrad)}
    return(-totgrad)

  }


  Hess2 <- function(b){
    time1 <- Sys.time()
    if(b.compressed){bmat <- vexpand(b,buse=mask)}
    if(!b.compressed){bmat <- matrix(b,nrow=nXs,ncol=nYS)}
    #so bmat is something like 52 by 7

    #new style for xss5
    b <- matrix(bmat,nc=1)
    #e <- TZ%*%b   #not needed
    dedy <- tZ%*%b

    #components of the Hessian come in (weighted) contribution
    #form

    #Part1 <- TZ.wt
    Part2 <- tZ.wt/as.vector(dedy)
    Hess <- -(t(TZ.wt)%*%TZ.wt+t(Part2)%*%Part2)
    if(!maxing)Hess <- -Hess
    time2 <- Sys.time()
    #cat("time in Hess2:", time2-time1,"dim(Hess)=",dim(Hess),"\n")
    return(Hess)

  }

  Hess2.unmasked <- function(b){
    #full unmasked Hessian uses unmasked TZ, tZ
    #not suitable for weighted estimation (if that is ever to be done)
    #new mask logic here March 19 2018
    time1 <- Sys.time()
    if(b.compressed){bmat <- vexpand(b,buse=mask)}
    if(!b.compressed){bmat <- matrix(b,nrow=nXs,ncol=nYS)}

    #new style for xss5
    b <- matrix(bmat,nc=1)
    #e <- TZ%*%b   #not needed
    dedy <- tZ%*%b

    #components of the Hessian come in (weighted) contribution
    #form
    #here we use unweighted because that is unmasked....
    #assume all wts 1


    #Part1 <- TZ.wt
    Part2 <- tZ/as.vector(dedy)
    Hess <- -(t(TZ)%*%TZ+t(Part2)%*%Part2)
    if(!maxing)Hess <- -Hess
    time2 <- Sys.time()
    #cat("time in Hess2.unmasked:", time2-time1,"dim(Hess)=",dim(Hess),"\n")
    return(Hess)

  }

  objfn1 <- function(b){
    if(b.compressed){bmat <- vexpand(b,buse=mask)}
    if(!b.compressed){bmat <- matrix(b,nrow=nXs,ncol=nYS)}

    #new style for xss5
    b2 <- matrix(bmat,nc=1) #use b2 to preserve b for report use
    #below
    #cat("bdim at line 584 is:",bdim,"\n")
    if(!silent){cat("dim TZ at line 405 is:",dim(TZ),"\n")
      cat("dim b2 at line 405 is:",dim(b2),"\n")
    }
    e <- TZ%*%b2
    dedy <- tZ%*%b2
    llfvec <- log(data_norm(e)*dedy)
    llf <- sum(llfvec)
    if(is.nan(llf)){return(NaN)}

    if(!report){
      #cat("pentype in objfn1 is",pentype,"\n")
      if(pentype==2){
        hpen <- penfunc3(bmat,Xs,Ygrid,SYgrid,sYgrid,llfvec,pickmat)
        #penval <- hpen$penval
        #penllf <- llf   -(gam*penval)
        #penllf <- (1-gam)*llf+gam*hpen$ellf
        ellf <- hpen$ellf
        penllf <-  tval*llf+ellf
        objval <-  -penllf #assuming optimizer prog is minimizing
      }
      if(pentype==3){
        hpen <- penfunc3(bmat,Xs,Ygrid,SYgrid,sYgrid,llfvec,pickmat)
        #penval <- hpen$penval
        #penllf <- llf   -(gam*penval)
        #penllf <- (1-gam)*llf+gam*hpen$ellf
        penllf <- llf-(delta/2)*(llf-hpen$ellf)^2
        objval <-  -penllf #assuming optimizer prog is minimizing
      }
      if(pentype==4){
        #L1 penalization of b by conforming lam
        hpen <- penfunc4(b,lam,gradmode=F)
        penval <- hpen$penval
        penllf <- llf -gam*penval
        objval <- -penllf
      }
      #new pentype Aug 18 2018, both pentype==3 &pentype==4
      if(pentype==5){
        #get the lasso penalty
        hpen4 <- penfunc4(b,lam,gradmode=F)
        penval4 <- hpen4$penval
        #get the ELLF penalty
        hpen3 <- penfunc3(bmat,Xs,Ygrid,SYgrid,sYgrid,
                          llfvec,pickmat)
        penval3 <- (llf-hpen3$ellf)^2
        penllf <- llf-(delta/2)*penval3-gam*penval4
        objval <- -penllf
      }
      if(pentype==6){
        hpen <- penfunc6(bmat,Xs,Ygrid,SYgrid,sYgrid)
        penval <- hpen$penval
        penllf <-  tval*llf+penval  #not really a "penalty"
        objval <-  -penllf #assuming optimizer prog is minimizing
      }


      if(is.null(pentype)|pentype==0){objval <- -sum(llfvec)}
      #print(paste(maxing,objval))
      if(maxing)return(-objval)
      if(!maxing)return(objval)
    }
    #so we are reporting
    bic <- -2*llf+length(b)*log(length(y))
    beta <- Xs%*%bmat  #this is useful in interpretation
    grad <- gr2(b)
    #cat("think about grad for doECOS\n")
    if(b.compressed){grad <- vexpand(grad,buse=mask)}
    if(!b.compressed){grad <- matrix(grad,nrow=nXs,ncol=nYS)}

    report <- list(e=e,bmat=bmat,llfvec=llfvec,llf=llf,
                   dedy=dedy,YS=YS,ys=ys,beta=beta,
                   y=y,x=x,xm=xm,bic=bic,lam=lam,Ys=Ys,
                   info=info,
                   #penalty=pendel+pengam,
                   grad=grad)
    return(report)
  }

  objfn2 <- function(b){
    objfn2.count <<- objfn2.count+1
    value <- objfn1(b)
    #print(c(objfn2.count,round(value,4)))
    if(is.nan(value)){
      return(list(value= -Inf))
    }
    gradient <- gr2(b)
    #hessian <- jacobian(gr2,bstart,method="Richardson")
    hessian <- Hess2(b)
    ans <- list(value=value,gradient=gradient,
                hessian=hessian)
    return(ans)
  }





  if(is.null(bstart)){
    #since(yform), y lives in YS[,2]
    muy <- mean(YS[,2])
    sdy <- sd(YS[,2])
    bstart <- rep(0,len=nXs*nYS)
    bstart[1] <- -muy/sdy
    bstart[nXs+1] <- 1/sdy
    bstart.look <<- bstart #might be useful
  }

  if(is.null(btarg)){
    #since(yform), y lives in YS[,2]
    muy <- mean(YS[,2])
    sdy <- sd(YS[,2])
    btarg <- rep(0,len=nXs*nYS)
    btarg[1] <- -muy/sdy
    btarg[nXs+1] <- 1/sdy
  }

  if(b.compressed){bstart <- vcontract(bstart,buse=mask)}

  #new handling of lam
  #even newer--rhs 14 Mar 2017

  lam.orig <- lam   #who knows? leave breadcrumbs
  lamspecified <- is.matrix(lam)&is.null(lam.vec)
  #lamspecified <- is.matrix(lam)#&is.null(lam.vec)
  if(!lamspecified){
    if(!is.matrix(lam)){
      lam <- matrix(lam,nr=nXs,nc=nYS)
    }
    if(is.matrix(lam)){
      for(i in 1:length(lam)) lam[i] <- lam.vec[6] #gen
      lam[1,] <- lam.vec[5] #row1
      lam[,1] <- lam.vec[3] #col1
      lam[,2] <- lam.vec[4] #col1
      lam[1,1] <- lam.vec[1] #int
      lam[1,2] <- lam.vec[2] #int
    }
  }
  if(lamspecified){cat("lam is directly specified\n")
    cat("range(lam) is:",range(lam),"\n")
    lam <- abs(lam) #just in case, must be true
    #lam.vec <- rep(0,6)#so that it is not null, which confuses
    #logic below

  }
  #if lamspecified, we come down to here.
  lam.mat <- lam  #this is the unmasked form of lam
  if(b.compressed&!is.null(lam))lam <- vcontract(lam,mask)
  #if(!is.null(mask)&!is.null(lam))lam <- vcontract(lam,mask)
  #lam is now masked and is not altered until the end, lam.mat will
  #be returned as lam

  report <- FALSE
  Hmode <- FALSE
  maxing <- FALSE
  if((algor=="maxNR")|(algor=="trust")|(algor=="maxBFGSR")){maxing <- TRUE}
  if(checkg){
    #g0.look <<- grad(func=objfn1,x=bstart,method="Richardson")
    #g0.look <<- grad(func=objfn1,x=bstart,method="simple")
    g2.look <<- gr2(bstart)
    #stop(489)
    cat("in checkg\n")
    g1.look <<- grx(bstart)
    stop(127)
  }
  if(checkH){
    Hmode <- TRUE
    H2.look <<- Hess2(bstart)
    Hmode <- FALSE
    #H0.look <<- hessian(func=objfn1,x=bstart,method="Richardson")
    H1.look <<- jacobian(gr2,bstart,method="Richardson")
    stop(155)
  }

  maxing <- FALSE


  #summary of pentypes Feb 2018==========
  #pentype==3 get ELLF and dedy.grid; in CVXR dedy>0 imposed on the grid
  #pentype==4 L1 penalization 'pure'



  #============CVXR algorithm calls======================
  #wholesale rewrite/simplification in prep for CVXR
  #doCVXR <- FALSE
  if(algor=="SCS"|algor=="ECOS"){
    TZfull <- TZ
    tZfull <- tZ
    bdim <- ncol(TZ)
    lam <- lam.mat
    lam.look <<- lam

    #handle R here
    if(!is.null(R)){
      #rotation in use
      #the rotation matrix is K by Kp or nXs by nXsp
      #store TZ,tZ,Xs for restoration later
      TZ.R0 <- TZ
      tZ.R0 <- tZ
      Xs.R0 <- Xs
      Xs <- Xs%*%R
      #Xs[,1] <- 1  #intercept
      #intercept will be handled by construction of R
      TZ <- tz_form(X=Xs,Y=YS)
      tZ <- tz_form(X=Xs,Y=Ys)
      bdim <- ncol(TZ)
      Kp <- ncol(Xs)
      nXs <- Kp
      #this is not right for raw Y
      btarg <- matrix(0,nr=Kp,nc=ncol(YS))
      btarg[1,2] <- 1.
      btarg <- as.vector(btarg)
      lam.R0 <- lam
      lam <- matrix(0,nr=Kp,nc=ncol(YS))
      for(i in 1:length(lam)) lam[i] <- lam.vec[6] #gen
      lam[1,] <- lam.vec[5] #row1
      lam[,1] <- lam.vec[3] #col1
      lam[,2] <- lam.vec[4] #col1
      lam[1,1] <- lam.vec[1] #int
      lam[1,2] <- lam.vec[2] #int
    }



    if(!silent){
      cat("bdim at line 584 is:",bdim,"\n")
      cat("dim TZ at line 584 is:",dim(TZ),"\n")
    }
    Kgauss <- log(1/sqrt(2*pi))
    #M is the 'mask matrix' M%*%b=0, M determined by mask
    if(!is.null(mask)){
      nmask <- bdim-sum(mask)
      M <- matrix(0,nr=nmask,nc=bdim)
      rownow <- 0
      #Where is Stu Feldman?
      for(i in 1:bdim){
        if(as.integer(mask[i])==0){
          rownow <- rownow+1
          M[rownow,i] <- 1
        }}
      M.look <<- M
    }
    #check that the constraint grid looks the way we think
    #print(dim(t(sYgrid)))
    if(maxit>0){
      #b <- Variable(rows=nXs,cols=nYS)
      b <- Variable(bdim)

      #Beta <- Variable(nobs*nYS)
      LLF <- function(b){
        e <- TZ%*%b
        dedy <- tZ%*%b
        llfvec <- -.5*e^2+log(dedy)+Kgauss
        sum(llfvec)
      }
      E <- function(b){
        e <- TZ%*%b
        sum(e)/nobs
      }
      enow <- function(b){
        e <- TZ%*%b
        return(e)
      }
      V <- function(b){
        e <- TZ%*%b
        sumsq <- p_norm(e,2)
        sumsq/nobs -1
      }
      #       LLF.grid2 <- function(b){
      #                    e <- TZgrid%*%b
      #                    dedy <- tZgrid%*%b
      #                    llfvec <- -.5*e^2+log(dedy)+Kgauss
      #                    sum(llfvec)
      #                    }
      dedy <- function(b){
        dedy <- tZ%*%b
        dedy
      }
      #        LLF.grid <- function(b){
      #                     Beta <- Xs%*%reshape_expr(b,nXs,nYS)
      #                     e <- Beta%*%t(SYgrid)
      #                     dedy <- Beta%*%t(sYgrid)
      #                     llfvec <- -.5*e^2+log(dedy)+Kgauss
      #                     sum(llfvec)
      #                     }
      #        LLF.Bgrid <- function(Beta){
      #                     BX <- reshape_expr(Beta,nobs,nYS)
      #                     e <- BX%*%t(SYgrid)
      #                     dedy <- BX%*%t(sYgrid)
      #                     llfvec <- -.5*e^2+log(dedy)+Kgauss
      #                     sum(llfvec)
      #                     }
      #        LLF.Beta <- function(Beta){
      #                     #specifying the LLF as fctn of Beta(X)
      #                     BX <- reshape_expr(Beta,nobs,nYS)
      #                     e <- BX%*%t(SY)
      #                     dedy <- BX%*%t(sY)
      #                     llfvec <- -.5*e^2+log(dedy)+Kgauss
      #                     sum(llfvec)
      #                     }

      dedy.grid <- function(b){
        Beta <- Xs%*%reshape_expr(b,nXs,nYS)
        x <- Beta%*%t(sYgrid)  #x is dedyhat abbrev.
        a <- -40.
        return(min(x))

      }

      lamx <- as.vector(lam)
      #reg <- gam*sum(lamx*abs(b-btarg))
      reg <- gam*sum(lamx*abs(b))
      #elastic <- egam*(sum((b-btarg)^2))
      elastic <- egam*p_norm((b-btarg),2)
      #elastic <- egam*((b[1]-btarg[1])^2+(b[nXs+1]-btarg[nXs+1])^2)

      #adding etarg behavior here July 8 2019 rhs
      if(is.null(etarg)){etgam <- 0
      etargpen<- 0}
      else{
        etargpen <- etgam*p_norm(enow(b)-etarg,2)
      }

      Mb <- function(b){M%*%b}
      obj <- (LLF(b)-reg-elastic-etargpen)
      #obj <- (.998*LLF.Beta(Beta)+.002*LLF.Bgrid(Beta)-reg)
      constraint.condns<- dedy.grid(b)>= cval
      location <- E(b)==0.0
      scale <- V(b)==0
      #Wb <- (reshape_expr(Beta,nobs,nYS)-Xs%*%reshape_expr(b,nXs,nYS))==0
      if(is.null(mask)){
        constr <-list(constraint.condns)
        constr2 <- list(scale)
        #constr2 <- list(Wb)
        if(pentype==4)prob <- Problem(Maximize(obj))
        #if(pentype==4)prob <- Problem(Maximize(obj),constr2)
        if(pentype==3)prob <- Problem(Maximize(obj),constr)
      }
      if(!is.null(mask)){
        eqconditions<- Mb(b)==0
        eqcond <- list(eqconditions)
        bothcond <- list(constraint.condns,eqconditions)
        allcond <- list(eqconditions,location,scale)
        if(pentype==4)prob <- Problem(Maximize(obj),allcond)
        if(pentype==3)prob <- Problem(Maximize(obj),bothcond)
      }
      ecosolve.verbose <- as.integer(!silent)
      arglist <- ECOSolveR::ecos.control(maxit=as.integer(maxit),verbose=ecosolve.verbose)
      problem.data <- get_problem_data(prob,solver="ECOS")
      pdata.look <<- problem.data
      pd.scs <- get_problem_data(prob,solver="SCS")
      arglist2 <- list(max_iters=maxit,eps=reltol)
      #let's try a direct call to the solver
      #new sequence, I am novice in try() and will fake time.CVXR
      if(algor=="ECOS"){
        time.CVXR=3.14159
        ecos.out1 <- try(
          solve(prob,solver="ECOS",verbose=!silent,MAXIT=as.integer(maxit))
        )
        if(class(ecos.out1)=="try-error"){
          return(ecos.out1)
          #error-handling here
        }
        result1 <- ecos.out1
        result1.look <<-  result1
        if(!silent)cat("at line 697 ...result1.look saved\n")
        if(!silent)cat("time.CVXR is",time.CVXR,"\n")
        directcall <- FALSE
        if(directcall){
          time.CVXR <- system.time(
            ecos.out1 <-  ECOSolveR::ECOS_csolve(
              c=problem.data[["c"]],
              G=problem.data[["G"]],
              h=problem.data[["h"]],
              dims=problem.data[["dims"]],
              A=problem.data[["A"]],
              b=problem.data[["b"]],control=arglist
            )
          )
          result1 <- unpack_results(prob,"ECOS",ecos.out1)
        } #directcall end...
      }
      if(algor=="SCS"){time.CVXR <- system.time(
        scs.out1 <- scs::scs(
          A = pd.scs[['A']],
          b = pd.scs[['b']],
          obj = pd.scs[['c']],
          cone = pd.scs[['dims']],
          control=arglist2
        )
      )
      result1 <- unpack_results(prob, "SCS", scs.out1)
      }

      result1.look <<- result1
      #cat("line727.")
      if(result1$status=="solver_error"){
        #solver failed, try to return gracefully
        ans <- list(status="failed")
        return(ans)
      }
      bfinal <- result1$getValue(b)
      bfinal.look <<- bfinal

      #print(result1$getValue(obj))
      #print(result1$getValue(reg))


      if(!is.null(mask)){#bfinal <- vexpand(bfinal,buse=mask)
        TZ <- TZfull
        tZ <- tZfull
      }
    }

    #undo coef rotation here
    if(!is.null(R)){
      bfinal <- matrix(bfinal,nc=nYS,byrow=F)
      bfinal <- as.vector(R%*%bfinal)
      bfinalR.look <<- bfinal
      TZ.R0 -> TZ
      tZ.R0 -> tZ
      Xs.R0 -> Xs
      nXs <- ncol(Xs)  #seems to be problem at report stage
      #restore lam
      lam.R0 -> lam
    }

  }
  systime2 <- Sys.time()
  #cat("At line 737 time spent in xnss2 is:",systime2-systime1,"\n")


  #============optim (etc.) algorithm calls======================

  #===homebrew1 central path method experimental Aug 26 2018===
  if(algor=="homebrew"){

    cat("executing homebrew...\n")
    #get objfn for initial point
    if(!is.null(mask)){
      bdim <- sum(mask)
      cat("model dimension is:",bdim,"\n")
    }
    pentype <- 6
    report <- TRUE
    init <- objfn1(bstart)
    init.look <<- init
    cat("initial llf=",init$llf,"\n")
    #compute Hunmasked
    Hunmasked <- fastH1(list(TZ=TZ, tZ=tZ, dedy=init$dedy))
    cat("Hunmasked computed\n")

    #define wrappers for checking hpen6 grad and Hessian
    vpen6 <- function(bmat){
      hpen <- penfunc6(bmat=bmat,Xs=Xs,Ygrid=ygrid,SYgrid=SYgrid,
                       sYgrid=sYgrid,gradmode=FALSE,Hessmode=FALSE,
                       quiet=FALSE)
      return(hpen$penval)
    }
    gpen6 <- function(bmat){
      hpen <- penfunc6(bmat=bmat,Xs=Xs,Ygrid=ygrid,SYgrid=SYgrid,
                       sYgrid=sYgrid,gradmode=TRUE,Hessmode=FALSE,
                       quiet=FALSE)
      return(as.vector(hpen$pengrad))
    }

    checkpgrad <- FALSE
    if(checkpgrad){
      #let's check the reported gradient
      numgrad <- numDeriv::grad(vpen6,init$bmat,method="Richardson")
      agrad <- gpen6(init$bmat)
      cat("The range of numgrad-agrad is",range(numgrad-agrad),"\n")
    }


    checkpHess <- TRUE
    if(checkpHess){
      #let's check the reported Hessian
      cat("checking Hessian...\n")
      hpen6 <- penfunc6(bmat=init$bmat,Xs=Xs,Ygrid=ygrid,SYgrid=SYgrid,
                        sYgrid=sYgrid,gradmode=TRUE,Hessmode=TRUE,
                        quiet=FALSE)
      bx <- (init$bmat)
      time1 <- proc.time()
      #numhess <- numDeriv::jacobian(gpen6,bx,method="Richardson")
      #cat("numDeriv jacobian time=",proc.time()-time1,"\n")
      #numhess.look <<- numhess
      ahess.look <<- hpen6$Hess
      time1 <- proc.time()
      numhess2 <- pracma::jacobian(gpen6,bx)
      cat("pracma jacobian time=",proc.time()-time1,"\n")
      numhess <- numhess2
      numhess2.look <<- numhess2
      cat("range of numhess-ahess is:",range(numhess-hpen6$Hess),"\n")
      #call penfunc6 to get barrier Hessian, barrier grad
      #  hpen6 <- penfunc6(bmat=init$bmat,Xs=Xs,Ygrid=ygrid,SYgrid=SYgrid,
      #                    sYgrid=sYgrid,gradmode=TRUE,Hessmode=TRUE,
      #                    quiet=FALSE)
      #stop(731)
    }


    #the initial call with report==TRUE includes a gr2() call
    #Gb <- gr2(bstart)
    #call penfunc6 to get penalty value and deriv
    hpen6 <- penfunc6(bmat=init$bmat,Xs=Xs,Ygrid=ygrid,SYgrid=SYgrid,
                      sYgrid=sYgrid,gradmode=TRUE,Hessmode=FALSE,
                      quiet=FALSE)
    hpen6.look <<- hpen6
    Hunmasked.look <<- Hunmasked

    #compress the Hessians before computing directions, etc.
    maskv <- as.vector(mask)
    H1 <- Hunmasked[maskv,maskv]
    #Hb <- hpen6$Hess[maskv,maskv]
    Hb <- numhess[maskv,maskv]
    Htot <- (tval*H1+1*Hb)
    Hinv <- solve(Htot)
    Gb <- init$grad[maskv]  #the relevant components of the grad
    report <- TRUE
    bnow <- bstart
    b.compressed <- TRUE
    d <- Hinv%*%Gb  #d seems to have the correct sign this way..
    for(i in 1:60){
      #d <- Hinv%*%Gb  #d seems to have the correct sign this way..
      bnow <- bnow+.02*d
      #cat("i and dim(bnow) are:",i,dim(bnow),"\n")
      hold1 <- objfn1(bnow)
      bmat <- vexpand(bnow,mask)
      hpen6 <- penfunc6(bmat=bmat,Xs=Xs,Ygrid=ygrid,
                        SYgrid=SYgrid, sYgrid=sYgrid,gradmode=FALSE,
                        Hessmode=FALSE, quiet=FALSE)
      Gb <- hold1$grad[maskv]
      llf <- hold1$llf
      penval <- hpen6$penval
      objval <- tval*llf+penval
      cat("i is ",i,llf,penval,objval,"\n")
    }


    stop(710)
  }




  if(algor=="optim"){
    time.optim <- system.time(
      ans.optim <- optim(par=bstart,fn=objfn1,gr=gr2,method="BFGS",control=list(REPORT=itprint,trace=as.integer(!silent),maxit=maxit,reltol=reltol))
      #ans.optim <- ucminf(par=bstart,fn=objfn1,gr=gr2,control=list(trace=as.integer(!silent),maxeval=maxit))
    )
    bfinal <- ans.optim$par
    bfinal.look <<- bfinal
  }

  if(algor=="nlminb"){
    stepmax <- 10
    nlm.ans <- nlminb(objective=objfn1,gradient=gr2,start=bstart,
                      control=list(trace=0,iter.max=maxit,eval.max=3*maxit,step.max=stepmax))
    nlm.ans.look <<- nlm.ans
    bfinal <- nlm.ans$par
    ans.optim <- nlm.ans
  }

  if(algor=="maxNR"){
    maxing <- TRUE
    #ans.optim <- maxNR(fn=objfn1,grad=gr2,start=bstart,print.level=1)
    ans.optim <- maxNR(fn=objfn1,grad=gr2,hess=Hess2,start=bstart,
                       print.level=0,iterlim=maxit,
                       reltol=reltol,tol=reltol)
    maxNR.ans.look <<- ans.optim
    bfinal <- ans.optim$estimate
  }
  #============main algorithm call ends======================

  systime2 <- Sys.time()
  #cat("At line 885 time spent in xnss2 is:",systime2-systime1,"\n")

  if(maxit<=0){bfinal <- bstart}


  report <- TRUE
  #bfinal is in contracted form (if mask is defined)
  #cat("calling objfn1 for reporting purposes\n")
  bfinal.look <<- bfinal
  if(!silent)cat("dim TZ tz are:",dim(TZ),dim(tZ),"\n")
  ans <- objfn1(bfinal)
  ans$status <- "solved"
  if(fastreturn){
    #when doing simulations this is might be a good place to quit...
    #although monotonicity is not yet checked
    #I am pretty sure the objfn1() call just above with report==TRUE means
    #that gradmat.look contains cols whose sum is grad.raw
    grad.raw <- colSums(gradmat.look)
    ans$grad.raw <- grad.raw
    bmat.final <- matrix(bfinal,nr=ncol(Xs),nc=nYS)
    hpen <- penfunc3(bmat.final,Xs,Ygrid,SYgrid,sYgrid,ans$llfvec,pickmat,gradmode=FALSE) #now it is checked
    ans$dedyrange <- hpen$dedyrange
    return(ans)}
  systime2 <- Sys.time()
  #cat("At line 895 time spent in xnss2 is:",systime2-systime1,"\n")
  #cat("completed objfn1 for reporting purposes\n")
  ans.look <<- ans
  if(b.compressed){bmat.final <- vexpand(bfinal,buse=mask)}
  else{ bmat.final <- bfinal}
  bmat.final <- matrix(bmat.final,nr=ncol(Xs),nc=nYS)
  bmat.final.look <<- bmat.final
  slug7 <- system.time(
    hpen <- penfunc3(bmat.final,Xs,Ygrid,SYgrid,sYgrid,ans$llfvec,pickmat,gradmode=FALSE)
  )
  if(!silent)cat("penfunc3 call at 895 takes:",slug7,"\n")
  systime2 <- Sys.time()
  #cat("At line 898 time spent in xnss2 is:",systime2-systime1,"\n")
  #gradient without penalization; lam scoping does not permit doing it
  #this way from report==T and objfn1() call
  lamkeep <- lam
  lamkeep.look <<- lam
  lam <- 0
  grad.bare <- gr2(bfinal)
  lam <- lamkeep
  if(b.compressed) {grad.bare <- vexpand(grad.bare,buse=mask)}
  if(!b.compressed){grad.bare <- matrix(grad.bare,nrow=nXs,ncol=nYS)}
  ans$grad.bare <- grad.bare
  ans.look2 <<- ans
  #cat("expanding compressed b at line 792\n")
  if(b.compressed){bfinal <- vexpand(bfinal,buse=mask)
  }
  #get grad.raw, which is neither masked nor penalized
  mask.hold <- mask
  mask <- NULL
  lam <- 0
  gam.keep <- gam
  delta.keep <- delta
  gam <- 0
  delta <- 0
  b.compressed.hold <- b.compressed
  b.compressed <- FALSE
  grad.raw <- gr2(bfinal)
  b.compressed <- b.compressed.hold
  grad.raw <- matrix(grad.raw,nrow=nXs,ncol=nYS)
  gam <- gam.keep
  delta <- delta.keep
  #something like this is unneccessary: gradmat.raw <- gradmat.look
  #because gradmat is always raw; penalty only assessed on grad vector
  ans$grad.raw <- grad.raw
  lam <- lamkeep
  mask <- mask.hold
  #cat("at line 808\n")

  if(algor=="SCS"){ans.optim <- list(description="CVXR used")
  ans$algor <- "CVXR"
  }
  if(algor=="ECOS"){ans.optim <- list(description="ECOS used")
  ans$algor <- "ECOS"
  }
  if(maxit>0)ans$ans.optim <- ans.optim
  ans$Xs <- Xs
  report <- FALSE
  ans.look3 <<- ans

  #at this point bfinal is the 'expanded' form of the coefs if !is.null(mask)
  #here is one way to handle getting the Hessian (and J eventually))
  if(secalc|HJreturn){
    dedy <- function(b){
      dedy <- tZ%*%b
      dedy
    }
    if(b.compressed){
      btemp <- vcontract(bfinal,buse=mask)
      #Hess2 uses gr2 and so forth, all of which can see 'mask'
      #hfull <- Hess2(btemp)
      #hfull.look <<- hfull
      gradmat <- gradmat.look #cheating a bit...need to fix
      gradmat.unmasked <- gradmat.look
      gradmat <- cm_contract(gradmat,buse=mask) #why is this necessary?
      Jfull <- t(gradmat)%*%gradmat #add in Jfull ops too
      J1 <- Jfull
      h1 <- hfull  #although hfull is not 'full' but masked
      # Hunmasked <- Hess2.unmasked(btemp) #mask should be handled OK
      Hunmasked <- fastH1(list(TZ=TZ, tZ=tZ, dedy=dedy(btemp)))
    }

    else {
      #h1 <- Hess2(bfinal)
      btemp <- bfinal
      Hunmasked <- fastH1(list(TZ=TZ, tZ=tZ, dedy=dedy(btemp)))
      h1 <- Hunmasked
      gradmat <- gradmat.look #cheating a bit...need to fix
      J1 <- t(gradmat)%*%gradmat #add in Jfull ops too
    }



    #just in case (sic) things bomb, let's have the components
    #h1.look <<- h1
    #J1.look <<- J1
    #computing se's here
    if(secalc){
      ans$segg <- sqrt(diag(solve(J1,tol=1.e-40)))
      Hinv <- solve(h1,tol=1.e-40)
      if(!maxing)ans$seh <- sqrt(diag(Hinv))
      if(maxing)ans$seh <- sqrt(diag(-Hinv))
      ans$sew <- sqrt(diag(Hinv%*%J1%*%Hinv))
      ans$Hinv <- Hinv
    }
    #now going to send back gradmat, Hess, Hinv
    ans$gradmat <- gradmat
    ans$H <- h1
    ans$J <- J1
    ans.look4 <<- ans
    if(b.compressed){ans$Hunmasked <- Hunmasked}
    if(b.compressed){ans$gradmat.unmasked <-gradmat.unmasked}

  }

  #compute an alternative BIC based on 'effective' number of parameters
  nbeff <- sum(abs(bfinal)>bthresh)
  ans$bic2 <- -2*ans$llf+nbeff*log(length(y))
  ans$nbeff <- nbeff

  #another measure of complexity
  # if(lam==0){lam <- matrix(0,nrow=nXs,ncol=nYS)}
  D <- sum(2/(1+exp(lam)))
  ans$D <- D
  ans$Dbic <- -2*ans$llf+D*log(length(y))

  ans$TZ <- TZ
  ans$tZ <- tZ

  #discard coefs below the bthresh and recalc likelihood
  #probably need to review this section rhs Feb 2018
  getcleanb <- FALSE
  if(algor=="optim")getcleanb <- FALSE
  if(getcleanb){
    print("getting cleanb")
    bfinal <- matrix(bfinal,nrow=nXs,ncol=nYS)
    cleanb <- bfinal
    cleanb[abs(bfinal)<bthresh] <- 0
    cleanb[,1:2] <- bfinal[,1:2] #tricky
    report <- TRUE
    print("evaluating llf at cleanb")
    print(length(cleanb))
    cleanb.2 <- cleanb
    if(!is.null(mask)){cleanb.2 <- vcontract(cleanb,mask)}
    cleanllf <- sum(objfn1(cleanb.2)$llfvec)
    ans$cleanllf <- cleanllf
    ans$cleanb <- cleanb
  }

  #ans$h2 <- h2
  if(!is.null(mask)){ans$mask <- mask}
  ans$this.call <- sys.call()
  ans.look5 <<- ans
  ans$nu <- 1000
  ans$dtype <- "normal"  #NEEDED for display server to work properly!!!
  ans$xorth <- xorth
  ans$yform <- TRUE  #always now Feb 2018
  ans$yorth <- yorth
  ans$addxint <- addxint
  ans$yS <- yS
  ans$iyknots <- iyknots
  ans$ydf <- ydf #Jun 1 2018....much heartache before....
  ans$algor <- algor
  ans$diary <- diary
  ans$lam.orig <- lam.orig
  ans$lam <- gam*lam.mat  #NB lam return modifed RHS Nov 7 2017
  ans$delta <- delta
  ans$bthresh <- bthresh
  ans$reltol <- reltol
  ans$secalc <- secalc
  ans$silent <- silent
  ans$HJreturn <- HJreturn
  ans$checkg <- checkg
  ans$checkH <- checkH
  ans$date <- date()
  #if(!is.null(infom)) ans$infom <- infom
  if(any(isfactor)){ ans$infom <- infom}

  if(algor=="optim")ans$time.OPTIM <- time.optim
  if(algor=="SCS")ans$time.CVXR <- time.CVXR
  if(algor=="ECOS")ans$time.ECOS <- time.CVXR

  hpen3 <- penfunc3(ans$bmat,Xs,Ygrid,SYgrid,sYgrid,ans$llfvec,
                    pickmat,gradmode=TRUE)
  ans$hpen3 <- hpen3
  ans$ellf <- hpen3$ellf
  ans$overfit <- hpen3$overfit
  ans$dedyrange <- hpen3$dedyrange
  ans$ellfgrad <- hpen3$ellfgrad
  ans$gam <- gam
  ans$egam <- egam
  ans$delta <- delta

  ans$SYgrid <- SYgrid #this will make program testing easier Aug 25 2018
  ans$sYgrid <- sYgrid
  ans$Ygrid <- Ygrid
  ans$Ysing <- Ysing
  ans$btarg <- matrix(btarg,nr=nXs)
  #ans$DEDY <- DEDY
  #ans$HY2 <- HY2
  beta <- ans$beta
  dedygrid <- beta%*%t(sYgrid)
  egrid <- beta%*%t(SYgrid)
  phiegrid <- data_norm(egrid)
  fgrid <- phiegrid*dedygrid
  ans$dedygrid <- dedygrid
  ans$egrid <- egrid
  ans$fgrid <- fgrid

  systime3 <- Sys.time()
  if(!silent)cat("At line 1098 time spent in xnss2 is:",systime3-systime1,"\n")
  return(ans)
}
