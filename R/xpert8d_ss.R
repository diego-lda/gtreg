#' @title Xpert8d_ss
#'
#' @description Run the optimisation of the GTR problem using the CVXR package.
#'
#' @param TYX This is T(X,Y), the Kroneker product of the known functions of X and Y, W(X) and S(Y).
#' @param tYX This is the derivative vector of T(X,Y). Equivalent to t(X,Y) in the paper.
#' @param Kscore
#' @param gam
#' @param egam
#' @param lam
#' @param lam.vec
#' @param maxit This is the maximum number of iterations that will be run.
#' @param algor This is the algorithm that will be used to run the optimisation.
#' @param reltol
#' @param feastol
#' @param abstol
#' @param quiet
#' @param zeros
#' @param doprimal Whether or not the primal version of the problem should be run. Default is false.
#' @param btarg
#' @param silent
#' @param nXs
#' @param nYS
#' @param weights
#' @param cval
#' @param pen
#' @param beta2
#' @param Xs
#' @param sYgrid
#' @param bounded
#' @param Cbound
#' @param threshold This is the threshold level at which weights are eliminated. Default is 1e-5.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{llf} \tab The sum of all the elements in \code{llfvec}. \cr
#'    \tab \cr
#'    \code{e} \tab A numeric vector resulting from \eqn{\hat{b}`T(X,Y)}. \cr
#'    \tab \cr
#'    \code{eta} \tab  A numeric vector resulting from \eqn{\hat{b}`t(X,Y)}. \cr
#'    \tab \cr
#'    \code{finalscore} \tab The score of the problem: \eqn{T(X,Y)'e+t(X,Y)'eta}. \cr
#'    \tab \cr
#'    \code{result} \tab result from the CVXR optimisation. \cr
#'    \tab \cr
#'    \code{ehat} \tab The same as \code{e} if the primal form is run. \cr
#'    \tab \cr
#'    \code{etahat} \tab The same as \code{eta} if the primal form is run. \cr
#'    \tab \cr
#'    \code{h} \tab The number of observations. \cr
#'    \tab \cr
#'    \code{llfvec} \tab The numeric vector resulting from \code{log(dnorm(e)*(-1/h))}. \cr
#'    \tab \cr
#'    \code{lamx} \tab The numeric vector containing the penalisation parameters. \cr
#'    \tab \cr
#'    \code{bmat} \tab The numeric vector containing the parameters of interest, \eqn{\hat{b}}. \cr
#' }
#'
#' @import CVXR
#' @export
#'
#' @examples
xpert8d_ss <- function(TYX,tYX,Kscore=0,gam=0,egam=0,lam=0,lam.vec=NULL,maxit=200,algor="ECOS",
                       reltol=1.e-04,feastol=1.e-04,abstol=1.e-04,quiet=F,zeros=NULL,
                       doprimal=FALSE,btarg=0,silent=F,nXs=NULL,nYS=NULL,weights=1,cval=.1,pen=NULL,beta2=F,Xs=NULL,sYgrid=NULL,bounded=F,Cbound=1e6,threshold=1e-5){

  # tol.vec <<- c(reltol,feastol,abstol,doprimal,maxit)
  # TYX=datmat.l1$TZ;tYX=datmat.l1$tZ;Kscore=gam0;doprimal=F;pen=NULL;beta2=F;bounded=F;Cbound=1e6
  # gam <- gam0 ; egam<- 0 ; btarg=0
  # zeros <- which(weights<threshold)
  # if(length(zeros)>0){
  #   weights[zeros]  <- 0
  #   weights[-zeros] <- 1/as.numeric(weights[-zeros])
  # }
  # if(length(zeros)==0){
  #   weights <- 1/as.numeric(weights)
  # }

  #M is the 'mask matrix' M%*%b=0, M determined by mask
  #check that the constraint grid looks the way we think
  #print(dim(t(sYgrid)))
  M <- diag(1,nXs*nYS,nXs*nYS)
  if(length(zeros)>0){
    M            <- matrix(0,nXs*nYS,nXs*nYS-length(zeros))
    M[-zeros,]   <- diag(1,nXs*nYS-length(zeros),nXs*nYS-length(zeros))
  }

  #============CVXR algorithm calls======================
  #wholesale rewrite/simplification in prep for CVXR
  #doCVXR <- FALSE
  if(!doprimal){
    if(algor=="SCS"|algor=="ECOS"){
      if(algor=="SCS"){ SCS.dims_to_solver_dict <- CVXR:::SCS.dims_to_solver_dict } # This function is used for direct call to SCS
      bdim <- ncol(TYX)
      nobs <- nrow(TYX)
      if(!quiet)cat("Problem dimensions are:",nobs,bdim,"\n")

      if(length(Kscore)>1){Kscore <- matrix(Kscore*as.numeric(weights),nc=1)}
      if(gam>0 && mean(Kscore) == 0){ Kscore <- gam }
      scorebounds <- FALSE
      if(max(Kscore)>0){scorebounds <- TRUE}
      Kgauss <- log(1/sqrt(2*pi))

      e <- Variable(nobs)
      h <- Variable(nobs)
      DLLF  <- function(e,h){
        llfvec <- .5*e^2-log(-h)-1+Kgauss
        sum(llfvec)
      }
      score <- function(e,h){
        grad <- t(TYX)%*%(e)+t(tYX)%*%(h)
        grad
      }
      etafunc <- function(h){
        eta    <- inv_pos(-h)
        etadev <- eta-sum_entries(eta)/nobs
        eta
      }

      obj   <- DLLF(e,h)
      if(scorebounds){
        cond1a <- score(e,h)>= -Kscore
        cond1b <- score(e,h)<= Kscore
        constr <- list(cond1a,cond1b)
      }
      if(!scorebounds){
        cond1  <- score(e,h)==0
        constr <-list(cond1)
      }
      if(bounded){
        cond1  <- score(e,h)==0
        cond2  <- norm2(etafunc(h)) <= Cbound
        constr <- list(cond1,cond2)
      }
      prob <- Problem(Minimize(obj),constr)

      #if you want to use CVXR interface and/or check rules
      # 'unpack results' METHOD SEEMS TO NO LONGER WORK WITH CVXR version > "0.99-7".
      if(packageVersion("CVXR") <= "0.99-7"){
        doDCPcheck <- FALSE
        if(doDCPcheck){
          ecos.out0 <- psolve(prob,"ECOS",ignore_dcp=TRUE)
          result1   <- unpack_results(prob,"ECOS",ecos.out0)
        }
      }

      #let's try a direct call to the solver
      if(algor=="ECOS"){
        arglist <- ECOSolveR::ecos.control(maxit=as.integer(maxit),verbose=1L,reltol=reltol,feastol=feastol,abstol=abstol)
        if(quiet)arglist <- ECOSolveR::ecos.control(maxit=as.integer(maxit),verbose=0L,reltol=reltol,feastol=feastol,abstol=abstol)
        problem.data <- get_problem_data(prob,solver="ECOS")
        if(packageVersion("CVXR") > "0.99-7"){
          ECOS_dims <- ECOS.dims_to_solver_dict(problem.data$data[["dims"]])
          time.CVXR <- system.time(
            ecos.out1 <-  ECOSolveR::ECOS_csolve(
              c = problem.data$data[["c"]],
              G = problem.data$data[["G"]],
              h = problem.data$data[["h"]],
              dims = ECOS_dims,
              A = problem.data$data[["A"]],
              b = problem.data$data[["b"]],
              control = arglist
            )
          )

          result1 <- unpack_results(prob, ecos.out1, problem.data$chain, problem.data$inverse_data)
        }
        else{
          time.CVXR <- system.time(
            ecos.out1 <-  ECOSolveR::ECOS_csolve(
              c = problem.data[["c"]],
              G = problem.data[["G"]],
              h = problem.data[["h"]],
              dims = problem.data[["dims"]],
              A = problem.data[["A"]],
              b = problem.data[["b"]],control=arglist
            )
          )

          result1 <- unpack_results(prob, "ECOS", ecos.out1)

        }

      }

      if(algor=="SCS"){
        pd.scs   <- get_problem_data(prob,solver="SCS")
        arglist2 <- list(max_iters=maxit,eps=reltol)
        if(packageVersion("CVXR") > "0.99-7"){
          SCS_dims  <- SCS.dims_to_solver_dict(pd.scs$data[["dims"]])
          time.CVXR <- system.time(
            scs.out1  <- scs::scs(
              A = pd.scs$data[['A']],
              b = pd.scs$data[['b']],
              obj = pd.scs$data[['c']],
              cone = SCS_dims,
              control=arglist2
            )
          )

          result1   <- unpack_results(prob, scs.out1, pd.scs$chain, pd.scs$inverse_data)

        }
        else{
          time.CVXR <- system.time(
            scs.out1 <- scs::scs(
              A = pd.scs[['A']],
              b = pd.scs[['b']],
              obj = pd.scs[['c']],
              cone = pd.scs[['dims']],
              control=arglist2
            )
          )

          result1   <- unpack_results(prob, "SCS", scs.out1)
        }

      }

      if(result1$status=="solver_error" || result1$status=="unbounded_inaccurate"){ #solver failed, try to return gracefully

        ans <- list(status="failed")

        return(ans)

      }


      if(!quiet)print(result1$getValue(obj))
      e <- result1$getValue(e)
      h <- result1$getValue(h)


      finalscore <- score(e,h)
      objval     <- DLLF(e,h)
      llfvec     <- log(dnorm(e)*(-1/h))
      llf        <- sum(llfvec)
      if(!scorebounds){
        shadow1 <- result1$getDualValue(cond1)
        bhat    <- -shadow1
      }

      if(scorebounds){
        shadow1a <- result1$getDualValue(cond1a)
        shadow1b <- result1$getDualValue(cond1b)
        bhat     <- shadow1a-shadow1b
      }

      ehat   <- -(TYX)%*%bhat
      etahat <- tYX%*%bhat

      ans    <- list(llf=llf,e=e,eta=-1/h,finalscore=finalscore,
                     result=result1,ehat=ehat,etahat=etahat,
                     h=h,objval=objval,llfvec=llfvec,
                     bmat=bhat)

      ans$this.call <- sys.call()
      ans$time.CVXR <- time.CVXR
      #stop(44)
      return(ans)
    }
  }


  if(doprimal){

    bdim <- ncol(TYX)
    nobs <- nrow(TYX)
    if(!quiet)cat("Problem dimensions are:",nobs,bdim,"\n")

    if(length(Kscore)>1){ Kscore <- matrix(Kscore,nc=1) }
    scorebounds <- FALSE
    if(max(Kscore)>0){ scorebounds <- TRUE }

    Kgauss <- log(1/sqrt(2*pi))
    if(nXs==2){
      xmin <- min(TYX[,2])
      xmax <- max(TYX[,2])
    }

    b   <- Variable(bdim)

    LLF <- function(b){
      e      <- TYX%*%b
      dedy   <- tYX%*%b
      llfvec <- -.5*e^2+log(dedy)+Kgauss
      sum(llfvec)
    }

    dedy <- function(b){
      dedy <- tYX%*%b
      dedy
    }


    if(length(beta2)>0 || nYS==2){

      beta2X <- function(b){

        Xs.now <- NULL
        if(is.list(Xs)){
          for(kk in 1:length(Xs)){
            if(is.array(Xs[[kk]])){
              for(jj in 1:dim(Xs[[kk]])[3]){
                Xs.now <- rbind(Xs.now,Xs[[kk]][,,jj])
              }
            }
            if(!is.array(Xs[[kk]])){
              Xs.now <- rbind(Xs.now,Xs[[kk]])
            }
          }
        }

        if(!is.list(Xs)){ Xs.now <- Xs }

        if(packageVersion("CVXR") > "0.99-7") {
          Beta <- Xs.now%*%reshape_expr(M%*%b,c(nXs,nYS))
        }
        else{ Beta <- Xs.now%*%reshape_expr(M%*%b,nXs,nYS) }

        return(min(Beta[,2]))
      }

    }


    if(length(pen)>0){

      if(nYS==2){ dedy.grid <- beta2X }

      if(nYS!=2){

        if(nXs !=2){

          dedy.grid <- function(b){

            Xs.now <- NULL
            if(is.list(Xs)){
              for(kk in 1:length(Xs)){
                if(is.array(Xs[[kk]])){
                  for(jj in 1:dim(Xs[[kk]])[3]){
                    Xs.now <- rbind(Xs.now,Xs[[kk]][,,jj])
                  }
                }
                if(!is.array(Xs[[kk]])){
                  Xs.now <- rbind(Xs.now,Xs[[kk]])
                }
              }
            }

            if(!is.list(Xs)){ Xs.now <- Xs }

            if(packageVersion("CVXR") > "0.99-7") {
              Beta <- Xs.now%*%reshape_expr(M%*%b,c(nXs,nYS))
            }
            else{ Beta <- Xs.now%*%reshape_expr(M%*%b,nXs*nYS) }

            x <- Beta%*%t(sYgrid)  #x is dedyhat abbrev.
            a <- -40.
            return(min(x))

          }

        }

        if(nXs==2){

          dedy.grid  <- function(b){

            if(packageVersion("CVXR") > "0.99-7") {
              BetaY <- sYgrid%*%reshape_expr(M%*%b,c(nXs,nYS))
            }
            else{
              BetaY <- sYgrid%*%t(reshape_expr(M%*%b,nXs,nYS))
            }
            x1 <- BetaY[,1] + BetaY[,2]*xmin#,BetaY[,1] + BetaY[,2]*xmax)
            x2 <- BetaY[,1] + BetaY[,2]*xmax#,BetaY[,1] + BetaY[,2]*xmax)

            return(min(x1,x2))
          }

        }

      }

    }


    score <- function(e,eta){

      grad <- t(TYX)%*%(e)+t(tYX)%*%(eta)
      grad

    }



    if(Kscore>0){ lamx <- as.vector(Kscore) }
    if(Kscore==0 && gam>0){

      lamspecified <- is.matrix(lam)&is.null(lam.vec)

      if(!lamspecified){

        if(!is.matrix(lam)){ lam <- matrix(lam,nr=nXs,nc=nYS) }

        if(is.matrix(lam)){
          for(i in 1:length(lam)) lam[i] <- lam.vec[6] #gen
          lam[1,]  <- lam.vec[5] #row1
          lam[,1]  <- lam.vec[3] #col1
          lam[,2]  <- lam.vec[4] #col1
          lam[1,1] <- lam.vec[1] #int
          lam[1,2] <- lam.vec[2] #int
        }
      }

      if(length(zeros)>0){
        lamx <- as.vector(lam)[-zeros]*weights
      }
      if(length(zeros)==0){
        lamx <- as.vector(lam)*weights
      }

    }

    if(Kscore==0 && gam==0){ lamx <- 0 }
    #gam     <- 1 #tmp patch
    reg     <- gam*sum(lamx*abs(b))
    elastic <- egam*p_norm((b-btarg),2)
    obj     <- (LLF(b)-reg-elastic)

    if(length(pen)==0 && bounded){
      constraint.condns  <- norm2(b) <= Cbound
      constr             <- list(constraint.condns)
    }
    if(length(pen)==0 && beta2){
      constraint.condns  <- beta2X(b) >= cval
      constr             <- list(constraint.condns)
    }
    if(length(pen)>0 && !bounded && !beta2){
      constraint.condns <- dedy.grid(b)>= cval
      constr            <- list(constraint.condns)
    }
    if(length(pen)>0 && bounded && !beta2){
      constraint.condns1 <- dedy.grid(b) >= cval
      constraint.condns2 <- norm2(b) <= Cbound
      constr             <- list(constraint.condns1,constraint.condns2)
    }
    if(length(pen)>0 && !bounded && beta2){
      constraint.condns1 <- dedy.grid(b) >= cval
      constraint.condns2 <- beta2X(b) >= cval
      constr             <- list(constraint.condns1,constraint.condns2)
    }
    #just do Lasso
    #   obj <- (LLF(b)-reg)

    if(length(pen)==0 && !bounded && !beta2){ prob <- Problem(Maximize(obj)) }
    if(length(pen)>0 || bounded || beta2){  prob <- Problem(Maximize(obj),constr) }
    ecosolve.verbose <- as.integer(!silent)
    if(algor=="SCS"){
      #maxit <-5000
      try( solver.out1 <- solve(prob,solver=algor,verbose=!silent,max_iters=as.integer(maxit),eps=reltol) )
      #solver.out1$status
    }
    if(algor=="ECOS"){
      try( solver.out1 <- solve(prob,solver=algor,verbose=!silent,MAXIT=as.integer(maxit),RELTOL=reltol,FEASTOL=reltol,ABSTOL=reltol) )
    }
    if(class(solver.out1)=="try-error"){
      return(solver.out1)       #error-handling here
    }
    result1 <- solver.out1
    result1.look <<-  result1
    cat("at line 697 ...result1.look saved\n")

    result1.look <<- result1
    if(result1$status=="solver_error" || result1$status=="unbounded" || result1$status=="unbounded_inaccurate" || result1$status=="infeasible"){ #solver failed, try to return gracefully

      ans <- list(status="failed")

      return(ans)

    }

    bfinal <- result1$getValue(b)
    bfinal.look <<- bfinal
    cat("done primal")

    llf        <- LLF(bfinal)
    e          <- TYX%*%bfinal
    eta        <- tYX%*%bfinal
    finalscore <- score(e,eta)
    ehat       <- e
    etahat     <- eta
    h          <- NULL
    llfvec     <- log(dnorm(e)*eta)

    #not sending objval back for now...
    ans <- list(llf=llf,e=e,eta=eta,finalscore=finalscore,
                result=result1,ehat=ehat,etahat=etahat,
                h=h,llfvec=llfvec,lamx=lamx,
                bmat=bfinal
    )

    return(ans)

  }

}
