#' TIC Function
#'
#' @description
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
#'
#' @examples
TIC.func <- function(res,TZ,tZ,Ginv=T,Adex=NULL,method="TIC"){

  if(length(Adex)==0){ Adex <- 1:ncol(TZ) }
  TZ       <- TZ[,Adex]
  tZ       <- tZ[,Adex]
  bmat     <- res$bmat[Adex]
  ifelse(length(Adex)==1, e <- TZ*bmat, e <- TZ%*%bmat)
  ifelse(length(Adex)==1, eta <- tZ*bmat, eta <- tZ%*%bmat)
  ans <- NULL

  if(method=="TIC" || method=="both"){

    if(min(eta)>0){

      llf      <- sum(log(dnorm(e)*eta))
      gradmat  <- -TZ*as.vector(e)+tZ*as.vector(1/eta)
      grad2    <- t(gradmat)%*%gradmat
      hessmat  <- -(t(TZ)%*%TZ+t(tZ/as.vector(eta))%*%(tZ/as.vector(eta)))
      Hinv     <- NULL
      ifelse(length(Ginv)==0, try(Hinv <- solve(-hessmat,tol=1e-40), silent=T), Hinv <- ginv(-hessmat,tol=1e-50))

      bias <- sum(diag(Hinv%*%grad2))

      ans  <- c(ans,-2*llf + 2*bias)

    }
    else{ ans <- c(ans,1e22) }

  }

  if(method=="BIC" || method=="both"){

    if(min(eta)>0){

      llf <- sum(log(dnorm(e)*eta))
      df  <- ncol(TZ)
      ans <- c(ans,-2*llf+df*log(length(y)))

    }
    else{ ans <- c(ans,1e22) }

  }

  return(ans)

}
