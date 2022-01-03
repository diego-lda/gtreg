#' gX6 Data Frame (Main)
#'
#' @description This is the main dataframe function for gX6.
#'
#' @param X
#' @param info
#' @param orth
#'
#' @return
#' @export
#'
#' @examples
gX6.data.frame <- function(X,info,orth=FALSE){
  nops <- length(info)

  W <- X[NULL]

  for(i in 1:nops){
    type <- info[[i]]$type  #operation type, make life simple
    if(type=="bare"){
      whold <- gX6_bare.data.frame(info=info[[i]],X=X)
    }
    if(type=="spline"){
      whold <- gX6_spline.data.frame(info=info[[i]],X=X,orth=orth)
    }
    if(type=="tensor"){
      w  <- list()
      nw <- NULL
      nt <- length(info[[i]])-1     # info[[i]] includes info[[i]]$type="tensor", hence "-1".

      for(j in 1:nt){
        info.now <- info[[i]][[j]]
        if(info.now$type=="bare"){w[[j]] <- gX6_bare.data.frame(info=info.now,X=X)}
        if(info.now$type=="spline"){w[[j]] <- gX6_spline.data.frame(info=info.now,X=X,orth=orth)}
        nw <- c(nw,ncol(w[[j]]))
      }
      nam   <- paste0("tens",i)

      nwhold <- nw[1]
      w1.now <- w[[1]]
      for(ll in 1:(nt-1)){
        jdex  <- 0
        whold <- as.data.frame(matrix(0,nr=nrow(w[[1]]),ncol=nwhold*nw[ll+1]))
        for(jj in 1:nwhold){
          for(kk in 1:nw[ll+1]){
            jdex <- jdex+1
            whold[jdex] <- w1.now[jj]*w[[ll+1]][kk]
          }
        }
        nwhold <- ncol(whold)
        w1.now <- whold
      }

      colnames(whold) <- lapply(1:dim(whold)[2L],FUN=function(i){paste0(nam,".",i)})
    }

    W <- cbind(W,whold)

  }
  return(W)
}
