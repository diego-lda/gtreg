#' gX6 Matrix (Main)
#'
#' @description This is the main matrix function version for gX6.
#'
#' @param X
#' @param info
#' @param orth
#'
#' @return
#' @export
#'
#' @examples
gX6.matrix <- function(X,info,orth=FALSE){
  nops <- length(info)  #number of operations to be carried out

  # js <- 1
  W <- NULL
  Xsingleton <- FALSE
  if(!is.matrix(X)){Xsingleton <- TRUE
  X <- matrix(X,nr=1)}


  for(i in 1:nops){
    type <- info[[i]]$type  #operation type, make life simple
    if(type=="bare"){
      whold <- gX6_bare.matrix(info=info[[i]],X=X)
    }
    if(type=="spline"){
      whold <- gX6_spline.matrix(info=info[[i]],X=X,orth=orth)
    }
    if(type=="tensor"){
      #for now we allow only 2-d tensor
      info1 <- info[[i]]$info1
      info2 <- info[[i]]$info2
      if(info1$type=="bare"){w1 <- gX6_bare.matrix(info=info1,X=X)}
      if(info1$type=="spline"){w1 <- gX6_spline.matrix(info=info1,X=X,orth=orth)}
      if(info2$type=="bare"){w2 <- gX6_bare.matrix(info=info2,X=X)}
      if(info2$type=="spline"){w2 <- gX6_spline.matrix(info=info2,X=X,orth=orth)}

      nam <- paste0("tens",i)

      #these 2 statements should not be necessary!!
      #fix? fix!
      if(!is.matrix(w1)){w1 <- matrix(w1,nr=1)}
      if(!is.matrix(w2)){w2 <- matrix(w2,nr=1)}
      nw1 <- ncol(w1)
      nw2 <- ncol(w2)
      jdex <- 0
      whold <- matrix(0,nr=nrow(w1),ncol=nw1*nw2)
      for(jj in 1:nw1){
        for(kk in 1:nw2){
          jdex <- jdex+1
          whold[,jdex] <- w1[,jj]*w2[,kk]
        }
      }
      colnames(whold) <- lapply(1:dim(whold)[2L],FUN=function(i){paste0(nam,".",i)})
    } #tensor type ends

    #print(dim(W));print(dim(whold))

    if(Xsingleton){whold <- matrix(whold,nr=1)}
    W <- cbind(W,whold)
  }
  return(W)
}
