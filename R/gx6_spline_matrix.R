#' gX6 Matrix (Spline)
#'
#' @description Transforms the data in X according to the instructions in info. Sees if there is a spline; if not just apply shift, scale.
#'
#' @param info
#' @param X
#' @param orth
#' @param xord
#'
#' @return
#' @export
#'
#' @examples
gX6_spline.matrix <- function(info,X,orth=FALSE,xord=4){
  if(is.null(info$knots)){
    print("barf in gX6.spline")
    stop(82)
  }
  knots1 <- info$knots
  knots <- expand.knots(knots1)
  jnow <- info$coords
  center <- info$center

  if(jnow > ncol(as.matrix(X))){
    print("warning: skipping coord outside X")
    return(NULL)
  }
  Xnow <- as.matrix(X)[,jnow]
  if(!is.null(colnames(X))) nam <- colnames(X)[jnow] else nam <- paste0("var",jnow)

  if(class(Xnow)=="Date" | class(knots)=="Date" | class(center)=="Date"){
    knots <- as.numeric(knots)
    center <- as.numeric(center)
    Xnow <- as.numeric(Xnow)
  }

  if(orth){
    spline <- SplineBasis(knots)
    spline <- orthogonalize(spline)
    tx <- EvalBasis(object=spline,x=Xnow)
  }else{
    # spline <- SplineBasis(knots)
    # tx <- EvalBasis(object=spline,x=Xnow)
    tx <- splineDesign(knots=knots,x=Xnow,outer.ok=TRUE,ord=xord) ## waaaay faster than EvalBasis.
    #splineDesign patched out rhs tmp Sep 22 2017 for pqR compat.
  }

  if(!is.null(center)){
    if(!orth) spline <- SplineBasis(knots)
    cen_ev <- EvalBasis(object=spline,x=center)
    for(i in 1:nrow(tx)){
      tx[i,] <- tx[i,]-cen_ev
    }
  }

  # spline <- SplineBasis(knots)
  # if(orth)spline <- orthogonalize(spline)
  # jnow <- info$coords
  #jnow is always (???) scalar
  # Xnow <- X[,jnow]
  #if(length(Xnow)==1){Xnow <- as.vector(Xnow)}
  #tx <- evaluate(object=spline,x=Xnow)
  # tx <- EvalBasis(object=spline,x=Xnow)

  colnames(tx) <- lapply(1:dim(tx)[2L],FUN=function(i){paste0(nam,".",i)})

  tx <- as.matrix(tx[,-1]) #assume an intercept will be added
  if(is.null(info$shift)&is.null(info$scale)){return(tx)}

  ntx <- ncol(tx)
  if(is.null(info$shift))shift <- rep(0,ntx)
  else shift <- info$shift
  if(is.null(info$scale))scale <- rep(1,ntx)
  else scale <- info$scale
  for(jj in 1:ncol(tx)){
    tx[,jj] <- shift[jj]+scale[jj]*tx[,jj]
  }
  return(tx)
}
