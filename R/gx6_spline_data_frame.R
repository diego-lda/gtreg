#' gX6 Data Frame (Spline)
#'
#' @description This is the spline dataframe function version for gX6.
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
gX6.spline.data.frame <- function(info,X,orth,xord=4){
  if(is.null(info$knots)){
    print("barf in gX6.spline")
    stop(82)
  }
  knots1 <- info$knots
  knots <- orthogonalsplinebasis::expand.knots(knots1)
  jnow <- info$coords
  center <- info$center

  Xnow <-X[jnow]
  if(!is.null(colnames(X))) nam <- colnames(X)[jnow] else nam <- paste0("var",jnow)

  if(class(Xnow[[1]])=="Date" | class(knots)=="Date" | class(center)=="Date"){
    knots <- as.numeric(knots)
    center <- as.numeric(center)
    Xnow[[1]] <- as.numeric(Xnow[[1]])
  }

  if(orth){
    spline <- SplineBasis(knots)
    spline <- orthogonalize(spline)
    tx <- eval_basis(object=spline,x=as.matrix(Xnow))
  }else{
    # spline <- SplineBasis(knots)
    # tx <- eval_basis(object=spline,x=Xnow)
    tx <- splineDesign(knots=knots,x=as.matrix(Xnow),outer.ok=TRUE,ord=xord) ## waaaay faster than eval_basis.
    #splineDesign patched out rhs tmp Sep 22 2017 for pqR compat.
  }

  if(!is.null(center)){
    if(!orth) spline <- SplineBasis(knots)
    cen_ev <- eval_basis(object=spline,x=center)
    for(i in 1:nrow(tx)){
      tx[i,] <- tx[i,]-cen_ev
    }
  }

  colnames(tx) <- lapply(1:dim(tx)[2L],FUN=function(i){paste0(nam,".",i)})
  tx <- as.data.frame(tx)

  tx <- tx[-1] #assume an intercept will be added
  if(is.null(info$shift)&is.null(info$scale)){return(tx)}

  ntx <- ncol(tx)
  if(is.null(info$shift))shift <- rep(0,ntx)
  else shift <- info$shift
  if(is.null(info$scale))scale <- rep(1,ntx)
  else scale <- info$scale
  for(jj in 1:ncol(tx)){
    tx[jj] <- shift[jj]+scale[jj]*tx[jj]
  }
  return(tx)
}
