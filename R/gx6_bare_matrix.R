#' gX6 Matrix (Bare)
#'
#' @description This is the bare matrix function version for gX6.
#'
#' @param info
#' @param X
#'
#' @return
#' @export
#'
#' @examples
gX6_bare.matrix <- function(info,X){
  jnow <- info$coords
  Xnow <- if(jnow=="_intercept") rep(1,nrow(X)) else as.numeric(X[,jnow])
  if(!is.null(info$scale)) Xnow <- scale*Xnow
  if(!is.null(info$shift)) Xnow <- Xnow + shift
  Xnow <- matrix(Xnow, nc=1)
  if(!is.null(info$center)){
    if(class(Xnow[[1]])=="Date"){
      Xnow <- as.Date(as.numeric(Xnow) - as.numeric(info$center))
    }else{
      Xnow <- Xnow - info$center
    }
  }
  return(Xnow)
}
