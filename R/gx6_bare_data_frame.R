#' gX6 Data Frame (Bare)
#'
#' @description This is the bare dataframe function version for gX6.
#'
#' @param info
#' @param X
#'
#' @return
#' @export
#'
#' @examples
gX6_bare.data.frame <- function(info,X){
  jnow <- info$coords
  Xnow <- if(jnow=="_intercept") data.frame(int=rep(1,nrow(X))) else X[jnow]
  if(is.factor(Xnow[[1]])) return(Xnow)
  if(!is.null(info$scale)) Xnow <- scale*Xnow
  if(!is.null(info$shift)) Xnow <- Xnow + shift
  if(!is.null(info$center)){
    if(class(Xnow[[1]])=="Date"){
      Xnow[1] <- as.Date(as.numeric(Xnow[[1]]) - as.numeric(info$center))
    }else{
      Xnow <- Xnow - info$center
    }
  }
  return(Xnow)
}
