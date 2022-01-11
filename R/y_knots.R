#' Y Knots Function
#'
#' @description This function generates the knots for the dependent variable.
#'
#' @param y
#' @param e0mode
#' @param ydf
#' @param yorder
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
y_knots <- function(y,e0mode=F,ydf,yorder,delta=1){
  if(e0mode){ iyknots <- quantile(y,probs=seq(0,1,length=ydf)) }
  if(!e0mode){
    yabsmax <- max(abs(y))
    knoteps <- 1.e-5
    iyknots    <- qnorm(seq(knoteps,1-knoteps,length=ydf))
    iyknots[1] <- -delta*yabsmax
    iyknots[length(iyknots)] <- delta*yabsmax
  }

  return(iyknots)
}
