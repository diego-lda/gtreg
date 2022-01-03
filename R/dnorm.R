#' DNorm (v2)
#'
#' @description This function normalises the data.
#'
#' @param x
#' @param mean
#' @param sd
#' @param log
#'
#' @return
#' @export
#'
#' @examples
dnorm2 <- function(x,mean=0,sd=1,log=FALSE){
  ans <- 0.3989422804014327*exp(-0.5*((x-mean)/sd)^2)/sd
  return(ans)
}
