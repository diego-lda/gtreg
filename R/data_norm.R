#' Data Normalisation
#'
#' @description This function normalises the data.
#'
#' @param x The data to be normalised.
#' @param mean The value to center on.
#' @param sd The standard deviation to target.
#' @param log Whether to log the data. NOT IMPLEMENTED
#'
#' @return It returns a normalised version of the data.
#' @export
#'
#' @examples
data_norm <- function(x,mean=0,sd=1,log=FALSE){
  ans <- 0.3989422804014327*exp(-0.5*((x-mean)/sd)^2)/sd
  return(ans)
}
