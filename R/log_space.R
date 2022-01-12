#' Log Space Function
#'
#' @description This function does a log transformation onto the log space.
#'
#' @param a
#' @param b
#' @param n
#'
#' @return
#' @export
#'
#' @examples
logspace   <- function( a, b, n){

  return(exp(log(10)*seq(a, b, length.out=n)))
}
