#' V Contract Function
#'
#' @description This gives the gradient pieces for b for the optimisation algorithm
#'
#' @param g
#' @param buse
#'
#' @return Gradient pieces.
#' @export
#'
#' @examples
vcontract <- function(g,buse){

  usevec <- as.logical(buse)
  g2 <- g[usevec]
  return(g2)
}
