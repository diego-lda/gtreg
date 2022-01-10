#' CM Contract Function
#'
#' @description This contracts a 'contributions' (n by big) matrix by the mask in buse, by removing the appropriate columns.
#' @param cm
#' @param buse
#'
#' @return
#' @export
#'
#' @examples
cmcontract <- function(cm,buse){

  usevec <- as.logical(buse)
  cm2 <- cm[,usevec]
  return(cm2)
}
