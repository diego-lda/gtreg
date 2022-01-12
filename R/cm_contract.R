#' Contribution Matrix Contraction Function
#'
#' @description This contracts a 'contributions' (n by big) matrix by the mask in buse, by removing the appropriate columns.
#' @param cm The contribution matrix
#' @param buse A matrix with the Betas to use.
#'
#' @return
#' @export
#'
#' @examples
cm_contract <- function(cm,buse){

  usevec <- as.logical(buse)
  cm2 <- cm[,usevec]
  return(cm2)
}
