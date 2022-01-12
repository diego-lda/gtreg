#' Evaluation of Basis Functions Helper
#'
#' @description This function helps to evluate the basis functions.
#'
#' @param x
#' @param M
#' @param knots
#' @param order
#'
#' @return
#' @export
#'
#' @examples
EB <- function(x,M,knots,order){
  if(x==knots[length(knots)-order+1]){ ind <- x<=knots } else { ind <- x<knots }
  if(all(ind)|all(!ind))  {
    if(x==knots[length(knots)-order+1])
      return(rep(1,order)%*%M[[length(M)]])
    else
      return(rep(0,dim(M)[2L]))
  }
  i<-.Internal(which(ind))[1L]-1L
  u<-(x-knots[i])/(knots[i+1]-knots[i])
  U<-u^(0:(order-1))
  return(U%*%M[[i-order+1L]])
}
