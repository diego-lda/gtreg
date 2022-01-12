#' Evaluation of Basis Functions.
#'
#' @description This function evluates the basis functions. It is based closely on evaluate method for spline basis in package orthogonalsplinebasis, but about 30% faster.
#'
#' @param object
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
eval_basis <- function(object,x,...) {
  stopifnot(is.numeric(x))
  dots <- list(...)
  mat <- object@Matrices
  knots <- object@knots
  order <- object@order
  M <-lapply(1:dim(mat)[3L],FUN=function(i) matrix(mat[,,i],nrow=order))
  x[toolow <- x < knots[order]] <- knots[order]
  x[toohigh <- x > knots[length(knots)-order+1]] <- knots[length(knots)-order+1]
  results <-t(vapply(x,eb_helper,numeric(dim(mat)[2L]),M=M,knots=knots,order=order))
  # Bit of a patchy solution
  results[(toolow | toohigh),] <- 0
  return(results)
}
