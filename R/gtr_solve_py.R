#' GTR Solver Using python Optimisation
#'
#' @description This function solves the optimisation problem using the python package.
#'
#' @param TYX
#' @param tYX
#'
#' @return The dictionary with the optimisation results.
#' @export
#'
#' @examples
gtr_solve_py <- function(TYX, tYX) {
  if (!reticulate::py_available()) {
    stop("Python is not available. Please configure Python using reticulate::use_python('/path/to/python').")
  }
  # Import your Python package
  gtregpy <- import("gtregpy")

  # Call the gtr_solve function from the Python package
  result <- gtregpy$gtr_solve(TYX, tYX)

  # Return the result to R
  return(result)
}
