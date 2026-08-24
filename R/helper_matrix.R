#' Internal: Trace of a matrix
#' 
#' Compute the trace of `.matrix`, i.e., the sum of the diagonal elements of `.matrix`.
#'
#' @param .matrix Square numeric matrix.
#'
#' @return A numeric scalar: value with the trace of `.matrix`.
#'
#' @keywords internal
trace <- function(.matrix) {
  if (!is.matrix(.matrix))
    stop2("Input is not a matrix")
  else if (nrow(.matrix) != ncol(.matrix))
    stop2("Matrix must be square!")

  sum(diag(.matrix))
}


#' Internal: Check if matrix is positive semi definite
#' 
#' Check if a matrix is positive semi definite (PSD). 
#' A matrix is PSD if all of its eigenvalues are greater than or
#' equal to zero.
#'
#' @param .matrix Numeric matrix to be checked.
#' @param .tolerance Tolerance for close to zero (negative) values.
#'
#' @return A logical scalar: `TRUE` if `X` is PSD, otherwise `FALSE`.
#'
#' @keywords internal
isPositiveSemiDefinite <- function(.matrix, .tolerance = 1e-08) {
  # This matches the old behaviour of `matrixcalc::is.positive.semi.definite()`
  # With some additional checks for NA/NaN/Inf values, and error handling of eigen()

  if (!is.matrix(.matrix))
      stop2("Input is not a matrix")
  else if (!isSymmetric(.matrix)) 
      stop2("Matrix is not a symmetric matrix")
  else if (!is.numeric(.matrix)) 
      stop2("Matrix is not a numeric matrix")

  ev <- tryCatch(
    eigen(.matrix, only.values = TRUE, symmetric = TRUE)$values, # faster with symmetric = TRUE
    error = function(e) NA
  )

  # If any (negative) values are very close to zero, they are not treated as
  # negative but as zero values. This was in the original code done by replacing
  # abs(ev) < zero.tol with 0. This is a simpler approach yielding the same results.
  !any(ev < -.tolerance | !is.finite(ev))
}
