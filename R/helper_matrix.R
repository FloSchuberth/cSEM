#' Internal: Trace of a matrix
#' 
#' Compute the trace of `X`, i.e., the sum of the diagonal elements of `X`.
#'
#' @param X Square numeric matrix.
#'
#' @return A numeric scalar: value with the trace of `X`.
#'
#' @keywords internal
trace <- function(X) { # trace of a matrix
  if (!is.matrix(X))
    stop2("X is not a matrix")
  else if (nrow(X) != ncol(X))
    stop2("X must be square!")

  sum(diag(X))
}


#' Internal: Check if matrix is positive semi definite
#' 
#' Check if a matrix is positive semi definite (PSD). 
#' A matrix is PSD if all of its eigenvalues are greater than or
#' equal to zero. One caveat, the package previously used
#' `matrixcalc::is.positive.semi.definite()` where negative values
#' close to zero `(abs(eigenvalues) < tol)` were not counted as negative.
#' This behaviour is replicated here as well.
#'
#' @param X Numeric matrix to be checked.
#' @param zero.tol Tolerance for close to zero (negative) values.
#'
#' @return A logical scalar: `TRUE` if `X` is PSD, otherwise `FALSE`.
#'
#' @keywords internal
isPositiveSemiDefinite <- function(X, zero.tol = 1e-08) {
  # This matches the old behaviour of `matrixcalc::is.positive.semi.definite()`
  # With some additional checks for NA/NaN/Inf values, and error handling of eigen()

  if (!is.matrix(X))
      stop2("X is not a matrix")
  else if (!isSymmetric(X)) 
      stop2("X is not a symmetric matrix")
  else if (!is.numeric(X)) 
      stop2("X is not a numeric matrix")

  ev <- tryCatch(
    eigen(X, only.values = TRUE, symmetric = TRUE)$values, # faster with symmetric = TRUE
    error = function(e) NA
  )

  # If any (negative) values are very close to zero, they are not treated as
  # negative but as zero values. This was in the original code done by replacing
  # abs(ev) < zero.tol with 0. This is a simpler approach yielding the same results.
  !any(ev < -zero.tol | !is.finite(ev))
}
