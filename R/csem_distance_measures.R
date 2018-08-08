#' Internal: Calculate matrix difference using distance measure
#'
#' Calculates the differences between symmetric matrices using a given 
#' distance measure. This is typically used to calculate the difference between
#' the model-implied and the empirical indicator covariance matrix.
#' 
#' `.matrices` must be a list of at least two matrices. If more than two matrices 
#' are supplied the arithmetic mean distance over all possible matrix 
#' distances is computed. Hence, supplying a large number of matrices will 
#' quickly become computationally challenging. Use with care.
#' 
#' Currently two distance measures are supported:
#' 
#' \describe{
#'   \item{`geodesic`}{The geodesic distance. Details here.}
#'   \item{`squared_euclidian`}{The squared euclidian distance. Details here.}
#' }
#' 
#' @usage calculateDistance(
#'   .matrices = args_default()$.matrices, 
#'   .distance = args_default()$.distance
#'   )
#' 
#' @inheritParams csem_arguments
#' 
#' @return A numeric vector of length one containing the (arithmetic mean) difference 
#'   between all matrices.
#' @keywords internal

calculateDistance <- function(
  .matrices = args_default()$.matrices, 
  .distance = args_default()$.distance
  ){
  ### Checks and errors ========================================================
  ## Check if list and at least of length 2
  if(!is.list(.matrices) && length(.matrices) < 2) {
    stop("`.matrices` must be a list of at least length two.", call. = FALSE)
  }
  # ## Check if all elements of matrices are of class matrix and symmetric
  # if(!all(sapply(.matrices, is.matrix))) {
  #   stop("All elements of `.matrices` must be matrices.", call. = FALSE)
  # }
  if(!all(sapply(.matrices, matrixcalc::is.symmetric.matrix))) {
    stop("All matrices in `.matrices` must be symmetric.", 
         call. = FALSE)
  }

  ### Calculation ==============================================================
  ## Combine all matrices into lists of two
  temp <- utils::combn(.matrices, 2, simplify = FALSE)
  
  ## Compute the distance measure for each group combination
  distances <- lapply(temp, function(x) {
    switch (.distance,
            "geodesic" = {dG(.matrix1 = x[[1]], .matrix2 = x[[2]])},
            "squared_euclidian" = {dL(.matrix1 = x[[1]], .matrix2 = x[[2]])}
    )
  })
  
  ## Compute the mean. Since dG and dL are defined as the mean of difference between
  ## two matrices (1/2 is part of the definition). This must be corrected if 
  ## more than two matrices are to be compared.
  
  if(length(.matrices) > 2) {
    out <- mean(2 * unlist(distances))
  } else {
    out <- unlist(distances)
  }
  return(out)
}

# Distance functions

### Squared euclidean distance

dL <- function(
  .matrix1 = args_default()$.matrix1,
  .matrix2 = args_default()$.matrix2
) {
  ## Check if dimensions are identical
  if(!identical(dim(.matrix1), dim(.matrix2))) {
    stop("`.matrix1` and `.matrix2` must have the same dimension.",
         call. = FALSE)
  }
  
  ## Calculate distance
  0.5 * sum((.matrix1 - .matrix2)[lower.tri(.matrix1, diag = FALSE)]^2)
}

### Geodesic distance

dG <- function(
  .matrix1 = args_default()$.matrix1,
  .matrix2 = args_default()$.matrix2
) {
  ## Check if dimensions are identical
  if(!identical(dim(.matrix1), dim(.matrix2))) {
    stop("`.matrix1` and `.matrix2` must have the same dimension.",
         call. = FALSE)
  }
  # not sure if logarithm naturalis is used or logarithm with base 10. 
  Eigen            <- eigen(solve(.matrix1) %*% .matrix2)
  logEigenvaluessq <- (log(Eigen$values, base = 10))^2   

  ## Calculate distance
  0.5 * sum(logEigenvaluessq)
}

### dML: the fitting function used in FIML

dML <- function(.object=args_default()$.object){
  
  nobs <- .object$Information$Number_of_observations
  S         <- .object$Estimates$Indicator_VCV
  p         <- dim(S)[1]
  Sigma_hat <- fit(.object)
  (nobs - 1)*(log(det(Sigma_hat)) + sum(diag(S %*% solve(Sigma_hat))) - log(det(S)) - p)
}