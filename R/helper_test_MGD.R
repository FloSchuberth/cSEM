#' Internal: Calculate matrix difference using distance measure
#'
#' Calculates the differences between symmetric matrices using a given 
#' distance measure. This is typically used to calculate the difference between
#' the model-implied and the empirical indicator covariance matrix.
#' 
#' `.matrices` must be a list of at least two matrices. If more than two matrices 
#' are supplied the arithmetic mean over all possible matrix 
#' distances is computed. Hence, supplying a large number of matrices will 
#' quickly become computationally challenging. Use with care.
#' 
#' Currently two distance measures are supported:
#' \describe{
#'   \item{`geodesic`}{The geodesic distance}
#'   \item{`squared_euclidian`}{The squared Euclidian distance}
#' }
#' 
#' @usage calculateDistance(
#'   .matrices = args_default()$.matrices, 
#'   .distance = args_default()$.distance
#'   )
#' 
#' @inheritParams csem_arguments
#' 
#' @return A numeric vector of length one containing the (arithmetic mean) of 
#'   the differences between all possible combinations of matrices.
#' @keywords internal

calculateDistance <- function(
  .matrices = args_default()$.matrices, 
  .distance = args_default()$.distance
){
  ### Checks and errors ========================================================
  ## Check if list and at least of length 2
  if(!is.list(.matrices) && length(.matrices) < 2) {
    stop2("`.matrices` must be a list of at least length two.")
  }

  ## Check if matrices are all symmetric
  if(!all(sapply(.matrices, matrixcalc::is.symmetric.matrix))) {
    stop2("All matrices in `.matrices` must be symmetric.")
  }
  
  ### Calculation ==============================================================
  ## Combine all matrices into lists of two
  temp <- utils::combn(.matrices, 2, simplify = FALSE)
  
  ## Compute the distance measure for each group combination
  distances <- lapply(temp, function(x) {
    switch (.distance,
            "geodesic" = {calculateDG(.matrix1 = x[[1]], .matrix2 = x[[2]])},
            "squared_euclidian" = {calculateDL(.matrix1 = x[[1]], .matrix2 = x[[2]])}
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


#' Parameter differences across groups
#' 
#' Calculates the difference between one or more paramater estimates based 
#' one different groups (data sets).
#'
#' @usage calculateParameterDifference(
#'   .object     = NULL,
#'   .comparison = args_default()$.comparison
#' )
#' 
#' @inheritParams csem_arguments
#' 
#' @return A list of length equal to the number of possible combinations of
#'   groups in .object (basically n choose k = 2), e.g., 3 if there are three 
#'   groups and 6 if there are 4 groups. Each list elements contains the
#'   values of the difference between the parameter estimates based on the
#'   data of group i and and group j.
#' @keywords internal
#' 
calculateParameterDifference <- function(
  .object     = NULL,
  .comparison = args_default()$.comparison
  ){
  
  ## Summarize
  x <- summarize(.object)
  
  if(inherits(.object, "cSEMResults_2ndorder")) {
    
    x22 <- x[[1]]$Second_stage$Information
    
  } else {

    x22  <- x[[1]]$Information
  }
  

  # Parse model that indicates which parameters should be compared
  # if no model indicating the comparisons is provided, all parameters are compared
  # This prevents the the test_MGD function to break down if no comparison model is supplied.
  
  if(is.null(.comparison)) {
    model_comp <- x22$Model  
  } else {
    model_comp <- parseModel(.comparison, .check_errors = FALSE)
  }
  
  # # Extract different types of constructs
  construct_type <- x22$Model$construct_type
  
  ### Extract names ============================================================
  # Extract names of the path to be tested
  temp <- outer(rownames(model_comp$structural), colnames(model_comp$structural), 
                FUN = function(x, y) paste(x, y, sep = " ~ "))
  
  names_path <- t(temp)[t(model_comp$structural) != 0]
  
  # Extract names of the loadings to be tested
  names_row <- rownames(model_comp$measurement)
  names_col <- colnames(model_comp$measurement)
  
  i <- intersect(names(which(construct_type == "Common factor")), names_row)
  
  temp <- rep(rownames(model_comp$measurement[i, , drop = FALSE]), times = rowSums(model_comp$measurement[i, , drop = FALSE]))
  names_loadings <- paste0(temp, " =~ ", colnames(model_comp$measurement[i, colSums(model_comp$measurement[i, , drop = FALSE]) != 0, drop = FALSE]))
  
  # Extract names of the weights to be tested
  i <- intersect(names(which(construct_type == "Composite")), names_row)
  temp <- rep(rownames(model_comp$measurement[i, , drop = FALSE]), times = rowSums(model_comp$measurement[i, , drop = FALSE]))
  names_weights <- paste0(temp, " <~ ", colnames(model_comp$measurement[i, colSums(model_comp$measurement[i, , drop = FALSE]) != 0, drop = FALSE]))
  
  ### Compute differences ======================================================
  if(inherits(.object, "cSEMResults_2ndorder")) {
    path_estimates  <- lapply(x, function(y) {y$Second_stage$Estimates$Path_estimates})
    loading_estimates <- lapply(x, function(y) {y$Second_stage$Estimates$Loading_estimates})
    weight_estimates <- lapply(x, function(y) {y$Second_stage$Estimates$Weight_estimates})
  } else {
    path_estimates  <- lapply(x, function(y) {y$Estimates$Path_estimates})
    loading_estimates <- lapply(x, function(y) {y$Estimates$Loading_estimates})
    weight_estimates <- lapply(x, function(y) {y$Estimates$Weight_estimates})
  }
  
  ## Select
  path_estimates  <- lapply(path_estimates, function(y) {
    y1 <- y[y$Name %in% names_path, "Estimate"]
    if(length(y1) != 0) {
      names(y1) <- names_path
    }
    y1
  })
  loading_estimates <- lapply(loading_estimates, function(y) {
    y1 <- y[y$Name %in% names_loadings, "Estimate"] 
    if(length(y1) != 0) {
      names(y1) <- names_loadings
    }
    y1
  })
  weight_estimates <- lapply(weight_estimates, function(y) {
    y1 <- y[y$Name %in% names_weights, "Estimate"]
    if(length(y1) != 0) {
      names(y1) <- names_weights
    }
    y1
  })
  
  ## Path model
  temp <- utils::combn(path_estimates, 2, simplify = FALSE)
  diff_path <- lapply(temp, function(y) y[[1]] - y[[2]])
  names(diff_path) <- sapply(temp, function(x) paste0(names(x)[1], '_', names(x)[2]))
  
  ## Loadings
  temp <- utils::combn(loading_estimates, 2, simplify = FALSE)
  diff_loadings <- lapply(temp, function(y) y[[1]] - y[[2]])
  names(diff_loadings) <- sapply(temp, function(x) paste0(names(x)[1], '_', names(x)[2]))
  
  ## Weights
  temp <- utils::combn(weight_estimates, 2, simplify = FALSE)
  diff_weights <- lapply(temp, function(y) y[[1]] - y[[2]])
  names(diff_weights) <- sapply(temp, function(x) paste0(names(x)[1], '_', names(x)[2]))

  # merge list together
  out <- mapply(function(x,y,z) c(x,y,z),
             x = diff_path,
             y = diff_loadings,
             z = diff_weights,
             SIMPLIFY = FALSE)

  return(out)
}



#' Multiple testing correction
#'
#' Adjust a given significance level .alpha to accomodate multiple testing. 
#' 
#' @usage adjustAlpha <- function(
#'  .alpha                 = args_default()$.alpha,
#'  .approach_alpha_adjust = args_default()$.approach_alpha_adjust,
#'  .nr_comparisons        = args_default()$.nr_comparisons
#' )
#' 
#' @inheritParams csem_arguments
#' 
#' @return A vector of (adjusted) significance levels.
#'
#' @keywords internal
adjustAlpha <- function(
  .alpha                 = args_default()$.alpha,
  .approach_alpha_adjust = args_default()$.approach_alpha_adjust,
  .nr_comparisons        = args_default()$.nr_comparisons
  ){
  
  if(.approach_alpha_adjust == 'none'){
    return(.alpha)
  }
  
  if(.approach_alpha_adjust == 'bonferroni'){
    return(.alpha/.nr_comparisons)
  }
  
}
