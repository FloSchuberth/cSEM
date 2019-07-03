#' Internal: Calculate matrix difference using distance measure
#'
#' Calculates the differences between symmetric matrices using a given 
#' distance measure. This is typically used to calculate the difference between
#' the model-implied and the empirical indicator covariance matrix.
#' C
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


#' Get parameter names
#' 
#' Based on a cSEmodelDetermines the names of the parameters to be used for testing 
#' 
#' @usage getParameterNames(.object = NULL, .model = NULL)
#' 
#' @inheritParams csem_arguments
#' @param .model A model in [lavaan model syntax][lavaan::model.syntax] indicating which 
#'   parameters (i.e, path (`~`), loadings (`=~`), or weights (`<~`)) should be
#'   compared across groups. Defaults to `NULL` in which case all parameters of the model
#'   are compared.
#'   
#' @return A list containing the names of the structural parameters,
#'   the loadings and weight to be compared.
#'   
#' @keywords internal
#' 
getParameterNames <- function(
  .object  = NULL,
  .model   = NULL
){
  
  ## Summarize
  x <- summarize(.object)
  
  if(inherits(.object, "cSEMResults_2ndorder")) {
    
    x12 <- x[[1]]$First_stage$Information
    x22 <- x[[1]]$Second_stage$Information
    
    # Extract different types of constructs
    construct_type <- c(x12$Model$construct_type, x22$Model$construct_type)
    construct_type <- construct_type[unique(names(construct_type))]
    measurement_org <- x22$Arguments_original$.model$measurement
    
  } else {
    
    x22  <- x[[1]]$Information
    # Extract different types of constructs
    construct_type <- x22$Model$construct_type
    measurement_org <- x22$Model$measurement
  }
  
  # Parse model that indicates which parameters should be compared
  # if no model indicating the comparisons is provided, all parameters are compared
  # This prevents the the test_MGD function to break down if no comparison model is supplied.
  
  if(is.null(.model)) {
    if(inherits(.object, "cSEMResults_2ndorder")) {
      model_comp <- x22$Arguments_original$.model
    } else {
      model_comp <- x22$Model 
    }
  } else {
    model_comp <- parseModel(.model, .check_errors = FALSE)
  }
  
  
  # Check whether the constructs specified in the comparison are equal to the constructs in the original model
  construct_type_comp=model_comp$construct_type[!is.na(model_comp$construct_type)]
  
  if(!all(names(construct_type_comp)%in%names(construct_type))){
    stop2("At least one construct appears in the comparison model and not in the original model.")
  }
  
  if(!all(construct_type_comp==construct_type[names(construct_type_comp)])){
    stop2("At least one construct's type in the comparison model differs from the original model.")
  }
  
  # Check wehther the indicators used in the comparison model also appear in the original model
  if(!all(colnames(model_comp$measurement)%in%colnames(measurement_org))){
    stop2("At least one indicator appears in the comparison model and not in the original model.")
  }

  ### Extract names ============================================================
  # Extract names of the path to be tested
  temp <- outer(rownames(model_comp$structural), colnames(model_comp$structural), 
                FUN = function(x, y) paste(x, y, sep = " ~ "))
  
  names_path <- t(temp)[t(model_comp$structural) != 0]
  
  if(length(names_path) == 0) {
    names_path <- NULL
  }
  
  ## Extract names of the loadings to be tested
  names_row <- rownames(model_comp$measurement)
  
  ## Select only concepts modeled as common factors. Reorder to be have the 
  ## same order as names_row.
  i <- intersect(names(which(construct_type == "Common factor")), names_row)
  i <- i[match(names_row, i)]
  i <- i[!is.na(i)]
  
  if(length(i) != 0) {
    temp <- rep(rownames(model_comp$measurement[i, , drop = FALSE]), 
                times = rowSums(model_comp$measurement[i, , drop = FALSE]))
    if(length(temp)!=0){
    temp_n <- model_comp$measurement[i, colSums(model_comp$measurement[i, , drop = FALSE]) != 0, drop = FALSE]
    names_loadings <- paste0(temp, " =~ ", colnames(temp_n)) 
    }else{
      names_loadings = NULL
    }
  } else {
    names_loadings <- NULL
  }
  
  ## Extract names of the weights to be tested
  # Select only concepts modeled as common factors. Reorder to be have the 
  # same order as names_row.
  i <- intersect(names(which(construct_type == "Composite")), names_row)
  i <- i[match(names_row, i)]
  i <- i[!is.na(i)]
  
  if(length(i) != 0) {
    temp <- rep(rownames(model_comp$measurement[i, , drop = FALSE]), 
                times = rowSums(model_comp$measurement[i, , drop = FALSE]))
    if(length(temp)!=0){
    temp_n <- model_comp$measurement[i, colSums(model_comp$measurement[i, , drop = FALSE]) != 0, drop = FALSE]
    names_weights <- paste0(temp, " <~ ", colnames(temp_n))
    }else{
      names_weights = NULL
    } 
  } else {
    names_weights <- NULL
  }
  
  ## Return as list
  out <- list(
    "names_path"     = names_path, 
    "names_weights"  = names_weights,
    "names_loadings" = names_loadings
    )
  return(out)
}

#' Parameter differences across groups
#' 
#' Calculates the difference between one or more paramater estimates based 
#' one different groups (data sets).
#'
#' @usage calculateParameterDifference(
#'   .object     = NULL,
#'   .model = args_default()$.model
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
  .model = args_default()$.model
  ){

  ## Summarize
  x <- summarize(.object)

  ## Get names to compare
  names <- getParameterNames(.object, .model = .model)
  names_path     <- names$names_path
  names_loadings <- names$names_loadings
  names_weights  <- names$names_weights
  
  ### Compute differences ======================================================
  if(inherits(.object, "cSEMResults_2ndorder")) {
    path_estimates  <- lapply(x, function(y) {y$Second_stage$Estimates$Path_estimates})
    loading_estimates <- lapply(x, function(y) {
      rbind(y$First_stage$Estimates$Loading_estimates, 
            y$Second_stage$Estimates$Loading_estimates)
      })
    weight_estimates <- lapply(x, function(y) {
      rbind(y$First_stage$Estimates$Weight_estimates, 
            y$Second_stage$Estimates$Weight_estimates)
      })
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

#' ANOVA F-test statistic 
#'
#' Calculates the ANOVA F-test statistic suggested by Sarstedt et al. (2011) 
#' 
#' @usage calculateFR <- function(.Parameter,
#' .id)
#' 
#' @inheritParams csem_arguments
#' 
#' @return A named scaler, the test statistic of the ANOVA F-test
#'
#' @keywords internal

calculateFR <- function(.Parameter,
                        .id){
  ParameterIdMatrix = cbind(.Parameter,.id)
  
  G <- length(unique(.id))
  
  Agbar <- sapply(1:G, function(x){
    mean(.Parameter[which(.id == unique(.id)[x])])
  })
  
  B <- nrow(ParameterIdMatrix)
  Abar <- mean(.Parameter)
  SSbetween <- G * B *(1/(G-1)) * sum((Agbar - Abar)^2) 
  
  SSwithin <- 1/(B-1) * sum((.Parameter - Agbar)^2)
  
  SSbetween/SSwithin
}