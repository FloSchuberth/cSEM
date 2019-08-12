#' Internal: Matrix difference
#'
#' Calculates the average of the differences between all possible pairs of 
#' (symmetric) matrices in a list using a given distance measure.
#' 
#' `.matrices` must be a list of at least two matrices. If more than two matrices 
#' are supplied the arithmetic mean of the differences between all possible pairs of 
#' (symmetric) matrices in a list is computed. Mathematically this is
#' n chose 2. Hence, supplying a large number of matrices will 
#' become computationally challenging.
#' 
#' Currently two distance measures are supported:
#' \describe{
#'   \item{`geodesic`}{(Default) The geodesic distance.}
#'   \item{`squared_euclidian`}{The squared Euclidian distance}
#' }
#' 
#' @usage calculateDistance(
#'   .matrices = NULL, 
#'   .distance = args_default()$.distance
#'   )
#' 
#' @inheritParams csem_arguments
#' 
#' @return A numeric vector of length one containing the (arithmetic) mean of 
#' the differences between all possible pairs of matrices supplied via `.matrices`.
#' 
#' @keywords internal

calculateDistance <- function(
  .matrices = NULL, 
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



#' Internal: Parameter names
#' 
#' Based on a model in [lavaan model syntax][lavaan::model.syntax], returns the 
#' names of the parameters of the structural
#' model, the measurement/composite model and the weight relationship. Used 
#' by [testMGD()] to extract the names of the parameters to compare across groups
#' according to the test proposed by \insertCite{Chin2010;textual}{cSEM}. 
#' 
#' @usage getParameterNames(
#'            .object  = args_default()$.object,
#'            .model   = args_default()$.model
#'   )
#' 
#' @inheritParams csem_arguments
#' @param .model A model in [lavaan model syntax][lavaan::model.syntax] indicating which 
#'   parameters (i.e, path (`~`), loadings (`=~`), or weights (`<~`)) should be
#'   compared across groups. Defaults to `NULL` in which case all parameters of the model
#'   are compared.
#'   
#' @return A list with elements `names_path`, `names_loadings`, and `names_weights`
#' containing the names of the structural parameters, the loadings, 
#' and the weight to compare across groups.
#' 
#' @references
#'   \insertAllCited{}
#'   
#' @keywords internal

getParameterNames <- function(
  .object  = args_default()$.object,
  .model   = args_default()$.model
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
    names_path_org <- x[[1]]$Second_stage$Estimates$Path_estimates$Name
    
  } else {
    
    x22  <- x[[1]]$Information
    # Extract different types of constructs
    construct_type <- x22$Model$construct_type
    
    measurement_org <- x22$Model$measurement
    names_path_org <- x[[1]]$Estimates$Path_estimates$Name
  }
  
  # Parse model that indicates which parameters should be compared.
  # If no model indicating the comparisons is provided, all parameters are compared.
  
  if(is.null(.model)) {
    if(inherits(.object, "cSEMResults_2ndorder")) {
      model_comp <- x22$Arguments_original$.model
    } else {
      model_comp <- x22$Model 
    }
  } else {
    model_comp <- parseModel(.model, .check_errors = FALSE)
  }
  
  # Check whether the constructs specified in the comparison are equal
  # to the constructs in the original model
  if(!all(rownames(model_comp$structural) %in% rownames(measurement_org))){
    stop2(
      "The following error occured in the `getParameterNames()` function:\n",
      "At least one construct appears in the comparison model but",
      " not in the original model."
      )
  }
  
  # Check whether the construct type specified in the comparison are equal
  # to the construct types in the original model
  construct_type_comp <- model_comp$construct_type[!is.na(model_comp$construct_type)]  
  if(!all(construct_type_comp == construct_type[names(construct_type_comp)])){
    stop2(
      "The following error occured in the `getParameterNames()` function:\n",
      "At least one construct's type in the comparison model differs ", 
      "from the original model.")
  }
  
  # Check wehther the indicators used in the comparison model also appear 
  # in the original model
  if(!all(colnames(model_comp$measurement) %in% colnames(measurement_org))){
    stop2(
      "The following error occured in the `getParameterNames()` function:\n",
      "At least one indicator appears in the comparison model and not in ",
      " the original model.")
  }
  
  # Check whether the indicators used are correctly assigned, i.e., 
  # as in the original model
  lapply(rownames(model_comp$measurement), function(x) {
    names_ind_comp <- colnames(model_comp$measurement[x,,drop=FALSE][,which(model_comp$measurement[x,,drop=FALSE]==1),drop=FALSE])
    names_ind_org  <- colnames(measurement_org[x,,drop=FALSE][,which(measurement_org[x,,drop=FALSE]==1),drop=FALSE])
    
    if(!all(names_ind_comp %in% names_ind_org)){
      stop2(
        "The following error occured in the `getParameterNames()` function:\n",
        "At least one indicator is not correctly assigned in the comparison model.")
    }
  })
  
  ### Extract names ============================================================
  # Extract names of the path to be tested
  temp <- outer(rownames(model_comp$structural), colnames(model_comp$structural), 
                FUN = function(x, y) paste(x, y, sep = " ~ "))
  
  names_path <- t(temp)[t(model_comp$structural) != 0]
  
  if(length(names_path) == 0) {
    names_path <- NULL
  }
  
  # Check whether the path to compare occur also in the original mode 
  if(!all(names_path %in% names_path_org)){
    stop2(
      "The following error occured in the `getParameterNames()` function:\n",
      "At least one of the paths specified for comparison does not ",
      " appear in the original model.")
  }
  
  # Extract names of the loadings to be tested
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
      names_loadings <- NULL
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
      names_weights <- NULL
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



#' Internal: Parameter differences across groups
#' 
#' Calculate the difference between one or more paramater estimates across
#' all possible pairs of groups (data sets) in `.object`.
#'
#' @usage calculateParameterDifference(
#'   .object     = args_default()$.object,
#'   .model      = args_default()$.model
#' )
#' 
#' @inheritParams csem_arguments
#' @param .model A model in [lavaan model syntax][lavaan::model.syntax] indicating which 
#'   parameters (i.e., path (`~`), loadings (`=~`), or weights (`<~`)) should be
#'   compared across groups. Defaults to `NULL` in which case all parameters of the model
#'   are compared.
#' 
#' @return A list of length equal to the number of possible pairs of
#'   groups in `.object` (mathematically, this is n choose 2, i.e., 3 if there are three 
#'   groups and 6 if there are 4 groups). Each list elements is itself a list of
#'   three. The first list element contains
#'   the difference between parameter estimates of the structural model, the second
#'   list element the difference between estimated loadings, and the third
#'   the difference between estimated weights.
#'   
#' @keywords internal

calculateParameterDifference <- function(
  .object  = args_default()$.object,
  .model   = args_default()$.model
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



#' Internal: Multiple testing correction
#'
#' Adjust a given significance level `.alpha` to accomodate multiple testing.
#' The following corrections are implemented:
#' \describe{
#'   \item{`none`}{(Default) No correction is done.}
#'   \item{`bonferroni`}{A Bonferroni correction is done, i.e., alpha is divided by the 
#'   number of comparisons `.nr_comparisons`.}
#' }
#' 
#' @usage adjustAlpha(
#'  .alpha                 = args_default()$.alpha,
#'  .approach_alpha_adjust = args_default()$.approach_alpha_adjust,
#'  .nr_comparisons        = args_default()$.nr_comparisons
#' )
#' 
#' @inheritParams csem_arguments
#' 
#' @return A vector of (possibly adjusted) significance levels.
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

#' Internal: ANOVA F-test statistic 
#'
#' Calculate the ANOVA F-test statistic suggested by 
#' \insertCite{Sarstedt2011;textual}{cSEM} in the OTG testing procedure.
#' 
#' @usage calculateFR(.resample_sarstedt)
#' 
#' @inheritParams csem_arguments
#' 
#' @return A named scaler, the test statistic of the ANOVA F-test
#'
#' @references
#'   \insertAllCited{}
#'   
#' @keywords internal

calculateFR <- function(.resample_sarstedt) {
  
  # In Sarstedt et al. (2011), it is assumed that you have the same number of bootstrap estimates per group. 
  # Consequently, the total number of estimates equals B*G.
  # Since, there might be different numbers of bootstrap estimates per group, we adjusted the calculation of the 
  # F test statistic.
  
  names_param <- colnames(.resample_sarstedt)[-ncol(.resample_sarstedt)]
  
  # Number of grousp
  G <- length(unique(.resample_sarstedt[, "group_id"]))
  
  # Total number of estimates 
  N <- nrow(.resample_sarstedt)
  
  # Order the estimates by the grouping variable
  x <- .resample_sarstedt[order(.resample_sarstedt[, "group_id"]), ]

  grand_means <- colMeans(x[, names_param, drop = FALSE])
  mean_group  <- aggregate(x[, names_param, drop = FALSE], by = list(x[, "group_id"]), FUN = mean)
  
  group_mean_matrix <- apply(mean_group[, -1, drop = FALSE], 2, function(y) rep(y, times = table(x[, "group_id"])))
  
  SS_between  <- rowSums((t(group_mean_matrix) - grand_means)^2)
  MSA_between <- SS_between/ (G -1)
  SS_within   <- colSums((x[, names_param, drop = FALSE] - group_mean_matrix)^2)
  MSA_within  <- SS_within/(N - G)
  # Note: in the original paper (B -1) was used. This is incorrect. It should
  # be B - G (G = the number of groups). See the F-statistic derived from ANOVA
  # for comparision.
  
  MSA_between/MSA_within 
  
}

#' Internal: Calculation of the CDF used in Henseler et al. (2009) 
#'
#' Calculates the probability that theta^1 is smaller than or equal to theta^2. 
#' See Equation (6) in \insertCite{Sarstedt2011;textual}{cSEM}.
#' 
#' @usage calculatePr(.resample_centered = NULL, .parameters_to_compare = NULL)
#' 
#' @inheritParams csem_arguments
#' 
#' @return A named vector
#'
#' @references
#'   \insertAllCited{}
#'   
#' @keywords internal
calculatePr <- function(
  .resample_centered       = NULL,
  .parameters_to_compare   = NULL
  ){
  
  # Remove names from .parameters_to_compare
  names(.parameters_to_compare) <- NULL
  
  group1 <- .resample_centered[[1]][,.parameters_to_compare,drop=FALSE]
  
  group2 <- .resample_centered[[2]][,.parameters_to_compare,drop=FALSE]
  
  ret <- sapply(.parameters_to_compare, function(x){
    # Matrix where each element of the first vector is substracted from the complete second vector 
    # theta^1 <= theta^2
    temp12 <- outer(group2[,x],group1[,x],"-")
    g1leqg2<-mean((1+sign(temp12))/2)
    # theta^2 <= theta^1 is equal to 1 - g1leqg2, therefore we do not need the following calculations
    # temp21 <- outer(group1[,x],group2[,x],"-") 
    # g2leqg1<-mean((1+sign(temp21))/2)
  
    # # The problem is that Henseler's approach is a one-sided test. 
    # # Therefore, the p-value is flipped if it larger than 0.5
    # # Return
    # if(g1leqg2<0.5){
    #   g1leqg2
    # }else{
    #   1-g1leqg2
    # }
    })
  names(ret) <- .parameters_to_compare
  ret
  
}