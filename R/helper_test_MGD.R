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
#'   \item{`squared_euclidean`}{The squared Euclidian distance}
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
            "squared_euclidean" = {calculateDL(.matrix1 = x[[1]], .matrix2 = x[[2]])}
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
    
    indicators <- x12$Model_original$indicators
    cons_exo <- x12$Model_original$cons_exo
     
    names_cor_exo_cons_org <- x[[1]]$Second_stage$Estimates$Exo_construct_correlation$Name
    names_cor_indicator_org <- c(x[[1]]$First_stage$Estimates$Indicator_correlation$Name,
                                 x[[1]]$Second_stage$Estimates$Indicator_correlation$Name)
    
    if(nrow(x[[1]]$First_stage$Estimates$Residual_correlation)!=0){
      temp <- x[[1]]$First_stage$Estimates$Residual_correlation$Name
    }else{
      temp <- NULL
    }
    
    if(nrow(x[[1]]$Second_stage$Estimates$Residual_correlation)!=0){
      temp1 <- x[[1]]$Second_stage$Estimates$Residual_correlation$Name
    }else{
      temp1 <- NULL
    }
    names_cor_measurement_error_org <- c(temp,temp1)
    
    
  } else {#if no second-order model
    
    x22  <- x[[1]]$Information
    # Extract different types of constructs
    construct_type <- x22$Model$construct_type
    
    measurement_org <- x22$Model$measurement
    names_path_org <- x[[1]]$Estimates$Path_estimates$Name
    
    indicators <- x22$Model$indicators
    cons_exo <- x22$Model$cons_exo
    
    names_cor_exo_cons_org <- x[[1]]$Estimates$Exo_construct_correlation$Name
    names_cor_indicator_org <- x[[1]]$Estimates$Indicator_correlation$Name
    
    if(nrow(x[[1]]$Estimates$Residual_correlation)!=0){
      names_cor_measurement_error_org <- x[[1]]$Estimates$Residual_correlation$Name
    }else{
      names_cor_measurement_error_org<-NULL
    }
  }
  
  # FOR FUTURE:
  # Perhaps it is better to transfer the following code to the testMGD function directly.
  # As getParameterNames should just return the names of a model and not internally select a model?
  
  # Parse model that indicates which parameters should be compared.
  # If no model indicating the comparisons is provided, all parameters are compared.
  # Correlation among the measurement errors and or indicators cannot be compared yet across groups as 
  # theta is not part of the csem output 

  
  # FOR FUTURE: Add all model parameters to the comparison if no model is specified 
  
  if(is.null(.model)) {
    if(inherits(.object, "cSEMResults_2ndorder")) {
      model_compare <- x22$Arguments_original$.model
      # add the correlations among the exogenous constructs
    } else {
      model_compare <- x22$Model
      # add the correlations among the exogenous constructs
      # add to model_cor_specified and expand
      # Attention if the some of the variables are already in that matrix
    }
  } else {
    model_compare <- parseModel(.model, .check_errors = FALSE)
  }
  

  # Check whether the constructs specified in the comparison are equal
  # to the constructs in the original model
  if(!all(rownames(model_compare$structural) %in% rownames(measurement_org))){
    stop2(
      "The following error occured in the `getParameterNames()` function:\n",
      "At least one construct appears in the comparison model but",
      " not in the original model."
      )
  }
  
  # Check whether the construct type specified in the comparison are equal
  # to the construct types in the original model
  construct_type_compare <- model_compare$construct_type[!is.na(model_compare$construct_type)]  
  if(!all(construct_type_compare == construct_type[names(construct_type_compare)])){
    stop2(
      "The following error occured in the `getParameterNames()` function:\n",
      "At least one construct's type in the comparison model differs ", 
      "from the original model.")
  }
  
  # Check wehther the indicators used in the comparison model also appear 
  # in the original model
  if(!all(colnames(model_compare$measurement) %in% colnames(measurement_org))){
    stop2(
      "The following error occured in the `getParameterNames()` function:\n",
      "At least one indicator appears in the comparison model and not in ",
      " the original model.")
  }
  
  # Check whether the indicators used are correctly assigned, i.e., 
  # as in the original model
  lapply(rownames(model_compare$measurement), function(x) {
    names_ind_compare <- colnames(model_compare$measurement[x,,drop=FALSE][,which(model_compare$measurement[x,,drop=FALSE]==1),drop=FALSE])
    names_ind_org  <- colnames(measurement_org[x,,drop=FALSE][,which(measurement_org[x,,drop=FALSE]==1),drop=FALSE])
    
    if(!all(names_ind_compare %in% names_ind_org)){
      stop2(
        "The following error occured in the `getParameterNames()` function:\n",
        "At least one indicator is not correctly assigned in the comparison model.")
    }
  })
  
  ### Extract names ============================================================
  # Extract names of the path to be tested
  temp <- outer(rownames(model_compare$structural), colnames(model_compare$structural), 
                FUN = function(x, y) paste(x, y, sep = " ~ "))
  
  names_path_compare <- t(temp)[t(model_compare$structural) != 0]
  
  if(length(names_path_compare) == 0) {
    names_path_compare <- NULL
  }
  
  # Check whether the path to compare occur also in the original mode 
  if(!all(names_path_compare %in% names_path_org)){
    stop2(
      "The following error occured in the `getParameterNames()` function:\n",
      "At least one of the paths specified for comparison does not ",
      " appear in the original model.")
  }
  
  
  # all correlated variables, i.e., that have been specified with ~~
  vars_correlated_compare <- rownames(model_compare$cor_specified)
  
  # Select only those that are indicators
  ind_correlated_compare <- intersect(vars_correlated_compare, indicators)
  
  # Extract names of the loadings and measurement error to be tested
  names_row <- rownames(model_compare$measurement)
  
  ## Select only concepts modeled as common factors. Reorder to be have the 
  ## same order as names_row.
  i <- intersect(names(which(construct_type == "Common factor")), names_row)
  i <- i[match(names_row, i)]
  i <- i[!is.na(i)]
  
  if(length(i) != 0) {
    temp <- rep(rownames(model_compare$measurement[i, , drop = FALSE]), 
                times = rowSums(model_compare$measurement[i, , drop = FALSE]))
    if(length(temp)!=0){
      # Loadings
    temp_n <- model_compare$measurement[i, colSums(model_compare$measurement[i, , drop = FALSE]) != 0, drop = FALSE]
    names_loadings_compare <- paste0(temp, " =~ ", colnames(temp_n)) 
    # Measurement error
    measurement_error_correlated <- intersect(ind_correlated_compare,colnames(temp_n))
    cor_measurement_error <- model_compare$cor_specified[ measurement_error_correlated, measurement_error_correlated]

    index <- which(cor_measurement_error == 1, arr.ind = TRUE)
    names_correlated_measurement_error_compare <- index
    
    if(nrow(names_correlated_measurement_error_compare) ==0 ){
      names_correlated_measurement_error_compare <- NULL
    }else{
      names_correlated_measurement_error_compare <- paste(rownames(cor_measurement_error)[index[,"row"]],
                                                    " ~~ ",
                                                    colnames(cor_measurement_error)[index[,"col"]],
                                                    sep="")
      
      # Select those that have the same order as in the original model
      names_correlated_measurement_error_compare <- names_correlated_measurement_error_compare[names_correlated_measurement_error_compare %in% names_cor_measurement_error_org]
    }
    
    }else{
      names_loadings_compare <- NULL
      names_correlated_measurement_error_compare <- NULL
    }
  } else {
    names_loadings_compare <- NULL
    names_correlated_measurement_error_compare <- NULL
  }
  
  ## Extract names of the weights and indicator correlations (composites) to be tested
  # Select only concepts modeled as composite. Reorder that it has the 
  # same order as names_row.
  ## Indicator correlation, i.e., variables that belong to a composite
  i <- intersect(names(which(construct_type == "Composite")), names_row)
  i <- i[match(names_row, i)]
  i <- i[!is.na(i)]
  
  if(length(i) != 0) {
    temp <- rep(rownames(model_compare$measurement[i, , drop = FALSE]), 
                times = rowSums(model_compare$measurement[i, , drop = FALSE]))
    if(length(temp)!=0){
      # weights
    temp_n <- model_compare$measurement[i, colSums(model_compare$measurement[i, , drop = FALSE]) != 0, drop = FALSE]
    names_weights_compare <- paste0(temp, " <~ ", colnames(temp_n))
    
    # indicator correlations
    ind_connected_to_composite <- colnames(temp_n)
    indicator_correlated <- intersect(ind_correlated_compare,ind_connected_to_composite)
    
    cor_indicator <- model_compare$cor_specified[indicator_correlated, indicator_correlated]
    
    
    index <- which(cor_indicator == 1, arr.ind = TRUE)
    names_correlated_indicator_compare <- index
    
    # In case that no indicator correlations are compared set it to NULL
    if(nrow(names_correlated_indicator_compare) ==0 ){
      names_correlated_indicator_compare <- NULL
    }else{
      names_correlated_indicator_compare <- paste(rownames(cor_indicator)[index[,"row"]],
                                            " ~~ ",
                                            colnames(cor_indicator)[index[,"col"]],
                                            sep="")
      
      # Select those that have the same order as in the original model
      names_correlated_indicator_compare <- names_correlated_indicator_compare[names_correlated_indicator_compare %in% names_cor_indicator_org]
    }
    }else{
      names_weights_compare <- NULL
      names_correlated_indicator_compare <- NULL
    } 
  } else {
    names_weights_compare <- NULL
    names_correlated_indicator_compare <- NULL
  }
  

  ## Exogenous construct correlations 
  cons_exo_correlated_compare <- intersect(vars_correlated_compare, cons_exo)

  # Consider only exogenous constructs
  cor_cons_exo <- model_compare$cor_specified[cons_exo_correlated_compare,cons_exo_correlated_compare]

  # Which are correlated
  index <- which(cor_cons_exo == 1,arr.ind = TRUE) 
  names_correlated_exo_cons_compare <- index
  
  # In case that no measurement error correlations are compared set it to NULL
  if(nrow(names_correlated_exo_cons_compare) == 0){
    names_correlated_exo_cons_compare <- NULL 
  }else{
    names_correlated_exo_cons_compare <-paste(rownames(cor_cons_exo)[index[,"row"]],
                              " ~~ ",
                              colnames(cor_cons_exo)[index[,"col"]], sep="")
    # Select those that have the same order as in the original model
    names_correlated_exo_cons_compare <- names_correlated_exo_cons_compare[names_correlated_exo_cons_compare %in% names_cor_exo_cons_org]
    
  }
  
  ## Return as list
  out <- list(
    "names_path"     = names_path_compare, 
    "names_weights"  = names_weights_compare,
    "names_loadings" = names_loadings_compare,
    "names_cor_exo_cons" = names_correlated_exo_cons_compare,
    "names_cor_indicator" = names_correlated_indicator_compare,
    "names_cor_measurement_error" = names_correlated_measurement_error_compare
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
  names_cor_exo_cons <- names$names_cor_exo_cons
  # names_cor_measurement_error <- names$names_cor_measurement_error
  
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
    cor_cons_exo_estimates <- lapply(x, function(y) {
      y$Second_stage$Estimates$Exo_construct_correlation})
    
  } else { # no second-order model
    path_estimates  <- lapply(x, function(y) {y$Estimates$Path_estimates})
    loading_estimates <- lapply(x, function(y) {y$Estimates$Loading_estimates})
    weight_estimates <- lapply(x, function(y) {y$Estimates$Weight_estimates})
    # all exogenous construct correlations
    cor_cons_exo_estimates <- lapply(x, function(y) {
      y$Estimates$Exo_construct_correlation
    })
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

  cor_cons_exo_estimates <- lapply(cor_cons_exo_estimates, function(y) {
    y1 <- y[y$Name %in% names_cor_exo_cons, "Estimate"]
    if(length(y1) != 0) {
      names(y1) <- names_cor_exo_cons
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

  # Exogenous construct correlations
  temp <- utils::combn(cor_cons_exo_estimates, 2, simplify = FALSE)
  cor_cons_exo <- lapply(temp, function(y) y[[1]] - y[[2]])
  names(cor_cons_exo) <- sapply(temp, function(x) paste0(names(x)[1], '_', names(x)[2]))
  
  # merge list together
  out <- mapply(function(w,x,y,z) c(w,x,y,z),
             w = diff_path,
             x = diff_loadings,
             y = diff_weights,
             z = cor_cons_exo,
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

#' Internal: Extract relevant parameters from several cSEMResults_multi
#' 
#' Extract the relevant parameters from a cSEMResult_multi object in `.object`.
#'
#' @usage getRelevantParameters(
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
#' @return A list of length equal to the number of groups in `.object`.
#'  Each list element is itself a list of three. The first list element contains
#'   the relevant parameter estimates of the structural model, the second
#'   list element the relevant estimated loadings, and the third
#'   the relevant estimated weights.
#'   
#' @keywords internal

getRelevantParameters <- function(
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
  names_cor_exo_cons <- names$names_cor_exo_cons
  # names_cor_measurement_error <- names$names_cor_measurement_error
  
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
    cor_cons_exo_estimates <- lapply(x, function(y) {
      y$Second_stage$Estimates$Exo_construct_correlation})
    
  } else { # no second-order model
    path_estimates  <- lapply(x, function(y) {y$Estimates$Path_estimates})
    loading_estimates <- lapply(x, function(y) {y$Estimates$Loading_estimates})
    weight_estimates <- lapply(x, function(y) {y$Estimates$Weight_estimates})
    # all exogenous construct correlations
    cor_cons_exo_estimates <- lapply(x, function(y) {
      y$Estimates$Exo_construct_correlation
    })
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
  
  cor_cons_exo_estimates <- lapply(cor_cons_exo_estimates, function(y) {
    y1 <- y[y$Name %in% names_cor_exo_cons, "Estimate"]
    if(length(y1) != 0) {
      names(y1) <- names_cor_exo_cons
    }
    y1
  })
  # Do not change the order
  list(path_estimates = path_estimates,
       loading_estimates = loading_estimates,
       weight_estimates = weight_estimates,
       cor_exo_cons_estimates = cor_cons_exo_estimates)
}



#' Internal: get structured cSEMTestMGD results
#'
#' Convenience function to summarize the results of all tests resulting from a 
#' call to [testMGD()] in a user-friendly way.
#' 
#' @usage structureTestMGDDecisions(.object)
#'
#' @inheritParams csem_arguments
#' 
#' @return A data.frame.
#' 
#' @keywords internal



structureTestMGDDecisions <- function(.object){
  
  ## Install dplyr, tidyr, and rlang if not already installed
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop2(
      "Package `dplyr` required. Use `install.packages(\"dplyr\")` and rerun.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop2(
      "Package `tidyr` (> 1.0.0) required. Use `install.packages(\"tidyr\")` and rerun.")
  }
  if (!requireNamespace("rlang", quietly = TRUE)) {
    stop2(
      "Package `rlang` required. Use `install.packages(\"rlang\")` and rerun.")
  }
  
  # Check if of class cSEMTestMGD
  if(!inherits(.object, "cSEMTestMGD")){
    stop2("The following error occured in the `structureTestMGDDecisions()` function:\n",
          "Object must be of class `cSEMResults`")
  }
  
  # Helper function ------------------------------------------------------------
  
  # Change CI alpha equivlanet (1-CI)
  changeCItoAlphaEquivalent <- function(.names){
    # .names <- c("99%","95%","90%")
    tmpNames <- as.numeric(gsub("%", "", .names))
    tmpNames <- 100-tmpNames
    tmpNames <- paste0(tmpNames, "%")
    return(tmpNames)
  }
  
  out <- dplyr::tibble()
  # Klesel ---------------------------------------------------------------------
  
  # Test only provides an overall decision 
  
  if("Klesel" %in% names(.object)){
    klesel <- .object$Klesel$Decision %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "alpha") %>%
      # longer format
      tidyr::pivot_longer(cols = c(.data$dG, .data$dL), 
                          names_to = "Distance_metric", 
                          values_to = "Decision") %>%
      tidyr::pivot_wider(names_from = .data$alpha, values_from = .data$Decision) %>%
      dplyr::mutate(Test = "Klesel", Comparison = "overall")
    
    out <- dplyr::bind_rows(out, klesel)
  }
  
  # Chin  ----------------------------------------------------------------------
  
  if("Chin" %in% names(.object)){
    
    # overall decision
    chin1 <- .object$Chin$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "p-value_correction") %>%
      dplyr::mutate(Test = "Chin", Comparison = "overall")
    
    chin2 <-  .object$Chin$Decision %>% 
      purrr::modify_depth(3, dplyr::bind_rows) %>%
      purrr::modify_depth(2, dplyr::bind_rows, .id = "Comparison") %>%
      purrr::modify_depth(1, dplyr::bind_rows, .id = "alpha") %>%
      dplyr::bind_rows(.id = "p-value_correction") %>%
      # Check all comparisons
      dplyr::group_by(.data$`p-value_correction`, .data$alpha) %>%
      dplyr::summarise_at(dplyr::vars(-.data$Comparison), all) %>%
      dplyr::ungroup() %>% 
      # put results into format 
      tidyr::pivot_longer(cols = (-c(.data$`p-value_correction`, .data$alpha)),
                          names_to = "Comparison", 
                          values_to = "Decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = .data$alpha, values_from = .data$Decision) %>%
      dplyr::mutate(Test = "Chin")
    
    out <- dplyr::bind_rows(out, chin1, chin2)
  }
  
  # Sarstedt  ------------------------------------------------------------------
  
  if("Sarstedt" %in% names(.object)){
    
    # overall result
    sarstedt1 <- .object$Sarstedt$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "p-value_correction") %>%
      dplyr::mutate(Test = "Sarstedt", Comparison = "overall")
    
    # decision based on single paths
    sarstedt2 <- .object$Sarstedt$Decision %>%
      purrr::modify_depth(2, dplyr::bind_rows) %>%
      purrr::map(dplyr::bind_rows, .id = "alpha") %>%
      dplyr::bind_rows(.id = "p-value_correction") %>%
      # put results into format 
      tidyr::pivot_longer(cols = (-c(.data$`p-value_correction`, .data$alpha)),
                          names_to = "Comparison", 
                          values_to = "Decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = .data$alpha, values_from = .data$Decision) %>%
      dplyr::mutate(Test = "Sarstedt")
    
    out <- dplyr::bind_rows(out, sarstedt1, sarstedt2)
  }
  
  # Keil  ----------------------------------------------------------------------
  
  if("Keil" %in% names(.object)){
    
    # overall decision
    keil1 <- .object$Keil$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "p-value_correction") %>%
      dplyr::mutate(Test = "Keil", Comparison = "overall")
    
    # decision based on single paths
    keil2 <- .object$Keil$Decision %>%
      purrr::modify_depth(3, dplyr::bind_rows) %>%
      purrr::modify_depth(2, dplyr::bind_rows, .id = "Comparison") %>%
      purrr::modify_depth(1, dplyr::bind_rows, .id = "alpha") %>%
      dplyr::bind_rows(.id = "p-value_correction") %>%
      dplyr::group_by(.data$`p-value_correction`, .data$alpha) %>%
      dplyr::summarise_at(dplyr::vars(-.data$Comparison),all) %>%
      dplyr::ungroup() %>%
      # put results into format 
      tidyr::pivot_longer(cols = (-c(.data$`p-value_correction`, .data$alpha)),
                          names_to = "Comparison", 
                          values_to = "Decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = .data$alpha, values_from = .data$Decision) %>%
      dplyr::mutate(Test = "Keil")
    
    out <- dplyr::bind_rows(out, keil1, keil2)
  }
  
  # Nitzl  ---------------------------------------------------------------------
  
  if("Nitzl" %in% names(.object)){
    
    # overall decision
    nitzl1 <- .object$Nitzl$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "p-value_correction") %>%
      dplyr::mutate(Test = "Nitzl", Comparison = "overall")
    
    # decision based on single paths
    nitzl2 <- .object$Nitzl$Decision %>%
      purrr::modify_depth(3, dplyr::bind_rows) %>%
      purrr::modify_depth(2, dplyr::bind_rows, .id = "Comparison") %>%
      purrr::modify_depth(1, dplyr::bind_rows, .id = "alpha") %>%
      dplyr::bind_rows(.id = "p-value_correction") %>%
      dplyr::group_by(.data$`p-value_correction`, .data$alpha) %>%
      dplyr::summarise_at(dplyr::vars(-.data$Comparison),all) %>%
      dplyr::ungroup() %>%
      # put results into format 
      tidyr::pivot_longer(cols = (-c(.data$`p-value_correction`, .data$alpha)),
                          names_to = "Comparison", 
                          values_to = "Decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = .data$alpha, values_from = .data$Decision) %>%
      dplyr::mutate(Test = "Nitzl")
    
    out <- dplyr::bind_rows(out, nitzl1, nitzl2)
  }
  
  # Henseler  ------------------------------------------------------------------
  
  if("Henseler" %in% names(.object)){
    
    # overall decision
    henseler1 <- .object$Henseler$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "p-value_correction") %>%
      dplyr::mutate(Test = "Henseler", Comparison = "overall")
    
    # decision based on single paths
    henseler2 <- .object$Henseler$Decision %>%
      purrr::modify_depth(3, dplyr::bind_rows) %>%
      purrr::modify_depth(2, dplyr::bind_rows, .id = "Comparison") %>%
      purrr::modify_depth(1, dplyr::bind_rows, .id = "alpha") %>%
      dplyr::bind_rows(.id = "p-value_correction") %>%
      dplyr::group_by(.data$`p-value_correction`, .data$alpha) %>%
      dplyr::summarise_at(dplyr::vars(-.data$Comparison),all) %>%
      dplyr::ungroup() %>%
      # put results into format 
      tidyr::pivot_longer(cols = (-c(.data$`p-value_correction`, .data$alpha)),
                          names_to = "Comparison", 
                          values_to = "Decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = .data$alpha, values_from = .data$Decision) %>%
      dplyr::mutate(Test = "Henseler")
    
    out <- dplyr::bind_rows(out, henseler1, henseler2)
  }
  
  # CI para  -------------------------------------------------------------------
  
  if("CI_para" %in% names(.object)){
    
    # overall decision
    CIpara1 <- .object$CI_para$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "alpha") %>%
      tidyr::pivot_longer(cols = dplyr::contains("CI"), 
                          names_to = "CI_type", 
                          values_to = "Decision") %>%
      tidyr::pivot_wider(names_from = .data$alpha, values_from = .data$Decision) %>%
      # Change % as equivlanet to alpha
      dplyr::rename_at(dplyr::vars(dplyr::contains("%")), changeCItoAlphaEquivalent) %>%
      dplyr::mutate(Test = "CI_para", Comparison = "overall")
    
    # decision based on single paths
    CIpara2 <-  .object$CI_para$Decision %>%
      purrr::modify_depth(3, dplyr::bind_rows) %>%
      purrr::modify_depth(3, ~ dplyr::select(., .data$Name, .data$Decision)) %>%
      purrr::modify_depth(2, dplyr::bind_rows, .id = "CI_type") %>%
      purrr::modify_depth(1, dplyr::bind_rows, .id = "Comparison") %>%
      dplyr::bind_rows(.id = "alpha") %>%
      tidyr::pivot_wider(names_from = .data$Name, values_from = .data$Decision) %>%
      dplyr::group_by(.data$alpha, .data$CI_type) %>%
      dplyr::summarise_at(dplyr::vars(-.data$Comparison), all) %>%
      dplyr::ungroup() %>%
      # put results into format 
      tidyr::pivot_longer(cols = (-c(.data$CI_type, .data$alpha)),
                          names_to = "Comparison", 
                          values_to = "Decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = .data$alpha, values_from = .data$Decision) %>%
      # Change % as equivalent to alpha
      dplyr::rename_at(dplyr::vars(dplyr::contains("%")), changeCItoAlphaEquivalent) %>%
      dplyr::mutate(Test = "CI_para")
    
    out <- dplyr::bind_rows(out, CIpara1, CIpara2)
    
  }
  # CI overlap  ----------------------------------------------------------------
  
  if("CI_overlap" %in% names(.object)){
    
    # overall decision
    CIoverlap1 <- .object$CI_overlap$Decision_overall %>% 
      purrr::map(dplyr::bind_rows) %>% 
      dplyr::bind_rows(.id = "alpha") %>%
      tidyr::pivot_longer(cols = dplyr::contains("CI"), 
                          names_to = "CI_type", 
                          values_to = "Decision") %>%
      tidyr::pivot_wider(names_from = .data$alpha, values_from = .data$Decision) %>%
      # Change % as equivalent to alpha
      dplyr::rename_at(dplyr::vars(dplyr::contains("%")), changeCItoAlphaEquivalent) %>%
      dplyr::mutate(Test = "CI_overlap", Comparison = "overall")
    
    # decision based on single paths
    CIoverlap2 <-  .object$CI_overlap$Decision %>%
      purrr::modify_depth(3, dplyr::bind_rows) %>%
      purrr::modify_depth(3, ~ dplyr::select(., .data$Name, .data$Decision)) %>%
      purrr::modify_depth(2, dplyr::bind_rows, .id = "CI_type") %>%
      purrr::modify_depth(1, dplyr::bind_rows, .id = "Comparison") %>%
      dplyr::bind_rows(.id = "alpha") %>%
      tidyr::pivot_wider(names_from = .data$Name, values_from = .data$Decision) %>%
      dplyr::group_by(.data$alpha, .data$CI_type) %>%
      dplyr::summarise_at(dplyr::vars(-.data$Comparison), all) %>%
      dplyr::ungroup() %>%
      # put results into format 
      tidyr::pivot_longer(cols = (-c(.data$CI_type, .data$alpha)),
                          names_to = "Comparison", 
                          values_to = "Decision") %>%
      # put alphas in cols
      tidyr::pivot_wider(names_from = .data$alpha, values_from = .data$Decision) %>%
      # Change % as equivalent to alpha
      dplyr::rename_at(dplyr::vars(dplyr::contains("%")), changeCItoAlphaEquivalent) %>%
      dplyr::mutate(Test = "CI_overlap")
    
    out <- dplyr::bind_rows(out, CIoverlap1, CIoverlap2)
  }
  
  # Summarize all results ------------------------------------------------------
  
  out <- dplyr::select(out, .data$Test, .data$Comparison, dplyr::contains("%"), 
                       .data$`p-value_correction`, .data$CI_type, 
                       .data$Distance_metric)
  
  return(out)
}


