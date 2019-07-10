#' Tests for multiple groups. 
#'
#' This function performs several permutation tests, i.e., the reference distribution 
#' of the test statistic is obtained by permutation.
#' 
#' The following test are implemented:
#' \describe{
#' \item{Approach suggested by \insertCite{Klesel2019;textual}{cSEM}}{
#'   The model-implied variance-covariance matrix (either indicator or construct) is
#'   compared across groups. 
#' 
#'   To measure the distance between the model-implied variance-covariance matrices, 
#'   the geodesic distance (dG) and the squared Euclidean distance (dL) are used.}
#' \item{Approach suggested by \insertCite{Chin2010;textual}{cSEM}}{
#'   Groups are compared in terms of parameter differences across groups.
#' }
#' }
#' 
#' @usage testMGD(
#'  .object                = args_default()$.object,
#'  .alpha                 = args_default()$.alpha,
#'  .approach_p_adjust     = args_default()$.approach_p_adjust,
#'  .approach_mgd          = args_default()$.approach_mgd,
#'  .model                 = args_default()$.model,
#'  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
#'  .R_permutation         = args_default()$.R_permutation,
#'  .R_bootstrap           = args_default()$.R_bootstrap,
#'  .saturated             = args_default()$.saturated,
#'  .seed                  = args_default()$.seed,
#'  .type_vcv              = args_default()$.type_vcv,
#'  .verbose               = args_default()$.verbose
#'  ) 
#' 
#' @inheritParams csem_arguments
#' @param .model A model in [lavaan model syntax][lavaan::model.syntax] indicating which 
#'   parameters (i.e, path (`~`), loadings (`=~`), or weights (`<~`)) should be
#'   compared across groups. Defaults to `NULL` in which case all parameters of the model
#'   are compared.
#'   
#' @inherit csem_test return
#'
#' @references
#'   \insertAllCited{}
#'   
#' @seealso [cSEMResults]
#'
#' @examples
#' \dontrun{
#' require(cSEM)
#' data(satisfaction)
#'
#' model <- "
#' # Structural model
#' QUAL ~ EXPE
#' EXPE ~ IMAG
#' SAT  ~ IMAG + EXPE + QUAL + VAL
#' LOY  ~ IMAG + SAT
#' VAL  ~ EXPE + QUAL
#'
#' # Measurement model
#'
#' EXPE <~ expe1 + expe2 + expe3 + expe4 + expe5
#' IMAG <~ imag1 + imag2 + imag3 + imag4 + imag5
#' LOY  =~ loy1  + loy2  + loy3  + loy4
#' QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
#' SAT  <~ sat1  + sat2  + sat3  + sat4
#' VAL  <~ val1  + val2  + val3  + val4
#' "
#' 
#' listData <- list(satisfaction[-3,], satisfaction[-5, ], satisfaction[-10, ])
#' out.cSEM <- csem(listData, model) 
#'
#' testMGD(.object = out.cSEM, .R = 20, .type_vcv= 'construct')
#' }
#'
#' @export

testMGD <- function(
  .object                = args_default()$.object,
  .alpha                 = args_default()$.alpha,
  .approach_p_adjust     = args_default()$.approach_p_adjust,
  .approach_mgd          = args_default()$.approach_mgd,
  .model                 = args_default()$.model,
  .handle_inadmissibles  = args_default()$.handle_inadmissibles,
  .R_permutation         = args_default()$.R_permutation,
  .R_bootstrap           = args_default()$.R_bootstrap,
  .saturated             = args_default()$.saturated,
  .seed                  = args_default()$.seed,
  .type_vcv              = args_default()$.type_vcv,
  .verbose               = args_default()$.verbose
){

  ## Match arguments
  match.arg(.handle_inadmissibles, args_default(.choices = TRUE)$.handle_inadmissibles)
  match.arg(.type_vcv, args_default(.choices = TRUE)$.type_vcv)
  
  ### Checks and errors ========================================================
  ## Check if at least two groups are present
  if(!inherits(.object, "cSEMResults_multi")) {
    stop2(
      "The following error occured in the testMGD() function:\n",
      "At least two groups required."
    )
  } 
  
  # Sarstedt et al. (2011) is not allowed to be used in combination with 
  # .handle_inadmissibles == "drop as permutation test statistic are be dropped
  # only because of estimations based on the permutated dataset are dropped.
  if(any(.approach_mgd %in% c("all", "Sarstedt")) & .handle_inadmissibles == "drop"){
    stop2(
      "The following error occured in the testMGD() function:\n",
      "Approach `'Sarstedt'` not supported if `.handle_inadmissibles == 'drop'`")
  }
  
  # If a non-linear model is used Klesel et al. approach cannot be used as 
  # it is not clear how to calculate the model-implied VCV
  # Extract model type information
  if(inherits(.object[[1]], "cSEMResults_2ndorder")){
    model_type <- .object[[1]]$Second_stage$Information$Model$model_type
  } else {
    model_type <- .object[[1]]$Information$Model$model_type
  }
  
  if(any(.approach_mgd %in% c("all", "Klesel")) & model_type == "Nonlinear"){
    stop2("The following error occured in the testMGD() function:\n",
          "The approach suggested by Klesel et al. (2019) cannot be applied",
          " to nonlinear models as cSEM currently cannot calculate",
          " the model-implied VCV matrix for such models.\n", 
          "Consider setting `.approach_mgd = c('Chin',  'Sarstedt')`")
  }
  
  ## Check if data for different groups is identical
  if(.verbose) {
    ## Check if any of the group estimates are inadmissible
    if(sum(unlist(verify(.object))) != 0) {
      warning2(
        "The following warning occured in the testMGD() function:\n",
        "Initial estimation results for at least one group are inadmissible.\n", 
        "See `verify(.object)` for details.")
    }
    
    if(TRUE %in% lapply(utils::combn(.object, 2, simplify = FALSE),
                        function(x){ identical(x[[1]], x[[2]])})){
      warning2(
        "The following warning occured in the testMGD() function:\n",
        "At least two groups are identical. Results may not be meaningful.")
    } 
  }
  
  ### Calculation of the test statistics========================================
  teststat <- list()
  
  ## Klesel et al. (2019) ------------------------------------------------------
  if(any(.approach_mgd %in% c("all", "Klesel"))) {
    ## Get the model-implied VCV
    fit <- fit(.object    = .object,
               .saturated = .saturated,
               .type_vcv  = .type_vcv)
    
    ## Compute test statistic
    temp <- c(
      "dG" = calculateDistance(.matrices = fit, .distance = "geodesic"),
      "dL" = calculateDistance(.matrices = fit, .distance = "squared_euclidian")
    )

    ## Save test statistic
    teststat[["Klesel"]] <- temp
  }
  
  ## Chin & Dibbern (2010) -----------------------------------------------------
  if(any(.approach_mgd %in% c("all", "Chin"))) {
    ## Compute and save test statistic
    teststat[["Chin"]] <- calculateParameterDifference(.object = .object, 
                                                       .model = .model)
  }
  
  ## Sarstedt et al. (2011) ----------------------------------------------------
  if(any(.approach_mgd %in% c("all", "Sarstedt"))) {
    
    ## Check if .object already contains resamples; if not, run bootstrap
    if(!inherits(.object, "cSEMResults_resampled")) {
      .object <- resamplecSEMResults(
        .object               = .object,
        .resample_method      = "bootstrap",
        .handle_inadmissibles = .handle_inadmissibles,
        .R                    = .R_bootstrap,
        .seed                 = .seed) 
    }
    
    ## Combine bootstrap results in one matrix
    ll <- lapply(.object, function(x) {
      if(inherits(.object, "cSEMResults_2ndorder")) {
        x <- x$Second_stage$Information$Resamples$Estimates$Estimates1
      } else {
        x <- x$Estimates$Estimates_resample$Estimates1
      }
      path_resamples    <- x$Path_estimates$Resampled
      loading_resamples <- x$Loading_estimates$Resampled
      weight_resamples  <- x$Weight_estimates$Resampled
      n                 <- nrow(path_resamples)
      
      list(
        "path_resamples"    = path_resamples, 
        "loading_resamples" = loading_resamples, 
        "weight_resamples"  = weight_resamples, 
        "n"                 = n)
    })
    
    ## Transpose, get id column, bind rows and columns
    ll <- purrr::transpose(ll)
    group_id <- rep(1:length(.object), unlist(ll$n))
    ll <- lapply(ll, function(x) do.call(rbind, x))
    
    all_comb <- cbind(ll$path_resamples, 
                      ll$loading_resamples, 
                      ll$weight_resamples, 
                      "group_id" = group_id)
    # all_comb contains all parameter estimate that could potentially be compared 
    # plus an id column indicating the group adherance of each row.
    
    ## Get the name of the parameters to be compared
    names_param <- unlist(getParameterNames(.object, .model = .model))
    
    ## Select relevant columns
    all_comb <- all_comb[, c(names_param, "group_id")]
    
    ## Add test statistic
    teststat[["Sarstedt"]] <- calculateFR(.resample_sarstedt = all_comb)
  }
  
  ### Permutation ==============================================================
  ## Preparation
  # Put data of each groups in a list and combine
  if(inherits(.object, "cSEMResults_2ndorder")) {
    # Data is saved in the first stage
    X_all_list  <- lapply(.object, function(x) x$First_stage$Information$Data)
    # Collect initial arguments (from the first object, but could be any other)
    arguments <- .object[[1]]$Second_stage$Information$Arguments_original
  } else {
    X_all_list  <- lapply(.object, function(x) x$Information$Data)
    # Collect initial arguments (from the first object, but could be any other)
    arguments <- .object[[1]]$Information$Arguments
  }
  X_all <- do.call(rbind, X_all_list)
  
  # Create a vector "id" to be used to randomly select groups (permutate) and
  # set id as an argument in order to identify the groups.
  id <- rep(1:length(X_all_list), sapply(X_all_list, nrow))
  arguments[[".id"]] <- "id"
  
  # Start progress bar if required
  if(.verbose){
    pb <- txtProgressBar(min = 0, max = .R_permutation, style = 3)
  }
  
  ## Create seed if not already set
  if(is.null(.seed)) {
    .seed <- sample(.Random.seed, 1)
  }
  ## Set seed
  set.seed(.seed)
  
  ## Calculate reference distribution
  ref_dist        <- list()
  n_inadmissibles  <- 0
  counter <- 0
  repeat{
    # Counter
    counter <- counter + 1
    
    # Permutate data
    X_temp <- cbind(X_all, id = sample(id))
    
    # Replace the old dataset by the new permutated dataset
    arguments[[".data"]] <- X_temp
    
    # Estimate model
    Est_temp <- do.call(csem, arguments)   
    
    # Check status
    status_code <- sum(unlist(verify(Est_temp)))
    
    # Distinguish depending on how inadmissibles should be handled
    if(status_code == 0 | (status_code != 0 & .handle_inadmissibles == "ignore")) {
      # Compute if status is ok or .handle inadmissibles = "ignore" AND the status is 
      # not ok
      
      ### Calculation of the test statistic for each resample ==================
      teststat_permutation <- list()
      
      ## Klesel et al. (2019) --------------------------------------------------
      if(any(.approach_mgd %in% c("all", "Klesel"))) {
        ## Get the model-implied VCV
        fit_temp <- fit(Est_temp, .saturated = .saturated, .type_vcv = .type_vcv)
        
        ## Compute test statistic
        temp <- c(
          "dG" = calculateDistance(.matrices = fit_temp, .distance = "geodesic"),
          "dL" = calculateDistance(.matrices = fit_temp, .distance = "squared_euclidian")
        )
        
        ## Save test statistic
        teststat_permutation[["Klesel"]] <- temp
      }
      
      ## Chin & Dibbern (2010) -------------------------------------------------
      if(any(.approach_mgd %in% c("all", "Chin"))) {
        ## Compute and save test statistic
        teststat_permutation[["Chin"]] <- calculateParameterDifference(
          .object = Est_temp, 
          .model  = .model)
      }
      ## Sarstedt et al. (2011) ------------------------------------------------
      if(any(.approach_mgd %in% c("all", "Sarstedt"))) {
        
        # Permutation of the bootstrap parameter estimates
        all_comb_permutation <- all_comb
        all_comb_permutation[ , "group_id"] <- sample(group_id)

        teststat_permutation[["Sarstedt"]] <- calculateFR(all_comb_permutation)
      }
      
      ref_dist[[counter]] <- teststat_permutation
      
    } else if(status_code != 0 & .handle_inadmissibles == "drop") {
      # Set list element to zero if status is not okay and .handle_inadmissibles == "drop"
      ref_dist[[counter]] <- NA
      
    } else {# status is not ok and .handle_inadmissibles == "replace"
      # Reset counter and raise number of inadmissibles by 1
      counter <- counter - 1
      n_inadmissibles <- n_inadmissibles + 1
    }
    
    # Break repeat loop if .R results have been created.
    if(length(ref_dist) == .R_permutation) {
      break
    } else if(counter + n_inadmissibles == 10000) { 
      # Stop if 10000 runs did not result in insufficient admissible results
      stop("Not enough admissible result.", call. = FALSE)
    }
    
    if(.verbose){
      setTxtProgressBar(pb, counter)
    }
    
  } # END repeat 
  
  # close progress bar
  if(.verbose){
    close(pb)
  }
  
  ### Postprocessing ===========================================================
  # Delete potential NA's
  ref_dist1 <- Filter(Negate(anyNA), ref_dist)
  
  # Order significance levels
  .alpha <- .alpha[order(.alpha)]
  
  ## Klesel et al. (2019) ------------------------------------------------------
  if(any(.approach_mgd %in% c("all", "Klesel"))) {
    
  # Collect permuation results and combine
  ref_dist_Klesel <- lapply(ref_dist1, function(x) x$Klesel)
  ref_dist_matrix_Klesel <- do.call(cbind, ref_dist_Klesel)
  
  # Extract test statistic
  teststat_Klesel <- teststat$Klesel
  
  # Calculation of p-values
  pvalue_Klesel <- rowMeans(ref_dist_matrix_Klesel > teststat_Klesel)
  
  # Decision 
  decision_Klesel <- lapply(.alpha, function(x) {
    pvalue_Klesel > x
  })
  
  names(decision_Klesel) <- paste0(.alpha * 100, "%")
  # TRUE = p-value > alpha --> not reject
  # FALSE = sufficient evidence against the H0 --> reject
  }
  
  ## Chin & Dibbern (2010) -----------------------------------------------------
  if(any(.approach_mgd %in% c("all", "Chin"))) {
    
    # Extract test statistic
    teststat_Chin <- teststat$Chin
    
    # Create list with matrices containing the reference distribution 
    # of the parameter differences
    ref_dist_Chin <- lapply(ref_dist1, function(x) x$Chin)
    
    # Transpose
    ref_dist_Chin_temp <- purrr::transpose(ref_dist_Chin)
    
    ref_dist_matrices_Chin <- lapply(ref_dist_Chin_temp, function(x) {
      temp <- do.call(cbind, x)
      temp_ind <- stats::complete.cases(temp)
      temp[temp_ind, ,drop = FALSE]
    })
    
    # Calculation of the p-values
    pvalue_Chin <- lapply(1:length(ref_dist_matrices_Chin), function(x) {
      # Share of values above the positive test statistic
      rowMeans(ref_dist_matrices_Chin[[x]] > abs(teststat_Chin[[x]])) +
        # share of values of the reference distribution below the negative test statistic 
        rowMeans(ref_dist_matrices_Chin[[x]] < (-abs(teststat_Chin[[x]])))
    })
    
    names(pvalue_Chin) <- names(ref_dist_matrices_Chin)
    
    # Adjust p-values
    padjusted_Chin <- lapply(as.list(.approach_p_adjust), function(x){
      # It is important to unlist the pvalues as pAdjust needs to now how many p-values
      # there are to do a proper adjustment
      pvector <- stats::p.adjust(unlist(pvalue_Chin),method = x)
      # Sort them back into list
      relist(flesh = pvector,skeleton = pvalue_Chin)
    })
    
    names(padjusted_Chin) <- .approach_p_adjust
    
    # Decision 
    decision_Chin <- lapply(padjusted_Chin, function(adjust_approach){ # over the different p adjustments
      temp <- lapply(.alpha, function(alpha){# over the different significance levels
        lapply(adjust_approach,function(group_comp){# over the different group comparisons
          # check whether the p values are larger than a certain alpha
          group_comp > alpha
        })
      })
      names(temp) <- paste0(.alpha*100, "%")
      temp
    })
    
    # Overall decision, i.e., was any of the test belonging to one significance level 
    # and one p value adjustment rejected
    decision_overall_Chin <- lapply(decision_Chin, function(decision_Chin_list){
      lapply(decision_Chin_list,function(x){
        all(unlist(x))
      })
    })
  }
  
  ## Sarstedt et al. (2011) ----------------------------------------------------
  if(any(.approach_mgd %in% c("all", "Sarstedt"))) {

    # Extract test statistic
    teststat_Sarstedt <- teststat$Sarstedt
    
    # Collect permuation results and combine to matrix
    ref_dist_Sarstedt <- lapply(ref_dist1, function(x) x$Sarstedt)
    ref_dist_matrix_Sarstedt <- do.call(cbind, ref_dist_Sarstedt)
    
    # Calculation of the p-value
    pvalue_Sarstedt <- rowMeans(ref_dist_matrix_Sarstedt > teststat_Sarstedt)
    
    # Adjust pvalues:
    padjusted_Sarstedt<- lapply(as.list(.approach_p_adjust), function(x){
      pvector <- stats::p.adjust(pvalue_Sarstedt, method = x)
    })
    names(padjusted_Sarstedt) <- .approach_p_adjust
    
    # Decision 
    decision_Sarstedt=lapply(padjusted_Sarstedt,function(p_value){
      temp=lapply(.alpha, function(alpha){
        p_value > alpha
      })
      names(temp) = paste(.alpha*100,"%",sep= '')
      temp
    })
    
    # Decision overall
    decision_overall_Sarstedt <- lapply(decision_Sarstedt, function(x){#over p-value adjustments
      lapply(x, function(xx){ #over different significant levels
        all(xx)
        })
    })
  }
  
  ### Return output ============================================================
  out <- list()
  
  ## Information
  out[["Information"]] <- list(
    "Number_admissibles"    = length(ref_dist1),
    "Total_runs"            = counter + n_inadmissibles,
    "Group_names"           = names(.object),
    "Number_of_observations"= sapply(X_all_list, nrow),
    "Approach"              = .approach_mgd,
    "Approach_p_adjust"     = .approach_p_adjust,
    "Seed"                  = .seed,
    "Alpha"                 = .alpha,
    "Permutation_values"    = list()
  )
  
  if(any(.approach_mgd %in% c("all", "Klesel"))) {
    out[["Klesel"]] <- list(
      "Test_statistic"     = teststat_Klesel,
      "P_value"            = pvalue_Klesel,
      "Decision"           = decision_Klesel
    )
    
    out[["Information"]][["Permutation_values"]][["Klesel"]] <- ref_dist_matrix_Klesel
  }
  
  if(any(.approach_mgd %in% c("all", "Chin"))) {
    out[["Chin"]] <- list(
      "Test_statistic"     = teststat_Chin,
      "P_value"            = pvalue_Chin,
      "P_value_adjusted"   = padjusted_Chin,
      "Decision"           = decision_Chin,
      "Decision_overall"   = decision_overall_Chin
    )
    
    out[["Information"]][["Permutation_values"]][["Chin"]] <- ref_dist_matrices_Chin
  }
  
  if(any(.approach_mgd %in% c("all", "Sarstedt"))) {
    out[["Sarstedt"]] <- list(
      "Test_statistic"     = teststat_Sarstedt,
      "P_value"            = pvalue_Sarstedt,
      "P_value_adjusted"   = padjusted_Sarstedt,
      "Decision"           = decision_Sarstedt,
      "Decision_overall"   = decision_overall_Sarstedt
    )
    
    out[["Information"]][["Permutation_values"]][["Sarstedt"]] <- ref_dist_matrix_Sarstedt
  }
  
  ## Remove the seed since it is set globally. Reset immediately by calling
  ## any kind of function that requires .Random.seed as this causes R to
  ## to create a new one.
  rm(.Random.seed, envir = .GlobalEnv)
  runif(1) # dont remove; this sets up a new .Random.seed
  
  class(out) <- "cSEMTestMGD"
  return(out)
}
