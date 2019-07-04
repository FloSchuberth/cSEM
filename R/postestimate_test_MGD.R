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
#'  .approach_alpha_adjust = args_default()$.approach_alpha_adjust,
#'  .approach_mgd          = args_default()$.approach_mgd,
#'  .model            = args_default()$.model,
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
  .approach_alpha_adjust = args_default()$.approach_alpha_adjust,
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
  if("Sarstedt" %in% .approach_mgd & .handle_inadmissibles == "drop"){
    stop2(
      "The following error occured in the testMGD() function:\n",
      "Approach `'Sarstedt'` not supported if `.handle_inadmissibles == 'drop'`")
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
  ## Get the model-implied VCV for Klesel et al. (2019)
  fit <- fit(.object = .object,
             .saturated = .saturated,
             .type_vcv = .type_vcv)
  
  ## Get bootstrapped parameter estimates for Sarstedt et al. (2011)
  if("Sarstedt" %in% .approach_mgd){
    
    # Check if .object already contains resamples; if not; run bootstrap
    if(!inherits(.object, "cSEMResults_resampled")) {
      .object <- resamplecSEMResults(
        .object               = .object,
        .resample_method      = "bootstrap",
        .handle_inadmissibles = .handle_inadmissibles,
        .R                    = .R_bootstrap) 
    }
    
    ## Combine bootstrap results in one matrix
    l <- lapply(.object, function(x) {
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
    
    l <- purrr::transpose(l)
    id_Sarstedt <- rep(1:length(.object), unlist(l$n))
    
    ll <- lapply(l, function(x) do.call(rbind, x))
    
    # Matrix that contains all parameter estimate that should compared 
    # plus an id variable indicating the group
    all_comb <- cbind(ll$path_resamples, ll$loading_resamples, ll$weight_resamples, id = id_Sarstedt)
  }
  
  ## Compute the test statistics 
  teststat <- list(
    #  Approach suggested by Klesel et al. (2019)
    "Klesel" = c(
      "dG" = calculateDistance(.matrices = fit, .distance = "geodesic"),
      "dL" = calculateDistance(.matrices = fit, .distance = "squared_euclidian")),
    # Approach suggested by Chin & Dibbern (2010)
    "Chin" = calculateParameterDifference(.object = .object, .model = .model)
    )
  
  # Approach suggested by Sarstedt et al. (2011) 
  if("Sarstedt" %in% .approach_mgd){
    
    ## Get the name of the parameters to be compared
    names_param <- getParameterNames(.object, .model = .model)
    
    # Calculate test statistic
    temp <- sapply(unlist(names_param), function(x) {
      calculateFR(
        .Parameter = all_comb[, x],
        .id        = all_comb[,ncol(all_comb)])
      })
    
    names(temp) <- unlist(names_param)
    ## Add test statistic
    teststat[["Sarstedt"]] <- temp
  } 
  
  ## Start Permutation
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
      
      # Calculate test statistic for permutation sample
      
      fit_temp <- fit(Est_temp, .saturated = .saturated, .type_vcv = .type_vcv)
      
      
      teststat_permutation <- list(
        "Klesel" = c(
          "dG" = calculateDistance(.matrices = fit_temp, .distance = "geodesic"),
          "dL" = calculateDistance(.matrices = fit_temp, .distance = "squared_euclidian")),
        "Chin" = calculateParameterDifference(.object=Est_temp,.model = .model))
      
      if("Sarstedt" %in% .approach_mgd){
        
        # Permutation of the bootstrap parameter estimates
        all_comb_permutation <- all_comb
        all_comb_permutation[,"id"] <- sample(id_Sarstedt)
        
        teststat_Sarstedt_permutation <- sapply(unlist(names_param),
                                        function(x){calculateFR(.Parameter = all_comb_permutation[,x],
                                                                .id = all_comb_permutation[,"id"])})
        names(teststat_Sarstedt_permutation) <- unlist(names_param)
        teststat_permutation[["Sarstedt"]] <- teststat_Sarstedt_permutation
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
      ## Stop if 10000 runs did not result in insufficient admissible results
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
  
  # Delete potential NA's
  ref_dist1 <- Filter(Negate(anyNA), ref_dist)
  
  # Order significance levels
  .alpha <- .alpha[order(.alpha)]
  
  ### Approach suggested by Klesel et al. (2019) -------------------------------
  ## Collect permuation results and combine
  ref_dist_Klesel <- lapply(ref_dist1, function(x) x$Klesel)
  ref_dist_matrix_Klesel <- do.call(cbind, ref_dist_Klesel)
  
  # Extract test statistic
  teststat_Klesel <- teststat$Klesel
  

  # Calculation of p-values
  pvalue_Klesel=rowMeans(ref_dist_matrix_Klesel> teststat_Klesel)
  
  # Decision 
  decision_Klesel <- lapply(.alpha,function(x){
    pvalue_Klesel > x
  })
  names(decision_Klesel) <- paste(.alpha*100,"%",sep= '')
  # TRUE = p-value > alpha --> not reject
  # FALSE = sufficient evidence against the H0 --> reject
 
  
  ### Approach suggested by Chin & Dibbern (2010) ------------------------------
  
  # Extract test statistic
  teststat_Chin <- lapply(teststat$Chin, function(x) x[!is.na(x)])
  
  # Create list with matrices containing the reference distribution of the parameter differences
  ref_dist_Chin <- lapply(ref_dist1,function(x) x$Chin)
  
  ref_dist_Chin_temp <- purrr::transpose(ref_dist_Chin)
  
  ref_dist_matrices_Chin <- lapply(ref_dist_Chin_temp, function(x) {
    temp <- do.call(cbind, x)
    temp_ind <- stats::complete.cases(temp)
    temp[temp_ind, ,drop = FALSE]
  })
  
  # Calculation of the p-values
  pvalue_Chin <- lapply(1:length(ref_dist_matrices_Chin),function(x){
    # share of values above the positive test statistic
    rowMeans(ref_dist_matrices_Chin[[x]]>abs(teststat_Chin[[x]]))+
    # share of values of the reference distribution below the negative test statistic 
      rowMeans(ref_dist_matrices_Chin[[x]]< (-abs(teststat_Chin[[x]])))
  })
  
  names(pvalue_Chin) <- names(ref_dist_matrices_Chin)
  
  padjusted_Chin <- lapply(.approach_alpha_adjust, function(x){
    # It is important to unlist the pvalues as pAdjust needs to now how many p-values
    # there are to do a proper adjustment
  pvector <- stats::p.adjust(unlist(pvalue_Chin),method = x)
    # Sort them back into list
    relist(flesh = pvector,skeleton = pvalue_Chin)
  })
  
  names(padjusted_Chin) <- .approach_alpha_adjust
  

  # Decision 
  decision_Chin=lapply(padjusted_Chin,function(adjust_approach){#over the different p adjustments
      temp=lapply(.alpha, function(alpha){#over the different significance levels
        lapply(adjust_approach,function(group_comp){#over the different group comparisons
          # check whether the p values are larger than a certain alpha
          group_comp > alpha
        })
      })
  names(temp) = paste(.alpha*100,"%",sep= '')
  temp
      })


  # Overall decision, i.e., was any of the test belonging to one significance level 
  # and one p value adjustment rejected
  decision_overall_Chin = lapply(decision_Chin, function(decision_Chin_list){
    lapply(decision_Chin_list,function(x){
      all(unlist(x))
    })
  })
  
  # Approach suggested by Sarstedt et al. (2011)-------------------------------------------
  if("Sarstedt" %in% .approach_mgd){

    # Extract test statistic
    teststat_Sarstedt <- teststat$Sarstedt
    
    # Collect permuation results and combine to matrix
    ref_dist_Sarstedt <- lapply(ref_dist1, function(x) x$Sarstedt)
    ref_dist_matrix_Sarstedt <- do.call(cbind, ref_dist_Sarstedt)
    
    # Calculation of the p-value
    pvalue_Sarstedt=rowMeans(ref_dist_matrix_Sarstedt> teststat_Sarstedt)
    
    # Adjust pvalues:
    padjusted_Sarstedt<- lapply(.approach_alpha_adjust, function(x){
      pvector <- stats::p.adjust(pvalue_Sarstedt,method = x)
    })
    names(padjusted_Sarstedt) <- .approach_alpha_adjust
    
    # Decision 
    decision_Sarstedt=lapply(padjusted_Sarstedt,function(p_value){
      temp=lapply(.alpha, function(alpha){
        p_value > alpha
      })
      names(temp) = paste(.alpha*100,"%",sep= '')
      temp
    })
    # TRUE -> no rejection 
    
    # Decision overall
    decision_overall_Sarstedt <- lapply(decision_Sarstedt, function(x){#over p-value adjustments
      lapply(x, function(xx){ #over different significant levels
        all(xx)
        })
    })
    
  }#End approach Sarstedt
  
  ### Return output ------------------------------------------------------------
  out <- list(
    "Klesel"=list(
      "Test_statistic"     = teststat_Klesel,
      "Critical_value"     = critical_values_Klesel, 
      "P_value"            = pvalue_Klesel,
      "Decision"           = decision_Klesel), 
    
    "Chin" = list(
      "Test_statistic"     = teststat_Chin,
      "Critical_value"     = critical_values_Chin, 
      "P_value"            = pvalue_Chin,
      "P_value_adjusted"   = padjusted_Chin,
      "Decision"           = decision_Chin,
      "Decision_overall"   = decision_overall_Chin
      ),
    "Information"        = list(
      "Number_admissibles"    = ncol(ref_dist_matrix_Klesel),
      "Total_runs"            = counter + n_inadmissibles,
      "Group_names"           = names(.object),
      "Number_of_observations"= sapply(X_all_list, nrow),
      "Permutation_values"      = list(
        "Klesel" = ref_dist_matrix_Klesel,
        "Chin"   = ref_dist_matrices_Chin),
      "Approach" = .approach_mgd,
      "Seed"     = .seed,
      "Alpha"    = .alpha
    )
  )
  
  if("Sarstedt" %in% .approach_mgd){
    out[["Sarstedt"]] <- list(
      "Test_statistic"   = teststat_Sarstedt,
      "Critical_value"     = critical_values_Sarstedt,
      "P_value"            = pvalue_Sarstedt,
      "P_value_adjusted"   = padjusted_Sarstedt,
      "Decision"           = decision_Sarstedt,
      "Decision_overall"   = decision_overall_Sarstedt
    )
    
    # Add the reference distirbutions
    out[["Information"]][["Permutation_values"]][["Sarstedt"]] <- ref_dist_matrix_Sarstedt
    
    # Order output
    out <- out[c("Klesel","Chin","Sarstedt","Information")]
  }

  
  ## Remove the seed since it is set globally. Reset immediately by calling
  ## any kind of function that requires .Random.seed as this causes R to
  ## to create a new one.
  rm(.Random.seed, envir = .GlobalEnv)
  runif(1) # dont remove; this sets up a new .Random.seed
  
  class(out) <- "cSEMTestMGD"
  return(out)
}
