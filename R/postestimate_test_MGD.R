#' Tests for multi-group comparisons
#'
#' This function performs several permutation tests, i.e., the reference distribution 
#' of the test statistic is obtained by permutation.
#' 
#' The following tests are implemented:
#' \describe{
#' \item{Approach suggested by \insertCite{Klesel2019;textual}{cSEM}}{
#'   The model-implied variance-covariance matrix (either indicator 
#'   `.type_vcv = "indicator"` or construct `.type_vcv = "construct"`) 
#'   is compared across groups. 
#' 
#'   To measure the distance between the model-implied variance-covariance matrices, 
#'   the geodesic distance (dG) and the squared Euclidean distance (dL) are used.
#'   If more than two groups are compared, the average distance over all groups
#'   is used.}
#' \item{Approach suggested by \insertCite{Sarstedt2011;textual}{cSEM}}{
#'   Groups are compared in terms of parameter differences across groups.
#'   \insertCite{Sarstedt2011;textual}{cSEM} tests if parameter k is equal
#'   across all groups. If several parameters are tested simulaneously 
#'   it is recommended to adjust the signficance level or the p-values (in \pkg{cSEM} correction is
#'   done by p-value). By default
#'   no multiple testing correction is done, however, several common
#'   adjustments are available via `.approach_p_adjust`. See 
#'   \code{\link[stats:p.adjust]{stats::p.adjust()}} for details.
#' }
#' \item{Approach suggested by \insertCite{Chin2010;textual}{cSEM}}{
#'   Groups are compared in terms of parameter differences across groups.
#'   \insertCite{Chin2010;textual}{cSEM} tests if parameter k is equal
#'   between two groups. If more than two groups are tested the parameter is compared 
#'   between all pairs of groups. In this case, it is recommended
#'   to adjust the signficance level or the p-values (in \pkg{cSEM} correction is
#'   done by p-value). If several parameters are tested simultaneously, correction
#'   is by group and number of parameters. By default
#'   no multiple testing correction is done, however, several common
#'   adjustments are available via `.approach_p_adjust`. See 
#'   \code{\link[stats:p.adjust]{stats::p.adjust()}} for details.
#' }
#' \item{Approach suggested by \insertCite{Keil2000;textual}{cSEM}}{
#'   Groups are compared in terms of parameter differences across groups.
#'   \insertCite{Keil2000;textual}{cSEM} tests if parameter k is equal
#'   between two groups. It is assumed, that the standard errors of the coefficients are 
#'   equal across groups. The calculation of the standard error of the parameter 
#'   difference is adjusted as proposed by \insertCite{Henseler2009;textual}{cSEM}.
#'   If more than two groups are tested the parameter is compared 
#'   between all pairs of groups. In this case, it is recommended
#'   to adjust the signficance level or the p-values (in \pkg{cSEM} correction is
#'   done by p-value). If several parameters are tested simultaneously, correction
#'   is by group and number of parameters. By default
#'   no multiple testing correction is done, however, several common
#'   adjustments are available via `.approach_p_adjust`. See 
#'   \code{\link[stats:p.adjust]{stats::p.adjust()}} for details.
#' }
#' \item{Approach suggested by \insertCite{Nitzl2010;textual}{cSEM}}{
#'   Groups are compared in terms of parameter differences across groups.
#'   Similarly to \insertCite{Keil2000;textual}{cSEM}, a single parameter k is tested
#'   whether it is equal between two groups. In contrast to \insertCite{Keil2000;textual}{cSEM},
#'   it is assumed, that the standard errors of the coefficients are 
#'   equal across groups \insertCite{Sarstedt2011}{cSEM}.
#'   If more than two groups are tested the parameter is compared 
#'   between all pairs of groups. In this case, it is recommended
#'   to adjust the signficance level or the p-values (in \pkg{cSEM} correction is
#'   done by p-value). If several parameters are tested simultaneously, correction
#'   is by group and number of parameters. By default
#'   no multiple testing correction is done, however, several common
#'   adjustments are available via `.approach_p_adjust`. See 
#'   \code{\link[stats:p.adjust]{stats::p.adjust()}} for details.
#' }
#' \item{Approach suggested by \insertCite{Henseler2009;textual}{cSEM}}{
#'   This approach is also known as PLS-MGA. It compares \insertCite{Henseler2009,Sarstedt2011}{cSEM}
#'   Groups are compared in terms of parameter differences across groups.
#'   A single parameter k is tested whether it is equal between two groups. 
#'   In this case, it is recommended
#'   to adjust the signficance level or the p-values (in \pkg{cSEM} correction is
#'   done by p-value). If several parameters are tested simultaneously, correction
#'   is by group and number of parameters. By default
#'   no multiple testing correction is done, however, several common
#'   adjustments are available via `.approach_p_adjust`. See 
#'   \code{\link[stats:p.adjust]{stats::p.adjust()}} for details.
#' }
#' }
#' 
#' Use `.approach_mgd` to choose the approach. By default all approaches are computed
#' (`.approach_mgd = "all"`).
#' 
#' By default, approaches based on parameter differences across groups compare
#' all parameters (`.parameters_to_compare = NULL`). To compare only
#' a subset of parameters provide the parameters in lavaan model syntax just like
#' providing the model. Note that the "model" provided to `.parameters_to_compare`
#' does not have to be an estimatable model! See the example below.
#' 
#' Note that compared to all other functions in \pkg{cSEM}, `.handle_inadmissibles`
#' defaults to `"replace"` to accomdate the Sarstedt et al. (2011) approach.
#' 
#' 
#' @usage testMGD(
#'  .object                = NULL,
#'  .alpha                 = 0.05,
#'  .approach_p_adjust     = "none",
#'  .approach_mgd          = c("all", "Klesel", "Chin", "Sarstedt", "Keil", "Nitzl", "Henseler"),
#'  .parameters_to_compare = NULL,
#'  .handle_inadmissibles  = c("replace", "drop", "none"),
#'  .R_permutation         = 499,
#'  .R_bootstrap           = 499,
#'  .saturated             = FALSE,
#'  .seed                  = NULL,
#'  .type_vcv              = c("indicator", "construct"),
#'  .verbose               = TRUE
#'  ) 
#' 
#' @inheritParams csem_arguments
#' @param .handle_inadmissibles Character string. How should inadmissible results 
#'   be treated? One of "*drop*", "*ignore*", or "*replace*". If "*drop*", all
#'   replications/resamples yielding an inadmissible result will be dropped 
#'   (i.e. the number of results returned will potentially be less than `.R`). 
#'   For "*ignore*" all results are returned even if all or some of the replications
#'   yieled inadmissible results (i.e. number of results returned is equal to `.R`). 
#'   For "*replace*" resampling continues until there are exactly `.R` admissible solutions. 
#'   Defaults to "*replace*" to accomodate all approaches.
#'   
#' @return A list of class `cSEMTestMGD`. Technically, `cSEMTestMGD` is a 
#'   named list containing the following list elements:
#'
#' \describe{
#'   \item{`$Information`}{Additional information.}
#'   \item{`$Klesel`}{A list with elements, `Test_statistic`, `P_value`, and `Decision`}
#'   \item{`$Chin`}{A list with elements, `Test_statistic`, `P_value`, `Decision`, and `Decision_overall`}
#'   \item{`$Sarstedt`}{A list with elements, `Test_statistic`, `P_value`, `Decision`, and `Decision_overall`}
#'   \item{`$Keil`}{A list with elements, `Test_statistic`, `P_value`, `Decision`, and `Decision_overall`}
#'   \item{`$Nitzl`}{A list with elements, `Test_statistic`, `P_value`, `Decision`, and `Decision_overall`}
#' }
#' @references
#'   \insertAllCited{}
#'   
#' @seealso [cSEMResults]
#'
#' @example inst/examples/example_testMGD.R
#' 
#' @export

testMGD <- function(
 .object                = NULL,
 .alpha                 = 0.05,
 .approach_p_adjust     = "none",
 .approach_mgd          = c("all", "Klesel", "Chin", "Sarstedt", "Keil", "Nitzl"),
 .parameters_to_compare = NULL,
 .handle_inadmissibles  = c("replace", "drop", "ignore"),
 .R_permutation         = 499,
 .R_bootstrap           = 499,
 .saturated             = FALSE,
 .seed                  = NULL,
 .type_vcv              = c("indicator", "construct"),
 .verbose               = TRUE
){

  # Check .approach_mgd argument choices
  diff <- setdiff(.approach_mgd, args_default(TRUE)$.approach_mgd)

  if(length(diff) != 0) {
    stop2(
      "The following error occured in the testMGD() function:\n",
      "Unknown approach: ", paste0(diff, collapse = ", "), ".",
      " Possible choices are: ",
      paste0(args_default(TRUE)$.approach_mgd, collapse = ", "))
  }
  
  ## Match arguments
  .approach_mgd         <- match.arg(.approach_mgd, several.ok = TRUE)
  .handle_inadmissibles <- match.arg(.handle_inadmissibles) # no reference to 
                           # args_default, because args_default()$.handle_inadmissibles
                           # has "drop" as default, but testMGD hast "replace".
  .type_vcv             <- match.arg(.type_vcv)
  
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
  
  # If a nonlinear model is used Klesel et al. approach cannot be used as 
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
  
  ## Check if any of the group estimates are inadmissible
  if(sum(unlist(verify(.object))) != 0) {
    warning2(
      "The following warning occured in the testMGD() function:\n",
      "Initial estimation results for at least one group are inadmissible.\n", 
      "See `verify(.object)` for details.")
  }
  
  ## Check if data for different groups is identical
  if(TRUE %in% lapply(utils::combn(.object, 2, simplify = FALSE),
                      function(x){ identical(x[[1]], x[[2]])})){
    warning2(
      "The following warning occured in the testMGD() function:\n",
      "At least two groups are identical. Results may not be meaningful.")
  } 
  
  if(.verbose) {
    cat(rule2("Several tests for multi-group comparisons",
              type = 3), "\n\n")
  }
  
  
  ## Get the name of the parameters to be compared. This is required for Sarstedt and Hensler
  names_param <- unlist(getParameterNames(.object, .model = .parameters_to_compare))
  
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
                                                       .model = .parameters_to_compare)
  }
  
  ## Sarstedt et al. (2011), Keil et al. (2000), Nitzl (2010) ------------------
  if(any(.approach_mgd %in% c("all", "Sarstedt", "Keil", "Nitzl", "Henseler"))) {
    
    ## Check if .object already contains resamples; if not, run bootstrap
    if(!inherits(.object, "cSEMResults_resampled")) {
      if(.verbose) {
        cat("Bootstrapping cSEMResults object ...\n\n")
      }
      .object <- resamplecSEMResults(
        .object               = .object,
        .resample_method      = "bootstrap",
        .handle_inadmissibles = .handle_inadmissibles,
        .R                    = .R_bootstrap,
        .seed                 = .seed) 
    }
    
    ## Combine bootstrap results in one matrix
    ll_org <- lapply(.object, function(y) {
      if(inherits(.object, "cSEMResults_2ndorder")) {
        x    <- y$Second_stage$Information$Resamples$Estimates$Estimates1
        nobs <- nrow(y$First_stage$Information$Data)
      } else {
        x    <- y$Estimates$Estimates_resample$Estimates1
        nobs <- nrow(y$Information$Data)
      }
      path_resamples    <- x$Path_estimates$Resampled
      loading_resamples <- x$Loading_estimates$Resampled
      weight_resamples  <- x$Weight_estimates$Resampled
      n                 <- nrow(path_resamples)
      
      # Calculation of the bootstrap SEs
      ses <- infer(.object=y,.quantity = "sd")
      
      path_se <- ses$Path_estimates$sd 
      loading_se <- ses$Loading_estimates$sd
      weight_se <- ses$Weight_estimates$sd
      
      # Calculation of the bias
      bias <- infer(.object=y,.quantity = "bias")
      
      path_bias <- bias$Path_estimates$bias
      loading_bias <- bias$Loading_estimates$bias
      weight_bias <- bias$Weight_estimates$bias
      # Return object
      list(
        "n"                 = n,
        "nObs"              = nobs,
        "para_all"          = cbind(path_resamples,loading_resamples,weight_resamples),
        "ses_all"           = c(path_se, loading_se, weight_se),
        "bias_all"          = c(path_bias, loading_bias, weight_bias)
       )
      
    })
    
    
    ## Hensler (2009) approach
    if(any(.approach_mgd %in% c("all", "Henseler"))) {
      # center the bootstrap sample Sarstedt et al. (2011) Eq. 4
      ll_centered <- lapply(ll_org, function(x){

        # Substract from each row of para_all the corresponding bias
        t(apply(x$para_all,1,function(row){
          row-x$bias_all
        }))
      })
      
      # Create pairs which should be compared
      pairs_centered <- utils::combn(ll_centered, 2, simplify = FALSE)
      
      
      prob_Henseler <- lapply(pairs_centered,function(x){
        calculatePr(.resample_centered = x,
                    .parameters_to_compare = names_param)
      })
    }
    
    
    ## Keil and Nitzl approach 
    if(any(.approach_mgd %in% c("all", "Nitzl", "Keil"))) {
      diff_para_Keil <- diff_para_Nitzl <- calculateParameterDifference(
        .object = .object, 
        .model  = .parameters_to_compare
        )
      
      # permute .objects
      object_permu <- utils::combn(.object, 2, simplify = FALSE)
      names(object_permu) <- sapply(object_permu, function(x) paste0(names(x)[1], '_', names(x)[2]))
 
      # Approach suggested by Keil 
      if(any(.approach_mgd %in% c("all", "Keil"))) {     
        
      teststat_Keil <- lapply(names(object_permu), function(x) {
       diff <- diff_para_Keil[[x]]
        ses1 <- ll_org[[names(object_permu[[x]][1])]]$ses_all
        ses2 <- ll_org[[names(object_permu[[x]][2])]]$ses_all
        
        n1<-ll_org[[names(object_permu[[x]][1])]]$nObs
        n2<-ll_org[[names(object_permu[[x]][2])]]$nObs
        
        # Calculation of the SE of the parameter difference as proposed by 
        # Henseler (2007a), Henseler et al. (2009)
        ses_total <- sqrt((n1-1)^2/(n1+n2-2)*ses1^2 + 
            (n2-1)^2/(n1+n2-2)*ses2^2)*sqrt(1/n1+1/n2) 
          
        test_stat <- diff/ses_total[names(diff)]
        list("teststat" = test_stat, "df" = n1 + n2 - 2)
      })
      names(teststat_Keil) <- names(object_permu)
      teststat[["Keil"]] <- teststat_Keil 
      }
      
      # Approach suggested by Nitzl (2010)
      if(any(.approach_mgd %in% c("all", "Nitzl"))) { 
        
        teststat_Nitzl <- lapply(names(object_permu), function(x){
          diff <- diff_para_Nitzl[[x]]
          ses1 <- ll_org[[names(object_permu[[x]][1])]]$ses_all
          ses2 <- ll_org[[names(object_permu[[x]][2])]]$ses_all
          
          n1<-ll_org[[names(object_permu[[x]][1])]]$nObs
          n2<-ll_org[[names(object_permu[[x]][2])]]$nObs
          
          # Calculation of the SE of the parameter difference as proposed by 
          # Henseler (2007a), Henseler et al. (2009)
          ses_total <- sqrt((n1-1)/(n1)*ses1^2 + 
                              (n2-1)/(n2)*ses2^2) 
          
          test_stat <- diff/ses_total[names(diff)]
          
          # calculation of the degrees of freedom
          numerator <- ((n1-1)/n1*ses1^2+(n2-1)/n2*ses2^2)^2
          denominator <- (n1-1)/n1^2*ses1^4+(n2-1)/n2^2*ses2^4
          df <- round(numerator/denominator-2)
          list("teststat" = test_stat, "df" = df)
        })
        names(teststat_Nitzl) <- names(object_permu)
        teststat[["Nitzl"]] <- teststat_Nitzl 
      }
      
    }
    ## Sarstedt approach
    if(any(.approach_mgd %in% c("all", "Sarstedt"))) {
      ## Transpose
      ll <- purrr::transpose(ll_org)
      group_id <- rep(1:length(.object), unlist(ll$n))

      all_comb <- cbind(do.call(rbind,ll$para_all), 
                        "group_id" = group_id)
      
      # all_comb contains all parameter estimate that could potentially be compared 
      # plus an id column indicating the group adherance of each row.
      
      
      ## Select relevant columns
      all_comb <- all_comb[, c(names_param, "group_id")]
      
      ## Add test statistic Sarstedt
      teststat[["Sarstedt"]] <- calculateFR(.resample_sarstedt = all_comb)
    } # END approach_mgd %in% c("all", "Sarstedt")
  } # END .approach_mgd %in% c("all", "Sarstedt", "Keil", "Nitzl")
  
  ### Permutation ==============================================================
  if(.verbose) {
    cat("Start permutation:\n\n")
  }
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
  
  # Save old seed and restore on exit! This is important since users may have
  # set a seed before calling testMGD, in which case the global seed would be
  # overwritten by testMGD if not explicitly restored
  old_seed <- .Random.seed
  on.exit({.Random.seed <<- old_seed})
  
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
          .model  = .parameters_to_compare)
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
    
    # Update progres bar
    if(.verbose){
      setTxtProgressBar(pb, counter)
    }
    
    # Break repeat loop if .R results have been created.
    if(length(ref_dist) == .R_permutation) {
      break
    } else if(counter + n_inadmissibles == 10000) { 
      # Stop if 10000 runs did not result in insufficient admissible results
      stop("Not enough admissible result.", call. = FALSE)
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
    decision_Sarstedt <- lapply(padjusted_Sarstedt,function(p_value){
      temp <- lapply(.alpha, function(alpha){
        p_value > alpha
      })
      names(temp) <- paste(.alpha*100,"%",sep= '')
      temp
    })
    
    # Decision overall
    decision_overall_Sarstedt <- lapply(decision_Sarstedt, function(x){#over p-value adjustments
      lapply(x, function(xx){ #over different significant levels
        all(xx)
        })
    })
  }
  
  ## Keil et al. 2000 ----------------------------------------------------------
  if(any(.approach_mgd %in% c("all", "Keil"))){

    # Calculate p-values
    pvalue_Keil <- lapply(teststat$Keil,function(x){
      p_value <- 2*(1-pt(abs(x$teststat),df = x$df))
      p_value
    })
    
    # Adjust p-values in case of more than one comparison
    padjusted_Keil<- lapply(as.list(.approach_p_adjust), function(x){
      pvector <- stats::p.adjust(unlist(pvalue_Keil),method = x)
      # Sort them back into list
      relist(flesh = pvector,skeleton = pvalue_Keil)
    })
    names(padjusted_Keil) <- .approach_p_adjust
    
    # Decision 
    decision_Keil <- lapply(padjusted_Keil, function(adjust_approach){ # over the different p adjustments
      temp <- lapply(.alpha, function(alpha){# over the different significance levels
        lapply(adjust_approach,function(group_comp){# over the different group comparisons
          # check whether the p values are larger than a certain alpha
          group_comp > alpha
        })
      })
      names(temp) <- paste0(.alpha*100, "%")
      temp
    })
    
    # One rejection leads to overall rejection
    decision_overall_Keil <- lapply(decision_Keil, function(decision_Keil_list){
      lapply(decision_Keil_list,function(x){
        all(unlist(x))
      })
    })
  }
  
  if(any(.approach_mgd %in% c("all", "Nitzl"))){
    
    # Calculate p-values
    pvalue_Nitzl <- lapply(teststat$Nitzl,function(x){
      p_value <- 2*(1-pt(abs(x$teststat),df = x$df))
      p_value
    })
    
    # Adjust p-values in case of more than one comparison
    padjusted_Nitzl<- lapply(as.list(.approach_p_adjust), function(x){
      pvector <- stats::p.adjust(unlist(pvalue_Nitzl),method = x)
      # Sort them back into list
      relist(flesh = pvector,skeleton = pvalue_Nitzl)
    })
    names(padjusted_Nitzl) <- .approach_p_adjust
    
    # Decision 
    decision_Nitzl <- lapply(padjusted_Nitzl, function(adjust_approach){ # over the different p adjustments
      temp <- lapply(.alpha, function(alpha){# over the different significance levels
        lapply(adjust_approach,function(group_comp){# over the different group comparisons
          # check whether the p values are larger than a certain alpha
          group_comp > alpha
        })
      })
      names(temp) <- paste0(.alpha*100, "%")
      temp
    })
    
    # One rejection leads to overall rejection
    decision_overall_Nitzl <- lapply(decision_Nitzl, function(decision_Nitzl_list){
      lapply(decision_Nitzl_list,function(x){
        all(unlist(x))
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
    "Permutation_seed"      = .seed,
    # "Bootstrap_seeds"        = bootstrap_seeds,
    "Alpha"                 = .alpha,
    "VCV_type"              = .type_vcv,
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
      "P_value"            = padjusted_Chin,
      "Decision"           = decision_Chin,
      "Decision_overall"   = decision_overall_Chin
    )
    
    out[["Information"]][["Permutation_values"]][["Chin"]] <- ref_dist_matrices_Chin
  }
  
  if(any(.approach_mgd %in% c("all", "Sarstedt"))) {
    out[["Sarstedt"]] <- list(
      "Test_statistic"     = teststat_Sarstedt,
      "P_value"            = padjusted_Sarstedt,
      "Decision"           = decision_Sarstedt,
      "Decision_overall"   = decision_overall_Sarstedt
    )
    
    out[["Information"]][["Permutation_values"]][["Sarstedt"]] <- ref_dist_matrix_Sarstedt
  }
  if(any(.approach_mgd %in% c("all", "Keil"))) {
    out[["Keil"]] <- list(
      "Test_statistic"     = purrr::transpose(teststat_Keil)$teststat,
      "P_value"            = padjusted_Keil,
      "Decision"           = decision_Keil,
      "Decision_overall"   = decision_overall_Keil,
      "df"                 = purrr::transpose(teststat_Keil)$df[[1]]   
    )
  }
  
  if(any(.approach_mgd %in% c("all", "Nitzl"))) {
    out[["Nitzl"]] <- list(
      "Test_statistic"     = purrr::transpose(teststat_Nitzl)$teststat,
      "P_value"            = padjusted_Nitzl,
      "Decision"           = decision_Nitzl,
      "Decision_overall"   = decision_overall_Nitzl,
      "df"                 = purrr::transpose(teststat_Nitzl)$df[[1]]
      
    )
  }
  
  
  class(out) <- "cSEMTestMGD"
  return(out)
}
