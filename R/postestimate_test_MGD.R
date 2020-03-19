#' Tests for multi-group comparisons
#'
#' This function performs several permutation tests, i.e., the reference distribution 
#' of the test statistic is obtained by permutation.
#' 
#' The following tests are implemented:
#' \describe{
#' \item{`.approach_mgd="Klesel"`: Approach suggested by \insertCite{Klesel2019;textual}{cSEM}}{
#'   The model-implied variance-covariance matrix (either indicator 
#'   (`.type_vcv = "indicator"`) or construct (`.type_vcv = "construct"`)) 
#'   is compared across groups. 
#' 
#'   To measure the distance between the model-implied variance-covariance matrices, 
#'   the geodesic distance (dG) and the squared Euclidean distance (dL) are used.
#'   If more than two groups are compared, the average distance over all groups
#'   is used.}
#' \item{`.approach_mgd="Sarstedt"`: Approach suggested by \insertCite{Sarstedt2011;textual}{cSEM}}{
#'   Groups are compared in terms of parameter differences across groups.
#'   \insertCite{Sarstedt2011;textual}{cSEM} tests if parameter k is equal
#'   across all groups. If several parameters are tested simultaneously  
#'   it is recommended to adjust the significance  level or the p-values (in \pkg{cSEM} correction is
#'   done by p-value). By default
#'   no multiple testing correction is done, however, several common
#'   adjustments are available via `.approach_p_adjust`. See 
#'   \code{\link[stats:p.adjust]{stats::p.adjust()}} for details. Note: the 
#'   test has some severe shortcomings. Use with caution.
#' }
#' \item{`.approach_mgd="Chin"`: Approach suggested by \insertCite{Chin2010;textual}{cSEM}}{
#'   Groups are compared in terms of parameter differences across groups.
#'   \insertCite{Chin2010;textual}{cSEM} tests if parameter k is equal
#'   between two groups. If more than two groups are tested for equality, parameter 
#'   k is compared between all pairs of groups. In this case, it is recommended
#'   to adjust the significance  level or the p-values (in \pkg{cSEM} correction is
#'   done by p-value) since this is essentially a multiple testing setup. 
#'   If several parameters are tested simultaneously, correction is by group 
#'   and number of parameters. By default
#'   no multiple testing correction is done, however, several common
#'   adjustments are available via `.approach_p_adjust`. See 
#'   \code{\link[stats:p.adjust]{stats::p.adjust()}} for details.
#' }
#' \item{`.approach_mgd="Keil"`: Approach suggested by \insertCite{Keil2000;textual}{cSEM}}{
#'   Groups are compared in terms of parameter differences across groups.
#'   \insertCite{Keil2000;textual}{cSEM} tests if parameter k is equal
#'   between two groups. It is assumed, that the standard errors of the coefficients are 
#'   equal across groups. The calculation of the standard error of the parameter 
#'   difference is adjusted as proposed by \insertCite{Henseler2009;textual}{cSEM}.
#'   If more than two groups are tested for equality, parameter k is compared 
#'   between all pairs of groups. In this case, it is recommended
#'   to adjust the significance  level or the p-values (in \pkg{cSEM} correction is
#'   done by p-value) since this is essentially a multiple testing setup. 
#'   If several parameters are tested simultaneously, correction
#'   is by group and number of parameters. By default
#'   no multiple testing correction is done, however, several common
#'   adjustments are available via `.approach_p_adjust`. See 
#'   \code{\link[stats:p.adjust]{stats::p.adjust()}} for details.
#' }
#' \item{`.approach_mgd="Nitzl"`: Approach suggested by \insertCite{Nitzl2010;textual}{cSEM}}{
#'   Groups are compared in terms of parameter differences across groups.
#'   Similarly to \insertCite{Keil2000;textual}{cSEM}, a single parameter k is tested
#'   for equality between two groups. In contrast to \insertCite{Keil2000;textual}{cSEM},
#'   it is assumed, that the standard errors of the coefficients are 
#'   unequal across groups \insertCite{Sarstedt2011}{cSEM}.
#'   If more than two groups are tested for equality, parameter k is compared 
#'   between all pairs of groups. In this case, it is recommended
#'   to adjust the significance  level or the p-values (in \pkg{cSEM} correction is
#'   done by p-value) since this is essentially a multiple testing setup. 
#'   If several parameters are tested simultaneously, correction
#'   is by group and number of parameters. By default
#'   no multiple testing correction is done, however, several common
#'   adjustments are available via `.approach_p_adjust`. See 
#'   \code{\link[stats:p.adjust]{stats::p.adjust()}} for details.
#' }
#' \item{`.approach_mgd="CI_param"`: Approach mentioned in \insertCite{Sarstedt2011;textual}{cSEM}}{
#'   This approach is based on the confidence intervals constructed around the 
#'   parameter estimates of the two groups. If the parameter of one group falls within 
#'   the confidence interval of the other group and/or vice versa, it can be concluded 
#'   that there is no group difference.   
#'   Since it is based on the confidence intervals `.approach_p_adjust` is ignored.
#' } 
#' 
#' \item{`.approach_mgd="CI_overlap"`: Approach mentioned in \insertCite{Sarstedt2011;textual}{cSEM}}{
#'   This approach is based on the confidence intervals (CIs) constructed around the 
#'   parameter estimates of the two groups. If the two CIs overlap, it can be concluded 
#'   that there is no group difference.   
#'   Since it is based on the confidence intervals `.approach_p_adjust` is ignored.
#' } 
#' 
#' 
#' }
#' 
#' Use `.approach_mgd` to choose the approach. By default all approaches are computed
#' (`.approach_mgd = "all"`).
#' 
#' By default, approaches based on parameter differences across groups compare
#' all parameters (`.parameters_to_compare = NULL`). To compare only
#' a subset of parameters provide the parameters in [lavaan model syntax][lavaan::model.syntax]  just like
#' the model to estimate. Take the simple model:
#' 
#' \preformatted{
#' model_to_estimate <- "
#' Structural model
#' eta2 ~ eta1
#' eta3 ~ eta1 + eta2
#' 
#' # Each concept os measured by 3 indicators, i.e., modeled as latent variable
#' eta1 =~ y11 + y12 + y13
#' eta2 =~ y21 + y22 + y23
#' eta3 =~ y31 + y32 + y33
#' "
#' }
#' If only the path from eta1 to eta3 and the loadings of eta1 are to be compared
#' across groups, write:
#' \preformatted{
#'to_compare <- "
#' Structural parameters to compare
#' eta3 ~ eta1
#' 
#' # Loadings to compare
#' eta1 =~ y11 + y12 + y13
#' "
#' }
#' Note that the "model" provided to `.parameters_to_compare`
#' does not have to be an estimable model! 
#' 
#' Note also that compared to all other functions in \pkg{cSEM} using the argument,
#'  `.handle_inadmissibles` defaults to `"replace"` to accomdate the Sarstedt et al. (2011) approach.
#' 
#' Argument `.R_permuation` is ignored for the `"Nitzl"` and the `"Keil"` approach. 
#' `.R_bootstrap` is ignored if  `.object` already contains resamples, 
#' i.e. has class `cSEMResults_resampled` and if only the `"Klesel"` or the `"Chin"`
#' approach are used.
#' 
#' The argument `.saturated` is used by `"Klesel"` only. If `.saturated = TRUE` 
#' the original structural model is ignored and replaced by a saturated model, 
#' i.e. a model in which all constructs are allowed to correlate freely. 
#' This is useful to test differences in the measurement models between groups 
#' in isolation.
#' 
#' @usage testMGD(
#'  .object                = NULL,
#'  .alpha                 = 0.05,
#'  .approach_p_adjust     = "none",
#'  .approach_mgd          = c("all", "Klesel", "Chin", "Sarstedt", 
#'                             "Keil", "Nitzl", "Henseler", "CI_para","CI_overlap"),
#'  .parameters_to_compare = NULL,
#'  .handle_inadmissibles  = c("replace", "drop", "ignore"),
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
#'   yielded inadmissible results (i.e. number of results returned is equal to `.R`). 
#'   For "*replace*" resampling continues until there are exactly `.R` admissible solutions. 
#'   Defaults to "*replace*" to accommodate all approaches.
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
#'   \item{`$Henseler`}{A list with elements, `Test_statistic`, `P_value`, `Decision`, and `Decision_overall`}
#' }
#' @references
#'   \insertAllCited{}
#'   
#' @seealso [csem()], [cSEMResults], [testMICOM()], [testOMF()]
#'
#' @example inst/examples/example_testMGD.R
#' 
#' @export

testMGD <- function(
 .object                = NULL,
 .alpha                 = 0.05,
 .approach_p_adjust     = "none",
 .approach_mgd          = c("all", "Klesel", "Chin", "Sarstedt", 
                            "Keil", "Nitzl","Henseler", "CI_para","CI_overlap"),
 .parameters_to_compare = NULL,
 .handle_inadmissibles  = c("replace", "drop", "ignore"),
 .R_permutation         = 499,
 .R_bootstrap           = 499,
 .saturated             = FALSE,
 .seed                  = NULL,
 .type_ci               = c("CI_percentile","CI_basic"),
 .type_vcv              = c("indicator", "construct"),
 .verbose               = TRUE
){

  ### Checks and errors ========================================================
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
  if(!all(.type_ci %in%  c("CI_standard_z", "CI_standard_t", "CI_percentile",
                      "CI_basic", "CI_bc", "CI_bca") )){
    stop2("The specified confidence interval in .typ.cv is not valid.\n",
          "Please choose one of the following: CI_standard_z, CI_standard_t,\n",
          "CI_percentile, CI_basic, CI_bc, CI_bca.")
  }
  

  # Check if at least two groups are present
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
          " to nonlinear models as cSEM currently can not calculate",
          " the model-implied VCV matrix for such models.\n", 
          "Consider setting `.approach_mgd = c('Chin',  'Sarstedt')`")
  }
  
  # Check if any of the group estimates are inadmissible
  if(sum(unlist(verify(.object))) != 0) {
    warning2(
      "The following warning occured in the testMGD() function:\n",
      "Initial estimation results for at least one group are inadmissible.\n", 
      "See `verify(.object)` for details.")
  }
  
  # Check if data for different groups is identical
  if(TRUE %in% lapply(utils::combn(.object, 2, simplify = FALSE),
                      function(x){ identical(x[[1]], x[[2]])})){
    warning2(
      "The following warning occured in the testMGD() function:\n",
      "At least two groups are identical. Results may not be meaningful.")
  } 
  
  # If Henseler, CI_para or CI_overlap are used with 
  # adjustment of p-value different than "none" return warning
  if(any(.approach_mgd %in% c("all", "Henseler","CI_para","CI_overlap")) & !all(.approach_p_adjust %in% "none")){
    warning2(
      "The following warning occured in the testMGD() function:\n",
      "Currently, there is no p-value adjustment possible for the approach suggested by\n",
      "Henseler (2007), CI_para, and CI_overlap. Every adjustment is ignored for this approach."
    )
  }
  
  if(.verbose) {
    cat(rule2("Several tests for multi-group comparisons",
              type = 3), "\n\n")
  }
  

  # Order significance levels
  .alpha <- .alpha[order(.alpha)]
  
  # Get the names of the parameters to be compared. 
  names_all_param <- getParameterNames(.object, .model = .parameters_to_compare)
  
  # getParameterNames() returns also measurement error and indicator correlations
  # Currently they cannot be handeled, therefore an error is returned. 
  # FOR FUTURE: MEasurement errors and indicator correlations can be added to the parameters allowed in the comparison  
  if(!is.null(names_all_param$names_cor_measurement_error)|!is.null(names_all_param$names_cor_indicator)){
    stop2("The following error occured in the testMGD() function:\n",
          "Currenlty it is not allowed to compare measurement error covariance",
          " and/or indicator covariances across groups.")
  }
  
  names_param <- unlist(names_all_param)
  
  ## Calculation of the test statistics========================================
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
  
  ## Sarstedt et al. (2011), Keil et al. (2000), Nitzl (2010), Henseler (2007) -----
  # All these approaches require bootstrap
  if(any(.approach_mgd %in% c("all", "Sarstedt", "Keil", "Nitzl", "Henseler", "CI_para","CI_overlap"))) {

    # Check if .object already contains resamples; if not, run bootstrap
    if(!inherits(.object, "cSEMResults_resampled")) {
      if(.verbose) {
        cat("Bootstrap cSEMResults object ...\n\n")
      }
      
      .object <- resamplecSEMResults(
        .object               = .object,
        .resample_method      = "bootstrap",
        .handle_inadmissibles = .handle_inadmissibles,
        .R                    = .R_bootstrap,
        .seed                 = .seed) 
    }
    
    ## Combine bootstrap results in one matrix
    bootstrap_results <- lapply(.object, function(y) {
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
      # Rename
      cor_exo_cons_resamples <- x$Exo_construct_correlation$Resampled
      n                 <- nrow(path_resamples)
      
      # Calculation of the bootstrap SEs
      ses <- infer(.object=y,.quantity = "sd")
      
      path_se <- ses$Path_estimates$sd 
      loading_se <- ses$Loading_estimates$sd
      weight_se <- ses$Weight_estimates$sd
      cor_exo_cons_se <- ses$Exo_construct_correlation$sd
      
      # Calculation of the bias
      bias <- infer(.object=y,.quantity = "bias")
      
      path_bias <- bias$Path_estimates$bias
      loading_bias <- bias$Loading_estimates$bias
      weight_bias <- bias$Weight_estimates$bias
      cor_exo_cons_bias <- bias$Exo_construct_correlation$bias
      
      # Structure output
      out<-list(
        "n"                 = n,
        "nObs"              = nobs,
        "para_all"          = cbind(path_resamples,loading_resamples,weight_resamples, cor_exo_cons_resamples),
        "ses_all"           = c(path_se, loading_se, weight_se, cor_exo_cons_se),
        "bias_all"          = c(path_bias, loading_bias, weight_bias, cor_exo_cons_bias)
        )
      
      # If comparison should be done via CIs 
      if(any(.approach_mgd %in% c("all", "CI_para", "CI_overlap"))){
        
        diff <- setdiff(.type_ci, args_default(TRUE)$.type_ci)
        
        if(length(diff) != 0) {
          stop2(
            "The following error occured in the testMGD() function:\n",
            "Unknown approach: ", paste0(diff, collapse = ", "), ".",
            " Possible choices are: ",
            paste0(args_default(TRUE)$.type_ci, collapse = ", "))
        }
        
        # calculation of the CIs
        cis <- infer(.object=y, .quantity=.type_ci, .alpha = .alpha)
        cis_temp <- purrr::transpose(cis)
        cis_ret<-lapply(cis_temp,function(x){
        path_ci <- x$Path_estimates
        loading_ci <- x$Loading_estimates
        weight_ci <- x$Weight_estimates
        cor_exo_cons_ci <- x$Exo_construct_correlation
        # deliberately name them xyz_estimates to be able to select them later via names.
        list(path_estimates = path_ci, loading_estimates = loading_ci,
             weight_estimates = weight_ci, cor_exo_cons_estimates = cor_exo_cons_ci)
        })
        names(cis_ret) <- names(cis_temp)
         
        out[["ci_all"]] <- cis_ret
      }
      
      # Return output
      return(out)
      
    })
    
    ## Keil, Nitzl, Henseler 
    if(any(.approach_mgd %in% c("all", "Nitzl", "Keil","Henseler"))) {
      
      # Calculate the difference for one parameter between two groups
      # Although Henseler approach does not compute a test statistic,
      # The difference is calculated to create a dummy test statistic list
      diff_para_Keil <- diff_para_Nitzl <- diff_para_Henseler<- calculateParameterDifference(
        .object = .object, 
        .model  = .parameters_to_compare
        )
      
      
      if(any(.approach_mgd %in% c("all","Henseler"))) {
        # Create dummy test statistic list containing NAs
        temp <- rep(NA,length(unlist(diff_para_Henseler)))
        teststat[["Henseler"]] <-relist(flesh = temp,skeleton = diff_para_Henseler)
      }
      
      # Build list that contains pairs of objects.
      object_permu <- utils::combn(.object, 2, simplify = FALSE)
      names(object_permu) <- sapply(object_permu, function(x) paste0(names(x)[1], '_', names(x)[2]))
 
      # Approach suggested by Keil 
      if(any(.approach_mgd %in% c("all", "Keil"))) {     
        
      teststat_Keil <- lapply(names(object_permu), function(x) {
       diff <- diff_para_Keil[[x]]
        ses1 <- bootstrap_results[[names(object_permu[[x]][1])]]$ses_all
        ses2 <- bootstrap_results[[names(object_permu[[x]][2])]]$ses_all
        
        # Get sample size per group
        n1<-bootstrap_results[[names(object_permu[[x]][1])]]$nObs
        n2<-bootstrap_results[[names(object_permu[[x]][2])]]$nObs
        
        # Calculation of the SE of the parameter difference as proposed by 
        # Henseler (2007a), Henseler et al. (2009) 
        # Assumption of equal variances
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
          ses1 <- bootstrap_results[[names(object_permu[[x]][1])]]$ses_all
          ses2 <- bootstrap_results[[names(object_permu[[x]][2])]]$ses_all
          
          # Get sample size per group
          n1<-bootstrap_results[[names(object_permu[[x]][1])]]$nObs
          n2<-bootstrap_results[[names(object_permu[[x]][2])]]$nObs
          
          # Calculation of the SE of the parameter difference as proposed by 
          # Henseler (2007a), Henseler et al. (2009) 
          # Assumption of unequal variances
          ses_total <- sqrt((n1-1)/(n1)*ses1^2 + 
                              (n2-1)/(n2)*ses2^2) 
          
          test_stat <- diff/ses_total[names(diff)]
          
          # calculation of the degrees of freedom
          numerator <- ((n1-1)/n1*ses1^2+(n2-1)/n2*ses2^2)^2
          denominator <- (n1-1)/n1^2*ses1^4+(n2-1)/n2^2*ses2^4
          df <- round(numerator/denominator-2)
          df <- df[names(diff)]
          list("teststat" = test_stat, "df" = df)
        })
        names(teststat_Nitzl) <- names(object_permu)
        teststat[["Nitzl"]] <- teststat_Nitzl 
      }
      
    }
    ## Sarstedt et al. approach
    if(any(.approach_mgd %in% c("all", "Sarstedt"))) {
      ## Transpose
      ll <- purrr::transpose(bootstrap_results)
      group_id <- rep(1:length(.object), unlist(ll$n))

      # all_comb contains all parameter estimate that could potentially be compared 
      # plus an id column indicating the group adherance of each row.
      all_comb <- cbind(do.call(rbind,ll$para_all), 
                        "group_id" = group_id)

      ## Select relevant columns
      all_comb <- all_comb[, c(names_param, "group_id")]
      
      ## Add test statistic Sarstedt
      teststat[["Sarstedt"]] <- calculateFR(.resample_sarstedt = all_comb)
    } # END approach_mgd %in% c("all", "Sarstedt")
  } # END .approach_mgd %in% c("all", "Sarstedt", "Keil", "Nitzl", "Henseler")
  
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
  
  # Permutation is only performed for approaches that require it
  if(any(.approach_mgd %in% c("all", "Klesel", "Chin", "Sarstedt"))) {
    
  # Start progress bar if required
  if(.verbose){
    cat("Start permutation:\n\n")
    pb <- txtProgressBar(min = 0, max = .R_permutation, style = 3)
  }
  
    # Save old seed and restore on exit! This is important since users may have
    # set a seed before, in which case the global seed would be
    # overwritten if not explicitly restored
    old_seed <- .Random.seed
    on.exit({.Random.seed <<- old_seed})
    
    ## Create seed if not already set
    if(is.null(.seed)) {
      set.seed(seed = NULL)
      # Note (08.12.2019): Its crucial to call set.seed(seed = NULL) before
      # drawing a random seed out of .Random.seed. If set.seed(seed = NULL) is not
      # called sample(.Random.seed, 1) would result in the same random seed as
      # long as .Random.seed remains unchanged. By resetting the seed we make 
      # sure that sample draws a different element everytime it is called.
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
  }
  
  # # Order significance levels
  # .alpha <- .alpha[order(.alpha)]
  
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
  # TRUE = p-value > alpha --> not reject
  # FALSE = sufficient evidence against the H0 --> reject
  decision_Klesel <- lapply(.alpha, function(x) {
    pvalue_Klesel > x
  })
  
  names(decision_Klesel) <- paste0(.alpha * 100, "%")
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
    # TRUE = p-value > alpha --> not reject
    # FALSE = sufficient evidence against the H0 --> reject
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
    # TRUE = p-value > alpha --> not reject
    # FALSE = sufficient evidence against the H0 --> reject
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
  
  ## Keil et al. (2000) --------------------------------------------------------
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
    # TRUE = p-value > alpha --> not reject
    # FALSE = sufficient evidence against the H0 --> reject
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
  
  ## Nitzl et al. (2010) -------------------------------------------------------
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
    # TRUE = p-value > alpha --> not reject
    # FALSE = sufficient evidence against the H0 --> reject
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
  
  ## Henseler (2009) approach  -------------------------------------------------
  if(any(.approach_mgd %in% c("all", "Henseler"))) {
    
    # Pseudo teststat for convenience
    teststat_Henseler <- teststat$Henseler
    
    # center the bootstrap sample Sarstedt et al. (2011) Eq. 4
    ll_centered <- lapply(bootstrap_results, function(x){
      
      # Substract from each row of para_all the corresponding bias
      t(apply(x$para_all,1,function(row){
        row-x$bias_all
      }))
    })
    
    names(ll_centered) <- names(bootstrap_results)
    
    # Create pairs which should be compared
    pairs_centered <- utils::combn(ll_centered, 2, simplify = FALSE)
    names(pairs_centered) <- sapply(pairs_centered, function(x) paste0(names(x)[1], '_', names(x)[2]))

    
    # Calculation of the probability
    pvalue_Henseler <- lapply(pairs_centered,function(x){
      calculatePr(.resample_centered = x,
                  .parameters_to_compare = names_param)
    })
    
    # Adjust p-value in case of multiple comparisons
    # Adjusting p-values is not straight forward. 
    # First it is a one-sided test, so we need to know the hypothesis
    # Just flipping the p-value is not reaaly an option as it might causes problems 
    # in situation where the order of the p-values is required for the correction
    # Therefore only the "none"method is applied
    padjusted_Henseler<- lapply(as.list("none"), function(x){
      pvector <- stats::p.adjust(unlist(pvalue_Henseler),method = x)
      # Sort them back into list
      relist(flesh = pvector,skeleton = pvalue_Henseler)
    })
    names(padjusted_Henseler) <- "none"

    # Decision is made:
    # The probability is compared to alpha and 1-alpha
    # In doing so, we assess both hypotheses theta^1 <= theta^2 and 
    # theta^1 => theta^2
    # TRUE = p-value > alpha & p-value < 1-alpha --> not reject
    # FALSE = sufficient evidence against the H0 --> reject
    decision_Henseler <- lapply(padjusted_Henseler, function(adjust_approach){ # over the different p adjustments
      temp <- lapply(.alpha, function(alpha){# over the different significance levels
        lapply(adjust_approach,function(group_comp){# over the different group comparisons
          # check whether the p values are larger than a certain alpha
          group_comp > alpha & group_comp < 1- alpha
        })
      })
      names(temp) <- paste0(.alpha*100, "%")
      temp
    })
    
    # It is not clear how an overall decision should be made  
    decision_overall_Henseler <- lapply(decision_Henseler, function(decision_Henseler_list){
      lapply(decision_Henseler_list,function(x){
        NA
      })
    })
  }
    
  
  # Comparison via confidence intervals
  if(any(.approach_mgd %in% c("all", "CI_para", "CI_overlap"))) {
    # Select CIs from the bootstrap_results
    cis <- lapply(bootstrap_results,function(x){
      x$ci_all
    })
    # Make group pairs of CIs
    cis_comp <- utils::combn(cis, 2, simplify = FALSE)
    names(cis_comp) <- sapply(cis_comp, function(x) paste0(names(x)[1], '_', names(x)[2]))

    # Select the relevant parameters per group and make group pairs
    # This is required fo both CI_para and CI_overlap
    # For the latter it is used to select the appropriate parameters
    param_per_group <- getRelevantParameters(.object = .object,
                                             .model = .parameters_to_compare)
    param_per_group <- purrr::transpose(param_per_group)
    param_comp <- utils::combn(param_per_group, 2, simplify = FALSE)
    names(param_comp) <- sapply(param_comp,
                                function(x)
                                  paste0(names(x)[1], '_', names(x)[2]))

    if (any(.approach_mgd %in% c("all", "CI_para"))) {
      
      # Investigate whether the estimate of one group is part of the CI of another group
      decision_ci_para <- lapply(.alpha, function(alpha) {
        # lapply over the comparisons
        tttt <- lapply(names(cis_comp), function(comp) {
            # lapply over the type of CI for the first group
            # As group1 and group2 have the same structure the selection does not matter
            tt = lapply(names(cis_comp[[comp]][[1]]), function(interval_type){
              # lapply over the different parameters
                ttt <- lapply(names(cis_comp[[comp]][[1]][[interval_type]]), 
                              function(param) {

                    lb <- paste0(100 * (1 - alpha), "%L")
                    ub <- paste0(100 * (1 - alpha), "%U")

                    # Estimate first group
                    temp_para1 <- param_comp[[comp]][[1]][[param]]
                    # CIs first group
                    temp_cis1 <- cis_comp[[comp]][[1]][[interval_type]][[param]]
                    temp_cis_selected1 <- temp_cis1[, names(temp_para1), drop = FALSE]
                    lb_temp1 <-temp_cis_selected1[lb,] 
                    names(lb_temp1) <- colnames(temp_cis_selected1[lb,,drop=FALSE])
                    ub_temp1 <-temp_cis_selected1[ub,] 
                    names(ub_temp1) <- colnames(temp_cis_selected1[ub,,drop=FALSE])
                    
                    # Estimate second group
                    temp_para2 <- param_comp[[comp]][[2]][[param]]
                    # CIs second group
                    temp_cis2 <- cis_comp[[comp]][[2]][[interval_type]][[param]]
                    temp_cis_selected2 <- temp_cis2[, names(temp_para2), drop = FALSE]
                    lb_temp2 <-temp_cis_selected2[lb,] 
                    names(lb_temp2) <- colnames(temp_cis_selected2[lb,,drop=FALSE])
                    ub_temp2 <-temp_cis_selected2[ub,] 
                    names(ub_temp2) <- colnames(temp_cis_selected2[ub,,drop=FALSE])

                    # # check whether parameter estimates of one group 
                    # falls within boundaries of the CIs of the other group
                    # TRUE parameter estimate of the other group falls within the CI
                    # => no group difference
                    # Otherwise, FALSE => group difference
                      decision<-
                        (lb_temp2 < temp_para1 & 
                           ub_temp2 > temp_para1) |
                        (lb_temp1 < temp_para2 & 
                           ub_temp1 > temp_para2)

                      # Structure output
                      out=data.frame("Estimate"=temp_para1,"lb"=lb_temp2,
                                     "ub"=ub_temp2,"Estimate"=temp_para2,"lb"=lb_temp1,
                                     "ub"=ub_temp1,"decision"=decision)
                      
                      # that can be solved more elegant
                      # If there is no parameter to compare the data.frame has zero rows
                      # If the if is removed there are problems as the number of columns does not match 
                      # the length of the name vector
                      if(nrow(out)!=0){
                      colnames(out)=c(paste0("Est_",names(cis_comp[[comp]])[1]),
                                      paste0("lb_",names(cis_comp[[comp]])[2]),
                                      paste0("ub_",names(cis_comp[[comp]])[2]),
                                      paste0("Est_",names(cis_comp[[comp]])[2]),
                                      paste0("lb_",names(cis_comp[[comp]])[1]),
                                      paste0("ub_",names(cis_comp[[comp]])[1]),
                                      "Decision")
                      }
                      out
                  })

                do.call(rbind,ttt)
            })
            names(tt) = names(cis_comp[[comp]][[1]])
            tt
        })
        names(tttt) <- names(cis_comp)
        tttt
      })
      names(decision_ci_para) <- paste0((1 - .alpha) * 100, "%")

      
      # Decision overall
      decision_overall_ci_para<-lapply(names(decision_ci_para), function(alpha) {
        t <- lapply(names(decision_ci_para[[alpha]]), function(comp) {
          tt <- lapply(names(decision_ci_para[[alpha]][[comp]]),function(interval_type){
            # Check whether there is one FALSE among the decisions
            all(decision_ci_para[[alpha]][[comp]][[interval_type]][,"Decision"])
          })
          names(tt) <- names(decision_ci_para[[alpha]][[comp]])
          tt
          })
          names(t) <- names(decision_ci_para[[alpha]])
          t
        })
      names(decision_overall_ci_para) <- names(decision_ci_para)
      
    }#end if .approach_mgd == CI_para
    
    if(any(.approach_mgd %in% c("all", "CI_overlap"))) {
      
      decision_ci_overlap <- lapply(.alpha, function(alpha) {
        # lapply over the group pairs
        ttt <- lapply(names(cis_comp), function(comp) {
          # lapply over the type of CI for the first group
          # As group1 and group2 have the same structure the selection does not matter
          t <- lapply(names(cis_comp[[comp]][[1]]), function(interval_type) {
            # lapply over the parameters that are compared
              tt <- lapply(names(cis_comp[[comp]][[1]][[interval_type]]), function(param) {
                lb <- paste0(100 * (1 - alpha), "%L")
                ub <- paste0(100 * (1 - alpha), "%U")

                # It does not matter which group is chosen to select the relevant parameters 
                # as in all groups the relevant parameters are the same
                para_rel <- names(param_comp[[comp]][[1]][[param]])

                # Select lower and upper bound of the CIs for the two groups
                lb1 <- cis_comp[[comp]][[1]][[interval_type]][[param]][lb,para_rel ]
                ub1 <- cis_comp[[comp]][[1]][[interval_type]][[param]][ub,para_rel ]
                lb2 <- cis_comp[[comp]][[2]][[interval_type]][[param]][lb,para_rel ]
                ub2 <- cis_comp[[comp]][[2]][[interval_type]][[param]][ub,para_rel ]
                
                # Check whether the boundaries of the CI of the first group fall
                # within the boundaries of the second group
                decision <- (lb2 < lb1 & lb1 < ub2) | 
                  (lb2 < ub1 & ub1 < ub2) | 
                  (lb1<lb2 & ub1>ub2) | 
                  (lb2<lb1 & ub2>ub1)

                # Structure output
                out <- data.frame(lb1,ub1,lb2,ub2,decision)
                
                # that can be solved more elegant
                # If there is no parameter to compare the data.frame has zero rows
                # If the if is removed there are problems as the number of columns does not match 
                # the length of the name vector
                if(nrow(out)!=0){
                colnames(out)=c(paste0("lb_",names(cis_comp[[comp]])[1]),
                  paste0("ub_",names(cis_comp[[comp]])[1]),
                  paste0("lb_",names(cis_comp[[comp]])[2]),
                  paste0("ub_",names(cis_comp[[comp]])[2]),
                  "Decision"
                )
                rownames(out)<-para_rel
                }
                out
                  })
            do.call(rbind,tt)
          })
          names(t) <- names(cis_comp[[comp]][[1]])
          t
        })
        names(ttt) <- names(cis_comp)
        ttt
      })
      names(decision_ci_overlap) <- paste0((1 - .alpha) * 100, "%")
      
      # Decision overall
      decision_overall_ci_overlap<-lapply(names(decision_ci_overlap), function(alpha) {
        t <- lapply(names(decision_ci_overlap[[alpha]]), function(comp) {
          tt <- lapply(names(decision_ci_overlap[[alpha]][[comp]]),function(interval_type){
            # Check whether there is one FALSE among the decisions
            all(decision_ci_overlap[[alpha]][[comp]][[interval_type]][,"Decision"])
          })
          names(tt) <- names(decision_ci_overlap[[alpha]][[comp]])
          tt
        })
        names(t) <- names(decision_ci_overlap[[alpha]])
        t
      })
      names(decision_overall_ci_overlap) <- names(decision_ci_overlap)
    }#end if: if one CI approach is requested
  }
  ### Return output ============================================================
  out <- list()
  
  ## General Information
  out[["Information"]] <- list(
    "Group_names"           = names(.object),
    "Number_of_observations"= sapply(X_all_list, nrow),
    "Approach"              = .approach_mgd,
    "Approach_p_adjust"     = .approach_p_adjust,
    "Alpha"                 = .alpha
  )

  ## Permutation-specific information
  if(any(.approach_mgd %in% c("all", "Klesel", "Chin", "Sarstedt"))) {
   out[["Information"]][["Information_permutation"]] <- list(
     "Number_admissibles"    = length(ref_dist1), 
     "Total_runs"            = counter + n_inadmissibles,
     "Permutation_seed"      = .seed,
     "Permutation_values"    = list(),
     "Handle_inadmissibles"  = .handle_inadmissibles
   ) 
  } else {
    out[["Information"]][["Information_permutation"]] <- list(
      "Number_admissibles"    = NA, 
      "Total_runs"            = NA,
      "Permutation_seed"      = NA,
      "Permutation_values"    = NA,
      "Handle_inadmissibles"  = NA
    ) 
  }
 
  ## Bootstrap-sprecific information
  if(any(.approach_mgd %in% c("all", "Sarstedt", "Keil", "Nitzl", "Henseler",
                              "CI_param","CI_overlap"))) {
    # Collect bootstrap information
    info_boot <-lapply(.object, function(x){
      if(inherits(.object, "cSEMResults_2ndorder")) {
        x <- x$Second_stage$Information$Resamples$Information_resample
      } else {
        x <- x$Information$Information_resample
      }
      
      list(
        "Number_admissibles"   = x$Number_of_admissibles,
        "Total_runs"           = x$Number_of_runs,
        "Bootstrap_seed"       = x$Seed,
        "Handle_inadmissibles" = x$Handle_inadmissibles
      )
    })
    
    names(info_boot) <- names(.object)
    info_boot <- purrr::transpose(info_boot)
    
    out[["Information"]][["Information_bootstrap"]] <- list(
      "Number_admissibles"    = info_boot$Number_admissibles, 
      "Total_runs"            = info_boot$Total_runs,
      "Bootstrap_seed"        = info_boot$Bootstrap_seed,
      "Handle_inadmissibles"  = info_boot$Handle_inadmissibles
    ) 
  }
   
# Information for Klesel et al. approach
  if(any(.approach_mgd %in% c("all", "Klesel"))) {
    out[["Klesel"]] <- list(
      "Test_statistic"     = teststat_Klesel,
      "P_value"            = pvalue_Klesel,
      "Decision"           = decision_Klesel,
      "VCV_type"           = .type_vcv
    )
    
    out[["Information"]][["Information_permutation"]][["Permutation_values"]][["Klesel"]] <- ref_dist_matrix_Klesel
  }
  
# Information for Chin & Dibbern approach
  if(any(.approach_mgd %in% c("all", "Chin"))) {
    out[["Chin"]] <- list(
      "Test_statistic"     = teststat_Chin,
      "P_value"            = padjusted_Chin,
      "Decision"           = decision_Chin,
      "Decision_overall"   = decision_overall_Chin
    )
    
    out[["Information"]][["Information_permutation"]][["Permutation_values"]][["Chin"]] <- ref_dist_matrices_Chin
  }

# Information for Sarstedt et al. approach  
  if(any(.approach_mgd %in% c("all", "Sarstedt"))) {
    out[["Sarstedt"]] <- list(
      "Test_statistic"     = teststat_Sarstedt,
      "P_value"            = padjusted_Sarstedt,
      "Decision"           = decision_Sarstedt,
      "Decision_overall"   = decision_overall_Sarstedt
    )
    
    out[["Information"]][["Information_permutation"]][["Permutation_values"]][["Sarstedt"]] <- ref_dist_matrix_Sarstedt
  }
  
# Information for Keil approach  
  if(any(.approach_mgd %in% c("all", "Keil"))) {
    out[["Keil"]] <- list(
      "Test_statistic"     = purrr::transpose(teststat_Keil)$teststat,
      "P_value"            = padjusted_Keil,
      "Decision"           = decision_Keil,
      "Decision_overall"   = decision_overall_Keil,
      "df"                 = purrr::transpose(teststat_Keil)$df[[1]]   
    )
  }
  
# Information for Nitzl approach
  if(any(.approach_mgd %in% c("all", "Nitzl"))) {
    out[["Nitzl"]] <- list(
      "Test_statistic"     = purrr::transpose(teststat_Nitzl)$teststat,
      "P_value"            = padjusted_Nitzl,
      "Decision"           = decision_Nitzl,
      "Decision_overall"   = decision_overall_Nitzl,
      "df"                 = purrr::transpose(teststat_Nitzl)$df[[1]]
      
    )
  }

# Information for Henseler approach  
  if(any(.approach_mgd %in% c("all", "Henseler"))) {
    out[["Henseler"]] <- list(
      "Test_statistic"     = teststat_Henseler,
      "P_value"            = padjusted_Henseler,
      "Decision"           = decision_Henseler,
      "Decision_overall"   = decision_overall_Henseler
    )
  }
  
  #Information CI_param
  if(any(.approach_mgd %in% c("all", "CI_para"))) {
    out[["CI_para"]] <- list(
      "Test_statistic"     = NA,
      "Decision"           = decision_ci_para,
      "Decision_overall"   = decision_overall_ci_para
    )
  }
  
  if(any(.approach_mgd %in% c("all", "CI_overlap"))) {
    out[["CI_overlap"]] <- list(
      "Test_statistic"     = NA,
      "Decision"           = decision_ci_overlap,
      "Decision_overall"   = decision_overall_ci_overlap
    )
  }
   
  class(out) <- "cSEMTestMGD"
  return(out)
}
