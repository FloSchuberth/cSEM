#' Test for overall model fit
#' 
#' \lifecycle{maturing}
#'
#' Bootstrap-based test for overall model fit originally proposed by \insertCite{Beran1985;textual}{cSEM}. 
#' See also \insertCite{Dijkstra2015;textual}{cSEM} who first suggested the test in 
#' the context of PLS-PM.
#' 
#' By default, `testOMF()` tests the null hypothesis that the population indicator 
#' correlation matrix equals the population model-implied indicator correlation matrix. 
#' Several discrepancy measures may be used. By default, `testOMF()` uses four distance
#' measures to assess the distance between the sample indicator correlation matrix
#' and the estimated model-implied indicator correlation matrix, namely the geodesic distance, 
#' the squared Euclidean distance, the standardized root mean square residual (SRMR),
#' and the distance based on the maximum likelihood fit function. 
#' The reference distribution for each test statistic is obtained by 
#' the bootstrap as proposed by \insertCite{Beran1985;textual}{cSEM}. 
#' 
#' It is possible to perform the bootstrap-based test using fit measures such 
#' as the CFI, RMSEA or the GFI if `.fit_measures = TRUE`. This is experimental. 
#' To the best of our knowledge the applicability and usefulness of the fit 
#' measures for model fit assessment have not been formally (statistically) 
#' assessed yet. Theoretically, the logic of the test applies to these fit indices as well. 
#' Hence, their applicability is theoretically justified. 
#' Only use if you know what you are doing.
#' 
#' If `.saturated = TRUE` the original structural model is ignored and replaced by
#' a saturated model, i.e., a model in which all constructs are allowed to correlate freely. 
#' This is useful to test misspecification of the measurement model in isolation.
#' 
#' @usage testOMF(
#'  .object                = NULL, 
#'  .alpha                 = 0.05,
#'  .fit_measures          = FALSE,
#'  .handle_inadmissibles  = c("drop", "ignore", "replace"), 
#'  .R                     = 499, 
#'  .saturated             = FALSE,
#'  .seed                  = NULL,
#'  .verbose               = TRUE
#' )
#' 
#' @inheritParams  csem_arguments
#' 
#' @return
#' A list of class `cSEMTestOMF` containing the following list elements:
#' \describe{
#'   \item{`$Test_statistic`}{The value of the test statistics.}
#'   \item{`$Critical_value`}{The corresponding  critical values obtained by the bootstrap.}
#'   \item{`$Decision`}{The test decision. One of: `FALSE` (**Reject**) or `TRUE` (**Do not reject**).}
#'   \item{`$Information`}{The `.R` bootstrap values; The number of admissible results;
#'                         The seed used and the number of total runs.}
#' }
#' 
#' @seealso [csem()], [calculateSRMR()], [calculateDG()], [calculateDL()], [cSEMResults],
#'   [testMICOM()], [testMGD()]
#' 
#' @references
#'   \insertAllCited{}
#'   
#' @example inst/examples/example_testOMF.R
#' 
#' @export

testOMF <- function(
  .object                = NULL, 
  .alpha                 = 0.05, 
  .fit_measures          = FALSE,
  .handle_inadmissibles  = c("drop", "ignore", "replace"),
  .R                     = 499, 
  .saturated             = FALSE,
  .seed                  = NULL,
  .verbose               = TRUE
) {
  
  # Implementation is based on:
  # Dijkstra & Henseler (2015) - Consistent Paritial Least Squares Path Modeling
  
  ## Match arguments
  .handle_inadmissibles <- match.arg(.handle_inadmissibles)
  
  if(inherits(.object, "cSEMResults_default")) {
    
    x11 <- .object$Estimates
    x12 <- .object$Information
    
    ## Collect arguments
    arguments <- x12$Arguments
    
  } else if(inherits(.object, "cSEMResults_multi")) {
    
    out <- lapply(.object, testOMF,
                  .alpha                = .alpha,
                  .handle_inadmissibles = .handle_inadmissibles,
                  .R                    = .R,
                  .saturated            = .saturated,
                  .seed                 = .seed,
                  .verbose              = .verbose
    )
    ## Return
    return(out)
    
  } else if(inherits(.object, "cSEMResults_2ndorder")) { 
    
    x11 <- .object$First_stage$Estimates
    x12 <- .object$First_stage$Information
    
    x21 <- .object$Second_stage$Estimates
    x22 <- .object$Second_stage$Information
    
    ## Collect arguments
    arguments <- x22$Arguments_original
    
    ## Append the 2ndorder approach to args
    arguments$.approach_2ndorder <- x22$Approach_2ndorder
    
  } else {
    stop2(
      "The following error occured in the testOMF() function:\n",
      "`.object` must be a `cSEMResults` object."
    )
  }
  
  # if(.verbose) {
  #   cat(rule2("Test for overall model fit based on Beran & Srivastava (1985)",
  #             type = 3), "\n\n")
  # }
  
  ### Checks and errors ========================================================
  ## Check if initial results are inadmissible
  if(sum(unlist(verify(.object))) != 0) {
    stop2(
      "The following error occured in the `testOMF()` function:\n",
      "Initial estimation results are inadmissible. See `verify(.object)` for details.")
  }
  
  # Return error if used for ordinal variables
  if(any(x12$Type_of_indicator_correlation %in% c("Polyserial", "Polychoric"))){
    stop2(
      "The following error occured in the `testOMF()` function:\n",
      "Test for overall model fit currently not applicable if polychoric or",
      " polyserial indicator correlation is used.")
  }
  
  ### Preparation ==============================================================
  ## Extract required information 
  X         <- x12$Data
  S         <- x11$Indicator_VCV
  Sigma_hat <- fit(.object,
                   .saturated = .saturated,
                   .type_vcv  = "indicator")
  
  ## Calculate test statistic
  if(.fit_measures) {
    teststat <- c(
      "dG"           = calculateDG(.object),
      "SRMR"         = calculateSRMR(.object),
      "dL"           = calculateDL(.object),
      "dML"          = calculateDML(.object),
      "Chi_square"   = calculateChiSquare(.object),
      "Chi_square_df"= calculateChiSquareDf(.object),
      "CFI"          = calculateCFI(.object),
      "GFI"          = calculateGFI(.object),
      "IFI"          = calculateIFI(.object),
      "NFI"          = calculateNFI(.object),
      "NNFI"         = calculateNNFI(.object),
      "RMSEA"        = calculateRMSEA(.object)
      # "RMS_theta"    = calculateRMSTheta(.object, .model_implied = FALSE),
      # "RMS_theta_mi" = calculateRMSTheta(.object, .model_implied = TRUE)
    )
  } else {
    teststat <- c(
      "dG"           = calculateDG(.object),
      "SRMR"         = calculateSRMR(.object),
      "dL"           = calculateDL(.object),
      "dML"          = calculateDML(.object)
    )
  }

  
  ## Transform dataset, see Beran & Srivastava (1985)
  S_half   <- solve(expm::sqrtm(S))
  Sig_half <- expm::sqrtm(Sigma_hat)
  
  X_trans           <- X %*% S_half %*% Sig_half
  colnames(X_trans) <- colnames(X)
  
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
  
  ### Start resampling =========================================================
  ## Calculate reference distribution
  ref_dist         <- list()
  n_inadmissibles  <- 0
  counter <- 0
  progressr::with_progress({
    progress_bar_csem <- progressr::progressor(along = 1:.R)
    repeat{
      # Counter
      counter <- counter + 1
      progress_bar_csem(message = sprintf("Bootstrap run = %g", counter))
      
      # Draw dataset
      X_temp <- X_trans[sample(1:nrow(X), replace = TRUE), ]
      
      # Replace the old dataset by the new one
      arguments[[".data"]] <- X_temp
      
      # Estimate model
      Est_temp <- if(inherits(.object, "cSEMResults_2ndorder")) {
        
        do.call(csem, arguments) 
        
      } else {
        
        do.call(foreman, arguments)
      }           
      
      # Check status (Note: output of verify for second orders is a list)
      status_code <- sum(unlist(verify(Est_temp)))

      # Distinguish depending on how inadmissibles should be handled
      if(status_code == 0 | (status_code != 0 & .handle_inadmissibles == "ignore")) {
        # Compute if status is ok or .handle inadmissibles = "ignore" AND the status is 
        # not ok
        
        if(inherits(.object, "cSEMResults_default")) {
          S_temp         <- Est_temp$Estimates$Indicator_VCV
        } else if(inherits(.object, "cSEMResults_2ndorder")) { 
          S_temp         <- Est_temp$First_stage$Estimates$Indicator_VCV
        }
        
        Sigma_hat_temp <- fit(Est_temp,
                              .saturated = .saturated,
                              .type_vcv  = "indicator")
        
       ref_dist[[counter]] <- if(.fit_measures) {
        c(
          "dG"           = calculateDG(Est_temp),
          "SRMR"         = calculateSRMR(Est_temp),
          "dL"           = calculateDL(Est_temp),
          "dML"          = calculateDML(Est_temp),
          "Chi_square"   = calculateChiSquare(Est_temp),
          "Chi_square_df"= calculateChiSquareDf(Est_temp),
          "CFI"          = calculateCFI(Est_temp),
          "GFI"          = calculateGFI(Est_temp),
          "IFI"          = calculateIFI(Est_temp),
          "NFI"          = calculateNFI(Est_temp),
          "NNFI"         = calculateNNFI(Est_temp),
          "RMSEA"        = calculateRMSEA(Est_temp)
          # "RMS_theta"    = calculateRMSTheta(Est_temp, .model_implied = FALSE),
          # "RMS_theta_mi" = calculateRMSTheta(Est_temp, .model_implied = TRUE)
        ) 
      } else {
        c(
          "dG"           = calculateDG(Est_temp),
          "SRMR"         = calculateSRMR(Est_temp),
          "dL"           = calculateDL(Est_temp),
          "dML"          = calculateDML(Est_temp)
        ) 
      }
        
      } else if(status_code != 0 & .handle_inadmissibles == "drop") {
        # Set list element to zero if status is not okay and .handle_inadmissibles == "drop"
        ref_dist[[counter]] <- NA
        
      } else {# status is not ok and .handle_inadmissibles == "replace"
        # Reset counter and raise number of inadmissibles by 1
        counter <- counter - 1
        n_inadmissibles <- n_inadmissibles + 1
      }
      
      # Break repeat loop if .R results have been created.
      if(length(ref_dist) == .R) {
        break
      }
    } # END repeat 
  }) # END with_progress
  
  # Delete potential NA's
  ref_dist1 <- Filter(Negate(anyNA), ref_dist)
  
  # Check if at least 3 admissible results were obtained
  n_admissibles <- length(ref_dist1)
  if(n_admissibles < 3) {
    stop2("The following error occured in the `testOMF()` functions:\n",
          "Less than 2 admissible results produced.", 
          " Consider setting `.handle_inadmissibles == 'replace'` instead.")
  }
  
  # Combine
  ref_dist_matrix <- do.call(cbind, ref_dist1) 
  ## Compute critical values (Result is a (d x p) matrix, where p is the number
  ## of quantiles that have been computed (1 by default) and d the number of
  ## distance/fit measures
  .alpha <- .alpha[order(.alpha)]
  
  ## Goodness of fit measures require the 1 - alpha quantile
  if(.fit_measures) {
    ref_dist_matrix_good <- ref_dist_matrix[c("GFI", "CFI", "IFI", "NFI", "NNFI"), ,drop = FALSE]
    critical_values_good <- matrixStats::rowQuantiles(ref_dist_matrix_good, 
                                                      probs = .alpha, drop = FALSE)
    
    ref_dist_matrix_bad <- ref_dist_matrix[c("dG", "dL", "SRMR", "dML", 
                                             "Chi_square", "Chi_square_df"
                                             # "RMS_theta", "RMS_theta_mi"), 
                                             ), , drop = FALSE]
    ## Badness-of-fit indices require the alpha quantile
    critical_values_bad <- matrixStats::rowQuantiles(ref_dist_matrix_bad, 
                                                     probs =  1 - .alpha, drop = FALSE)
    
    ## Compare critical value and teststatistic.
    # Dont reject when teststat is >= alpha% quantile
    decision_good <- teststat[rownames(critical_values_good)] >= critical_values_good # a logical (d x p) matrix with each column
    # Dont reject when teststat is <= alpha% quantile
    decision_bad  <- teststat[rownames(critical_values_bad)] <= critical_values_bad # a logical (d x p) matrix with each column
    # representing the decision for one significance level. 
    # TRUE  = no evidence against the H0 --> not reject
    # FALSE = evidence against the H0    --> reject
    
    # Combine
    critical_values <- rbind(critical_values_bad, critical_values_good)
    decision        <- rbind(decision_bad, decision_good)
  } else {
    
    critical_values <- matrixStats::rowQuantiles(ref_dist_matrix, 
                                                 probs =  1-.alpha, drop = FALSE)
    
    ## Compare critical value and teststatistic
    decision <- teststat < critical_values # a logical (d x p) matrix with each column
    # representing the decision for one significance level. 
    # TRUE  = no evidence against the H0 --> not reject
    # FALSE = evidence against the H0    --> reject
  }
  
  # Return output
  out <- list(
    "Test_statistic"     = teststat[rownames(critical_values)],
    "Critical_value"     = critical_values, 
    "Decision"           = decision, 
    "Information"        = list(
      "Bootstrap_values"   = ref_dist,
      "Number_admissibles" = ncol(ref_dist_matrix),
      "Seed"               = .seed,
      "Total_runs"         = counter + n_inadmissibles
    )
  )
  
  ## Set class and return
  class(out) <- "cSEMTestOMF"
  return(out)
}