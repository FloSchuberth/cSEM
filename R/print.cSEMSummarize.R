#' `cSEMSummarize` method for `print()`
#'
#' The [cSEMSummary] method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cSEMResults], [summarize()]
#'
#' @export
#' @keywords internal
print.cSEMSummarize <- function(x, .full_output = TRUE, ...) {
  
  ## Check the class
  if(inherits(x, "cSEMSummarize_2ndorder")) {
    x11 <- x$First_stage$Estimates
    x12 <- x$First_stage$Information
    
    x21 <- x$Second_stage$Estimates
    x22 <- x$Second_stage$Information
    
    # Correlations
    construct_cor <- x21$Exo_construct_correlation
    res_cor  <- x11$Residual_correlation
    indi_cor <- x11$Indicator_correlation 
  } else {
    
    x21 <- x$Estimates
    x22 <- x$Information
    
    # Correlation
    construct_cor <- x21$Exo_construct_correlation
    res_cor  <- x21$Residual_correlation
    indi_cor <- x21$Indicator_correlation 
  }
  
  cat2(
    rule2(type = 2), "\n",
    rule2("Overview"), 
    "\n"
  )
  
  ### Overview -----------------------------------------------------------------
  ## General information + resample information
  printSummarizeOverview(x)
  
  ## Construct details
  cat2("\n\n\tConstruct details:\n\t","------------------")
  
  printSummarizeConstructDetails(x)
  
  ### Estimates ----------------------------------------------------------------
  cat2("\n\n", rule2("Estimates"), "\n\n")
  
  ## Confidence intervals
  # Get the column names of the columns containing confidence intervals
  ci_colnames <- colnames(x21$Path_estimates)[-c(1:6)]
  
  # Are there more confidence intervals than the default (the 95% percentile CI)
  # Inform the user to use xxx instead.
  if(length(ci_colnames) > 2) {
    cat2(
      "By default, only one confidence interval supplied to `.ci` is printed.\n",
      "Use `xxx` to print all confidence intervals (not yet implemented)."
    )
    ci_colnames <- ci_colnames[1:2]
    cat("\n\n")
  }
  
  ## Path estimates
  cat2("Estimated path coefficients:\n============================")
  printSummarizePathCorrelation(x, .ci_colnames = ci_colnames)
  
  ## Loadings and Weights
  printSummarizeLoadingsWeights(x, .ci_colnames = ci_colnames)
  
  ## Exogenous construct correlation
  if(.full_output && nrow(construct_cor) != 0) {
    cat2("\n\nEstimated construct correlations:\n=================================")
    printSummarizePathCorrelation(x, .ci_colnames = ci_colnames, 
                                  .what = "Construct correlation")
  }
  
  ## Residual correlation
  if(.full_output && nrow(res_cor) != 0) {
    cat2("\n\nEstimated measurement error correlations:\n=========================================")
    printSummarizePathCorrelation(x, .ci_colnames = ci_colnames, 
                                  .what = "Residual correlation")
  }
  
  ## Indicator correlation
  if(.full_output && nrow(indi_cor) != 0) {
    cat2("\n\nEstimated indicator correlations:\n=================================")
    printSummarizePathCorrelation(x, .ci_colnames = ci_colnames, 
                                  .what = "Indicator correlation")
  }
  
  if(.full_output && x22$Model$model_type == "Linear") {
    ### Effects ----------------------------------------------------------------
    cat2("\n\n", rule2("Effects"), "\n\n")
    ## Path estimates
    cat2("Estimated total effects:\n========================")
    
    printSummarizePathCorrelation(x, .ci_colnames = ci_colnames, 
                                  .what = "Total effect")
    
    cat2("\n\nEstimated indirect effects:\n===========================")
    
    printSummarizePathCorrelation(x, .ci_colnames = ci_colnames, 
                                  .what = "Indirect effect")
  }
  
  cat2("\n", rule2(type = 2))
}