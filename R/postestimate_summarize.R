#' Summarize model
#'
#' Summarize the model (TODO). 
#' 
#' Summary (TODO)
#'
#' @usage summarize(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @export
#'
summarize <- function(
  .object                = NULL, 
  .csem_resample         = args_default()$.csem_resample,
  .alpha                 = args_default()$.alpha,
  .bias_corrected        = args_default()$.bias_corrected,
  .statistic             = args_default()$.statistic,
  .verbose               = args_default()$.verbose
  ) {
  UseMethod("summarize")
}

#' @describeIn summarize (TODO)
#' @export

summarize.cSEMResults_default <- function(
  .object                = NULL, 
  .csem_resample         = args_default()$.csem_resample,
  .alpha                 = args_default()$.alpha,
  .bias_corrected        = args_default()$.bias_corrected,
  .statistic             = args_default()$.statistic,
  .verbose               = args_default()$.verbose
  ) {
  
  x1  <- .object$Estimates
  x2  <- .object$Information

  ### Structure output =========================================================
  ## Path estimates ------------------------------------------------------------
  # Get construct type for relevant variables
  type <- rep(x2$Model$construct_type, times = rowSums(x2$Model$structural))
  
  # Build names
  temp <- outer(rownames(x1$Path_estimates), colnames(x1$Path_estimates), 
                FUN = function(x, y) paste(x, y, sep = " ~ "))
  
  path_estimates <- data.frame(
    "Name"           = t(temp)[t(x2$Model$structural) != 0],
    "Construct_type" = type,
    "Estimate"       = t(x1$Path_estimates)[t(x2$Model$structural) != 0 ],
    "Std_err"        = NA,
    "t_stat"         = NA,
    "p_value"        = NA,
    stringsAsFactors = FALSE)
  
  ## Loading estimates ---------------------------------------------------------
  # Get construct type for relevant variables
  type <- rep(x2$Model$construct_type, times = rowSums(x2$Model$measurement))
  
  # Build names
  temp <- rep(rownames(x1$Loading_estimates), times = rowSums(x2$Model$measurement))
  temp <- paste0(temp, " =~ ", colnames(x1$Loading_estimates))
  
  loading_estimates <- data.frame(
    "Name"           = temp,
    "Construct_type" = type,
    "Estimate"       = x1$Loading_estimates[x2$Model$measurement != 0 ],
    "Std_err"        = NA,
    "t_stat"         = NA,
    "p_value"        = NA,
    stringsAsFactors = FALSE)
  
  ## Weight estimates ----------------------------------------------------------
  temp <- rep(rownames(x1$Weight_estimates), times = rowSums(x2$Model$measurement))
  temp <- paste0(temp, " <~ ", colnames(x1$Weight_estimates))
  
  weight_estimates <- data.frame(
    "Name"           = temp,
    "Construct_type" = type, 
    "Estimate"       = x1$Weight_estimates[x2$Model$measurement != 0 ],
    "Std_err"        = NA,
    "t_stat"         = NA,
    "p_value"        = NA,
    stringsAsFactors = FALSE)
  
  ## Inner weight estimates ----------------------------------------------------
  if(x2$Arguments$.approach_weights == "PLS-PM") {
    i <- rownames(x1$Inner_weight_estimates)
    D <- x2$Model$structural[i, i , drop = FALSE] + t(x2$Model$structural[i, i, drop = FALSE])
    
    temp <- outer(i, i, FUN = function(x, y) paste(x, y, sep = " -- "))
    type <- rep(x2$Model$construct_type, times = colSums(D))
    
    inner_weight_estimates <- data.frame(
      "Name"           = t(temp)[t(D) != 0],
      "Construct_type" = type, 
      "Estimate"       = t(x1$Inner_weight_estimates)[t(D) != 0 ], 
      stringsAsFactors = FALSE)
  }

  ## Construct scores ----------------------------------------------------------
  construct_scores <- as.data.frame(x1$Construct_scores)

  ## Inference =================================================================
  
  if(!is.null(.csem_resample)) {
    if(class(.csem_resample) != "cSEMResults_resampled") {
      stop2(".csem_results must be an object of class `cSEMResults_resampled`.",
            " See ?resamplecSEMResults for details.")
    }
  
    infer_out <- infer(
      .csem_resample  = .csem_resample,
      .alpha          = .alpha,
      .bias_corrected = .bias_corrected,
      .statistic      = .statistic
    )
  
  temp   <- infer_out$Path_estimates
  t_temp <- t(x1$Path_estimates)[t(x2$Model$structural) != 0 ] / temp$Sd
  
  path_estimates["Std_err"] <- temp$Sd
  path_estimates["t_stat"]  <- t_temp
  path_estimates["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
  
  temp   <- infer_out$Loading_estimates
  t_temp <- t(x1$Loading_estimates)[t(x2$Model$measurement) != 0 ] / temp$Sd
  
  loading_estimates["Std_err"] <- temp$Sd
  loading_estimates["t_stat"]  <- t_temp
  loading_estimates["p_value"] <- 2*pnorm(abs(t_temp), lower.tail = FALSE)
  
  ## CI
  CIs <- c("CI_standard_z", "CI_standard_t", "CI_percentile", "CI_t_intervall",
           "CI_Bc", "CI_Bca")
  
  # Which statistics should be reported
  stats <- intersect(.statistic, CIs)
  # if(anyNA(x2) & any(.statistic %in% c("all", "CI_t_intervall", "CI_Bca"))) {
  #   stop2(paste0("`", intersect(.statistic, c("all", "CI_t_intervall", "CI_Bca")), 
  #                "`", collapse = ", "), 
  #         " requires the resampling distribution for each resample.")
  # }
  }
  
  ## Modify relevant .object elements ------------------------------------------
  .object$Estimates$Path_estimates    <- path_estimates
  .object$Estimates$Loading_estimates <- loading_estimates
  .object$Estimates$Weight_estimates  <- weight_estimates
  .object$Estimates$Inner_weight_estimates <- inner_weight_estimates
  .object$Estimates$Construct_scores <- construct_scores
  
  
  ## Set class for printing and return
  class(.object) <- "cSEMSummarize_default"
  return(.object)
}

#' @describeIn summarize (TODO)
#' @export

summarize.cSEMResults_multi <- function(.object) {
  
 lapply(.object, summarize.cSEMResults_default)
}

#' @describeIn summarize (TODO)
#' @export

summarize.cSEMResults_2ndorder <- function(.object) {
  
  ## Run summarize for each stage
  x <- lapply(.object, summarize.cSEMResults_default)
  
  x21 <- x$Second_stage$Estimates
  x22 <- x$Second_stage$Information
  
  ### Modify second stage estimates ============================================
  ## Path estimates 
  # Only second order path model estimates are relevant. Delete the "_temp" 
  # suffix.
  x21$Path_estimates$Name <- gsub("_temp", "", x21$Path_estimates$Name)

  ## Loading estimates 
  # Loadings for first stage (=all loadings except those of the 2nd order constructs)
  # Loadings for second stage (all loadings for 2nd order constructs)
  i <- x21$Loading_estimates$Name[!grepl("_temp", x21$Loading_estimates$Name)]
  x21$Loading_estimates <- x21$Loading_estimates[x21$Loading_estimates$Name %in% i, , drop = FALSE]
  
  ## Weights estimates 
  # Wieghts for first stage (=all weights except those of the 2nd order constructs)
  # Weights for second stage (all weights for 2nd order constructs)
  i <- x21$Weight_estimates$Name[!grepl("_temp", x21$Weight_estimates$Name)]
  x21$Weight_estimates <- x21$Weight_estimates[x21$Weight_estimates$Name %in% i, , drop = FALSE]
  
  ## Inner weight estimates
  # Rename. Delete the "_temp" suffix.
  # x21$Inner_weight_estimates$Name <- gsub("_temp", "", x21$Inner_weight_estimates$Name)
  
  ## Construct scores
  # colnames(x21$Construct_scores) <- gsub("_temp", "", colnames(x21$Construct_scores))
  
  ## Set class for printing and return
  out <- list("First_stage" = x$First_stage, 
              "Second_stage" = list("Estimates" = x21, "Information" = x22))
  
  class(out) <- "cSEMSummarize_2ndorder"
  return(out)
}
