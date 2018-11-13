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
summarize <- function(.object, .inference = "Bootstrap") {
  UseMethod("summarize")
}

#' @describeIn summarize (TODO)
#' @export

summarize.cSEMResults_default <- function(.object, .inference = "Bootstrap") {

  if(.inference == "Bootstrap") {
    boot_out <- infer(.object, 
                      .method = .inference, 
                      .runs = 100, 
                      .handle_inadmissibles = "drop",
                      .verbose = TRUE)
  }
  x1  <- .object$Estimates
  x2  <- .object$Information

  ### Structure output =========================================================
  ## Path estimates ------------------------------------------------------------
  # Get construct type for relevant variables
  type <- rep(x2$Model$construct_type, times = rowSums(x2$Model$structural))
  
  # Build names
  temp <- outer(rownames(x1$Path_estimates), colnames(x1$Path_estimates), 
                FUN = function(x, y) paste(x, y, sep = " ~ "))
  
  temp_boot <- boot_out$Path_estimates
  t_temp    <- t(x1$Path_estimates)[t(x2$Model$structural) != 0 ] / temp_boot$Sd
  
  path_estimates <- data.frame(
    "Name"           = t(temp)[t(x2$Model$structural) != 0],
    "Construct_type" = type,
    "Estimate"       = t(x1$Path_estimates)[t(x2$Model$structural) != 0 ],
    "Std_Err"        = temp_boot$Sd,
    "t_statistic"    = t_temp,
    "p_value"        = 2*pnorm(abs(t_temp), lower.tail = FALSE),
    stringsAsFactors = FALSE)
  
  ## Loading estimates ---------------------------------------------------------
  # Get construct type for relevant variables
  type <- rep(x2$Model$construct_type, times = rowSums(x2$Model$measurement))
  
  # Build names
  temp <- rep(rownames(x1$Loading_estimates), times = rowSums(x2$Model$measurement))
  temp <- paste0(temp, 
                 ifelse(x2$Model$construct_type[temp] == "Composite", " <~ ", " =~ "),  
                 colnames(x1$Loading_estimates))
  
  temp_boot <- boot_out$Loading_estimates
  t_temp    <- t(x1$Loading_estimates)[t(x2$Model$measurement) != 0 ] / temp_boot$Sd
  
  loading_estimates <- data.frame(
    "Name"           = temp,
    "Construct_type" = type,
    "Estimate"       = x1$Loading_estimates[x2$Model$measurement != 0 ], 
    "Std_Err"        = temp_boot$Sd,
    "t_statistic"    = t_temp,
    "p_value"        = 2*pnorm(abs(t_temp), lower.tail = FALSE),
    stringsAsFactors = FALSE)
  
  ## Weight estimates ----------------------------------------------------------
  temp <- rep(rownames(x1$Weight_estimates), times = rowSums(x2$Model$measurement))
  temp <- paste0(temp, " -- ", colnames(x1$Weight_estimates))
  
  weight_estimates <- data.frame(
    "Name"           = temp,
    "Construct_type" = type, 
    "Estimate"       = x1$Weight_estimates[x2$Model$measurement != 0 ], 
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
