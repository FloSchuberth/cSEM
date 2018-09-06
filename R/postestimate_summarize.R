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
summarize <- function(.object) {
  
  if(attr(.object, "single") == FALSE) {
    stop("`summarize()` not applicable to multiple groups or data sets.\n",
         "Use lapply(.object, summarize) instead.",
         call. = FALSE)
  }
  
  x1  <- .object$Estimates
  x2  <- .object$Information
  
  ### Structure output
  ## Path coefficients
  temp <- outer(rownames(x1$Path_estimates), colnames(x1$Path_estimates), 
                FUN = function(x, y) paste(x, y, sep = " ~ "))
  
  path_estimates <- data.frame(
    "Name" = temp[x2$Model$structural != 0],
    "Estimate" = x1$Path_estimates[x2$Model$structural != 0 ], 
    stringsAsFactors = FALSE)
  
  ## Loading estimates
  temp <- rep(rownames(x1$Loading_estimates), 
              times = apply(x1$Loading_estimates, 1, function(w) sum(w != 0)))
  temp <- paste0(temp, 
                 ifelse(x2$Model$construct_type[temp] == "Composite", " <~ ", " =~ "),  
                 colnames(x1$Loading_estimates))
  
  loading_estimates <- data.frame(
    "Name" = temp,
    "Estimate" = x1$Loading_estimates[x2$Model$measurement != 0 ], 
                               stringsAsFactors = FALSE)
  
  ## Loading estimates
  temp <- rep(rownames(x1$Weight_estimates), 
              times = apply(x1$Weight_estimates, 1, function(w) sum(w != 0)))
  temp <- paste0(temp, " -- ", colnames(x1$Weight_estimates))
  
  weight_estimates <- data.frame(
    "Name" = temp,
    "Estimate" = x1$Weight_estimates[x2$Model$measurement != 0 ], 
    stringsAsFactors = FALSE)
  
  ## Modify relevant .object elements
  
  .object$Estimates$Loading_estimates <- loading_estimates
  .object$Estimates$Weight_estimates  <- weight_estimates
  .object$Estimates$Path_estimates    <- path_estimates
  
  
  ## Set class for printing and return
  class(.object) <- "cSEMSummarize"
  return(.object)
}
