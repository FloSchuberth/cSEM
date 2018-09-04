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
  
  ## Structure output
  loading_estimates <- formatEstimates(.object$Estimates$Loading_estimates)
  weight_estimates  <- formatEstimates(.object$Estimates$Weight_estimates)
  path_estimates    <- formatEstimates(.object$Estimates$Path_estimates)
  
  ## Modify relevant .object elements
  
  .object$Estimates$Loading_estimates <- loading_estimates
  .object$Estimates$Weight_estimates  <- weight_estimates
  .object$Estimates$Path_estimates    <- path_estimates
  
  
  ## Set class for printing and return
  class(.object) <- "cSEMSummarize"
  return(.object)
}

formatEstimates <- function(x) {
  a <- paste0(rep(rownames(x), times = apply(x, 1, function(w) sum(w != 0))),
              " =~ ", colnames(x))
  data.frame("Name" = a,
             "Estimate" = unlist(t(x)[t(x) != 0 ]), stringsAsFactors = FALSE)
} 