#' @usage cca(
#'   .data             = NULL,
#'   .model            = NULL,
#'   .id               = NULL,
#'   .approach_weights = c("PLS", "SUMCOR", "MAXVAR", "SSQCOR", "MINVAR", "GENVAR", "GSCA", 
#'                         "fixed", "unit"),
#'   ...) 
#'   
#' @name csem   
#' @rdname csem
#' @export

cca <- function(
  .data             = NULL, 
  .model            = NULL, 
  .id               = NULL,
  .approach_weights = c("PLS", "SUMCOR", "MAXVAR", "SSQCOR", "MINVAR", "GENVAR", "GSCA",
                        "fixed", "unit"),
  ...
  ) {
  
  ## Match arguments
  .approach_weights <- match.arg(.approach_weights)
  
  args_used <- c(as.list(environment(), all.names = TRUE), list(...))
  args_used["..."] <- NULL
  # Dont estimate the structural model
  args_used[[".estimate_structural"]] <- FALSE

  # For compatibility with csem()
  args_used[[".approach_paths"]] <- args_default()$.approach_paths
  
  # Handle arguments
  args        <- handleArgs(args_used)
  args_needed <- args[intersect(names(args), names(as.list(formals(foreman))))] 
  
  do.call(csem, args_needed)
}

