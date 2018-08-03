#' Composite based SEM
#'
#' Estimates linear and nonlinear structural equation models using a
#' composite based approach.
#'
#' The `csem()` function is a wrapper around the more general [foreman()] function.
#' It is designed for quick, easy, and flexible use by providing the user with
#' defaults choices for all relevant arguments except the mandatory `.data`
#' and `.model` argument.
#'
#' \subsection{Data and model:}{
#' The `.data` and `.model` arguments are required. Data must be
#' provided as either a `matrix` or a `data.frame` with column names matching
#' the indicator names used in the model description of the measurement model.
#' Alternativly a named or unamed list of matrices or `data.frame`s may be provided
#' in which case estimation is repeated for each data set. 
#' The data provided via `.data` may contain a non numeric column whose column name 
#' must be provided to `.id`. Values of this column are interpreted as group 
#' identifiers and `csem()` will split the data by levels of that column and run
#' the estimation by calling [foreman()] for each level separately.
#'
#' To provide a model use the \code{\link[lavaan:model.syntax]{lavaan model syntax}}
#' with two notable extensions/changes. First: the "`<~`" operator in `cSEM` is
#' used to define a composite instead of a formative common factor. Second:
#' the "`.`" is used to indicate interactions between constructs as in e.g.,
#' `construct1.construct2`. Alternativly a standardized (possibly incomplete)
#' [cSEMModel]-list may be supplied.
#' }
#'
#' \subsection{Weights and path coefficients:}{
#' By default weights are estimated using the partial least squares algorithm (*PLS*).
#' Alternative approaches include all of *Kettenring's criteria*, "*fixed weights*"
#' or "*unit weight*". *Generalized Structured Component Analysis* (*GSCA*) may
#' also be chosen as a weighing approach although technically GSCA obtains weight
#' and structural coefficient estimates simultaneously. Hence, setting
#' `.approach_weights = "GSCA"` automatically sets `.approach_path = "GSCA"` (and
#' vice-versa).
#'
#' For PLS composite-indicator and composite-composite correlations are properly
#' rescaled using *PLSc* by yielding consistent estimates for the factor loadings, construct 
#' correlations, and path coefficients if any of the constructs involved is 
#' modelled as a common factor. If no disattenuation should be done, 
#' set `.disattenuate = FALSE`. 
#' }
#'
#' \subsection{Non linear models:}{
#' If the model is nonlinear `csem()` estimates a polynomial structural equation model
#' using a non-iterative method of moments approach described in
#' \insertCite{Dijkstra2014}{cSEM}. Non linear terms include interactions and
#' exponential terms. The latter is described in model syntax as an
#' "interaction with itself", e.g., `x_1^3 = x1.x1.x1`. Currently only exponential
#' terms up to a power of three (i.e. three-way interactions) are allowed.
#'
#' The current version of the package allows two kinds of estimation:
#' estimation of the reduced form equation (`.approach_nl = "reduced"`) and 
#' sequential estimation (`.approach_nl = "sequential"`). The latter does not 
#' not allow for multivariate normality of all exogenous variables, i.e., 
#' the latent variables and the error terms.
#'
#' Distributional assumptions are kept to a minimum (an i.i.d. sample from a 
#' population with finite moments for the relevant order); for higher order models, 
#' that go beyond interaction, we work in this version with the assumption that
#' as far as the relevant moments are concerned certain combinations of 
#' measurement errors behave as if they were Gaussian.
#' }
#' @usage csem(
#'   .data                = NULL,
#'   .model               = NULL,
#'   .id                  = NULL,
#'   .approach_weights    = c("PLS", "SUMCOR", "MAXVAR", "SSQCOR", "MINVAR", "GENVAR", "GSCA", "fixed", "unit"),
#'   .approach_path       = c("OLS", "2SLS", "3SLS"),
#'   .approach_nl         = c("sequential", "replace"),
#'   .disattenuate        = TRUE,
#'   .PLS_weight_scheme_inner = c("centroid", "factorial", "path"),
#'   .PLS_modes           = NULL,
#'   .estimate_structural = TRUE,
#'   .reliabilities       = NULL
#'   ...)
#'
#' @param .data A `data.frame` or a `matrix` containing the raw data. Addionally,
#'   a list of `data.frame`s or `matrices` may be specified in which case estimation
#'   is repeated for each data set. See details.
#' @inheritParams csem_arguments
#'
#' @inherit foreman return
#'
#' @references
#'   \insertAllCited{}
#'
#' @seealso [cca], [foreman]
#'
#' @examples
#' \dontrun{
#' # still to implement
#' }
#'
#' @export

csem <- function(
  .data                    = NULL,
  .model                   = NULL,
  .id                      = NULL,
  .approach_weights        = c("PLS", "SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR",
                               "GSCA", "fixed", "unit"),
  .approach_paths          = c("OLS", "2SLS", "3SLS"),
  .approach_nl             = c("sequential", "replace"),
  .disattenuate            = TRUE,
  .PLS_weight_scheme_inner = c("centroid", "factorial", "path"),
  .PLS_modes               = NULL,
  .estimate_structural     = TRUE,
  .reliabilities           = NULL,
  ...
  ) {

  ## Match arguments (for meaningful error messages)
  .approach_weights        <- match.arg(.approach_weights)
  .approach_paths          <- match.arg(.approach_paths)
  .approach_nl             <- match.arg(.approach_nl)
  .PLS_weight_scheme_inner <- match.arg(.PLS_weight_scheme_inner)

  ## Collect and handle arguments
  # Note: all.names = TRUE is neccessary for otherwise arguments with a leading
  #       dot will be omitted!
  args_used <- c(as.list(environment(), all.names = TRUE), list(...))
  args_used["..."] <- NULL
  args        <- handleArgs(args_used)
  args_needed <- args[intersect(names(args), names(as.list(formals(foreman))))] 
  
  ## Select cases
  if(!is.null(.id)) {
    if(is.matrix(.data)) {
      .data <- as.data.frame(.data)
    }
    if(length(.id) != 1) {
      stop("`.id` must be a character string or an integer identifying one single column.",
           call. = FALSE)
    }
    
    data_split <- split(.data, f = .data[, .id])
    out <- lapply(data_split, function(x) {
      if(is.numeric(.id)) {
        args_needed[[".data"]] <- x[, -.id]
      } else {
        args_needed[[".data"]] <- x[, -which(names(x) == .id)]
      }
      do.call(foreman, args_needed)
      
    })
  } else if(any(class(.data) == "list")) {
    
    out <- lapply(.data, function(x) {
      
      args_needed[[".data"]] <- x
      do.call(foreman, args_needed)
    })
    if(is.null(names(.data))) {
      names(out) <- paste0("Data_", 1:length(out))
    } else {
      names(out) <- names(.data)
    }
  } else {
    
    out <- do.call(foreman, args_needed)
  }
  
  ## Set class for output
  class(out) <- "cSEMResults"
  return(out)
}
