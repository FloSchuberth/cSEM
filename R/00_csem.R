#' Composite based SEM and CCA
#'
#' Estimate linear and nonlinear structural equation models using a 
#' composite based approach or conduct confirmatory composite analysis (CCA).
#'
#' `csem()` estimates linear and nonlinear structural equation models using a 
#' composite based approach like PLS, GSCA or unit weights. Technically, `csem()` is a wrapper 
#' around the more general [foreman()] function designed for quick and flexible 
#' use by providing the user with default options except for 
#' the mandatory `.data` and `.model` argument. 
#' 
#' `cca()` performs CCA. CCA and SEM differ in that the former allows all 
#' constructs to vary freely. Hence, `cca()` is technically a simple convenience wrapper
#' around `csem(..., .estimate_structural = FALSE)`.
#'
#' \subsection{Data and model:}{
#' The `.data` and `.model` arguments are required. Data must be
#' provided as either a `matrix` or a `data.frame` with column names matching
#' the indicator names used in the model description of the measurement model.
#' Alternatively, a list of matrices or `data.frame`s may be provided
#' in which case estimation is repeated for each data set. 
#' The data provided via `.data` may contain a non numeric column whose column name 
#' must be provided to `.id`. Values of this column are interpreted as group 
#' identifiers and `csem()` or `cca()` will split the data by levels of that column and run
#' the estimation for each level separately.
#'
#' To provide a model use the \code{\link[lavaan:model.syntax]{lavaan model syntax}}
#' with two notable extensions/changes. First: the "`<~`" operator in `cSEM` is
#' used to define a composite instead of a formative common factor. Second:
#' the "`.`" is used to indicate interactions between constructs as in e.g.,
#' `construct1.construct2`. Alternatively a standardized (possibly incomplete)
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
#' rescaled using *PLSc* \insertCite{Dijkstra2015}{cSEM} by default. *PLSc* yields
#' consistent estimates for the factor loadings, construct correlations, 
#' and path coefficients if any of the constructs involved is 
#' modeled as a common factor. Disattenuation my be suppressed by setting 
#' `.disattenuate = FALSE`. Note however that quantities in this case are inconsistent 
#' estimates for their construct level counterparts if any of the constructs is
#' the model as a common factor.
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
#'   .data             = NULL,
#'   .model            = NULL,
#'   .id               = NULL,
#'   .approach_weights = c("PLS", "SUMCOR", "MAXVAR", "SSQCOR", "MINVAR", "GENVAR", "GSCA", 
#'                         "fixed", "unit"),
#'   .approach_path    = c("OLS", "2SLS", "3SLS"),
#'   ...)
#'
#' @param .data A `data.frame` or a `matrix` containing the raw data. Additionally,
#'   a list of `data.frame`(s) or `matrices` is accepted in which case estimation
#'   is repeated for each data set. See details.
#' @inheritParams csem_arguments
#' @param ... Further arguments to be passed down to lower level functions of `csem()`
#'   or `cca()`. Type \code{\link[cSEM:args_default]{args_default(.only_dots = TRUE)}} 
#'   or \code{\link[cSEM:args_default]{args_default(.only_dots = TRUE, .which_fun = "cca")}}
#'   for a complete list of accepted `...` arguments for the respective function.
#'
#' @inherit csem_results return
#'
#' @references
#'   \insertAllCited{}
#'
#' @seealso [args_default], [cSEMArguments], [cSEMResults], [foreman] 
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#' 
#' @export

csem <- function(
  .data                    = NULL,
  .model                   = NULL,
  .id                      = NULL,
  .approach_weights        = c("PLS", "SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR",
                               "GSCA", "fixed", "unit"),
  .approach_paths          = c("OLS", "2SLS"),
  ...
  ) {
  ## Match arguments
  .approach_weights <- match.arg(.approach_weights)
  .approach_paths    <- match.arg(.approach_paths)
  
  ## Collect and handle arguments
  # Note: all.names = TRUE is neccessary for otherwise arguments with a leading
  #       dot will be omitted!
  args_used <- c(as.list(environment(), all.names = TRUE), list(...))
  args_used["..."] <- NULL
  args        <- handleArgs(args_used)
  args_needed <- args[intersect(names(args), names(as.list(formals(foreman))))] 
  
  ## Check data
  if(!any(class(.data) %in% c("data.frame", "matrix", "list"))) {
    stop("Data must be provided as a `matrix`, a `data.frame` or a `list`. ", ".data has class: ", 
         class(.data), call. = FALSE)
  }
  ## Select cases
  if(!is.null(.id) && !is.list(.data)) {

    if(length(.id) != 1) {
      stop("`.id` must be a character string or an integer identifying one single column.",
           call. = FALSE)
    }
    
    if(is.matrix(.data)) {
      .data <- as.data.frame(.data)
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
  ## Has estimation been done based on a single data 
  attr(out, "single") <- ifelse(any(class(.data) == "list")| !is.null(.id), FALSE, TRUE) 
  return(out)
}
