#' Composite based SEM
#'
#' Estimates linear and nonlinear structural equation models using a
#' composite based approach.
#'
#' The `csem` function is a wrapper around the more general [foreman] function.
#' It is designed for quick, easy, and flexible use by providing the user with
#' defaults choices for all relevant arguments except the mandatory `.data`
#' and `.model` argument.
#'
#' @section Data and model:
#' To get started the `.data` and `.model` arguments are required. Data must be
#' provided as either a `matrix` or a `data.frame` with column names matching
#' the variable/indicator names used in the model description of the measurement model.
#' Alternativly a named or unamed list of matrices or `data.frame`s may be provided
#' in which case estimation is repeated for each data set. 
#' The data provided via `.data` may contain a non numeric column whose column name 
#' must be provided to `.id`. Values of this column are interpreted as group identifiers
#' and `csem()` will split the data by levels of that column and run the estimation by calling 
#' [foreman()] for each level separately.
#'
#' To provide a model use the \code{\link[lavaan:model.syntax]{lavaan model syntax}}
#' with two notable extensions/changes. First: the "`<~`" operator in `cSEM` is
#' used to define a composite instead of a formative common factor. Second:
#' the "`.`" is used to indicate interactions between constructs as in e.g.,
#' `construct1.construct2`.
#'
#' @section Weights:
#' By default weights are estimated using *PLS*. Alternative approaches include
#' all of *Kettenring's approaches*, "*fixed weights*"
#' or "*unit weight*". *Generalized Structured Component Analysis* (*GSCA*) may
#' also be chosen as a weighing approach although technically GSCA obtains weight
#' and structural coefficient estimates simultaneously. Hence, setting
#' `.approach_weights = "GSCA"` automatically sets `.approach_path = "GSCA"` (and
#' vice-versa).
#'
#' Composite-indicator and composite-composite correaltions are properly rescaled
#' to get consistent estimates for the factor loadings and construct correlations
#' unless `.disattenuate` is set explicitly set to `FALSE`.
#' Consistent estimates are calculated for the entries of the moment equations
#' that define the structural parameters.
#' The solutions to the moment equations are reported as consistent estimates
#' for the structural parameters.
#'
#' If the model is nonlinear [csem] estimates a polynomial structural equation model
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
#'
#' @usage csem(
#'   .data                = NULL,
#'   .model               = NULL,
#'   .approach_weights    = c("PLS", "SUMCOR", "MAXVAR", "SSQCOR", "MINVAR", "GENVAR", "GSCA", "fixed", "unit"),
#'   .approach_path       = c("OLS", "2SLS", "3SLS"),
#'   .approach_nl         = c("sequential", "replace"),
#'   .disattenuate        = TRUE,
#'   .PLS_weight_scheme_inner = c("centroid", "factorial", "path"),
#'   .PLS_mode            = NULL,
#'   .estimate_structural = TRUE,
#'   .reliabilities       = NULL
#'   ...)
#'
#' @param .data A `data.frame` or a `matrix` containing the raw data. Addionally,
#'   a list of `data.frame`s or `matrices` may be specified in which case estimation
#'   is repeated for each data set.
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
  .PLS_mode                = NULL,
  .estimate_structural     = TRUE,
  .reliabilities           = NULL,
  ...
  ) {
  
  ## Collect and handle arguments
  args_used   <- as.list(match.call())[-1]
  args        <- handleArgs(args_used)
  args_needed <- args[intersect(names(args), names(as.list(formals(foreman))))] 
  ## 
  
  if(!is.null(.id)) {
    
    data_split <- split(.data, f = .data[, .id])
    out <- lapply(data_split, function(x) {
      
      args_needed[[".data"]] <- x[, -which(names(x) == .id)]
      do.call(foreman, args_needed)
    })
  } else if(class(.data) == "list") {
    
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
