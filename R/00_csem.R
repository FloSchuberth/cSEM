#' Composite based SEM and CCA
#'
#' Estimate linear and nonlinear structural equation models using a 
#' composite based approach or conduct confirmatory composite analysis (CCA).
#'
#' `csem()` estimates linear and nonlinear structural equation models using a 
#' composite based approach like PLS-PM, GSCA or unit weights. Technically, `csem()` is a wrapper 
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
#' By default weights are estimated using the partial least squares algorithm (*PLS-PM*).
#' Alternative approaches include all of *Kettenring's criteria*, "*fixed weights*"
#' or "*unit weight*". *Generalized Structured Component Analysis* (*GSCA*) may
#' also be chosen as a weighing approach although technically GSCA obtains weight
#' and structural coefficient estimates simultaneously. Hence, setting
#' `.approach_weights = "GSCA"` automatically sets `.approach_paths = "GSCA"` (and
#' vice-versa).
#'
#' For PLS-PM composite-indicator and composite-composite correlations are properly
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
#'   .approach_2ndorder= c("3stage", "repeated_indicators"),
#'   .approach_weights = c("PLS-PM", "SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR", "GSCA", 
#'                         "fixed", "unit"),
#'   .approach_path    = c("OLS", "2SLS", "3SLS"),
#'   ...
#'   )
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
  .data                  = NULL,
  .model                 = NULL,
  .id                    = NULL,
  .approach_2ndorder     = c("3stage", "repeated_indicators"),
  .approach_paths        = c("OLS", "2SLS"),
  .approach_weights      = c("PLS-PM", "SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR",
                             "GSCA", "fixed", "unit"),
  .resample_method       = c("none", "bootstrap", "jackknife"),
  .resample_method2      = c("none", "bootstrap", "jackknife"),
  .R                     = 499,
  .R2                    = 199,
  .handle_inadmissibles  = c("drop", "ignore", "replace"),
  .user_funs             = NULL,
  .eval_plan             = c("sequential", "multiprocess"),
  ...
  ) {
  ## Match arguments
  .approach_2ndorder    <- match.arg(.approach_2ndorder)
  .approach_paths       <- match.arg(.approach_paths)
  .approach_weights     <- match.arg(.approach_weights)
  .resample_method      <- match.arg(.resample_method)
  .resample_method2     <- match.arg(.resample_method2)
  .handle_inadmissibles <- match.arg(.handle_inadmissibles)
  .eval_plan            <- match.arg(.eval_plan)
  
  ## Collect and handle arguments
  # Note: all.names = TRUE is neccessary for otherwise arguments with a leading
  #       dot will be omitted!
  args_used <- c(as.list(environment(), all.names = TRUE), list(...))
  args_used["..."] <- NULL
  args        <- handleArgs(args_used)
  args_needed <- args[intersect(names(args), names(as.list(formals(foreman))))]
  
  ## Parse model
  model <- parseModel(.model)
  
  ## Modify model if model contains second order constructs
  if(any(model$construct_order == "Second order")) {
    model1 <- convertModel(
      .csem_model        = model, 
      .approach_2ndorder = args$.approach_2ndorder,
      .stage             = "first")
    ## Update model
    model1$construct_order <- model$construct_order
    args_needed[[".model"]] <- model1
  } else {
    args_needed[[".model"]] <- model
  }
    
  ## Check data
  if(!any(class(.data) %in% c("data.frame", "matrix", "list"))) {
    stop("Data must be provided as a `matrix`, a `data.frame` or a `list`. ", 
         ".data has class: ", 
         class(.data), call. = FALSE)
  }
  ## Select cases
  if(!is.null(.id) && !inherits(.data, "list")) {

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
  # See the details section of ?UseMethod() to learn how method dispatch works
  # for objects with multiple classes
  if(any(class(.data) == "list") | !is.null(.id)) {
    
    # Sometimes the original (unstandardized) pooled data set is required 
    # (e.g. for permutation and permutation based tests).
    # By convention the original, pooled dataset is therefore added to the first 
    # element of "out" (= results for the first group/dataset)! 
    # If ".data" was a list of data they are combined to on pooled dataset
    # If ".data" was originallay pooled and subsequently split by ".id"
    # The original unsplit data is returned.
    
    out[[1]]$Information$Data_pooled <- if(any(class(.data) == "list")) {
      data_pooled <- do.call(rbind, .data)
      data_pooled[, "id"] <- rep(names(out), times = sapply(.data, nrow))
      data_pooled
    } else {
      .data
    }

    class(out) <- c("cSEMResults", "cSEMResults_multi")
    
    ## Resample if requested:
    if(.resample_method != "none") {
      out <- lapply(out, function(x) {
        resamplecSEMResults(
          .object               = x,
          .resample_method      = .resample_method,
          .resample_method2     = .resample_method2,
          .R                    = .R,
          .R2                   = .R2,
          .handle_inadmissibles = .handle_inadmissibles,
          .user_funs            = .user_funs,
          .eval_plan            = .eval_plan  
        )
      })
    }
    
  } else if(any(model$construct_order == "Second order") && 
            args$.approach_2ndorder == "3stage") {
    
    ### Second step
    # Note: currently only data supplied as a list or grouped data is not allowed
    out2 <- calculate2ndOrder(model, out)
    out <- list("First_stage" = out, "Second_stage" = out2)
    ## Append original arguments needed as they are required by e.g. testOMF.
    
    out$Second_stage$Information$Arguments_original <- args_needed
    
    class(out) <- c("cSEMResults", "cSEMResults_2ndorder")
    
  } else {
    
    class(out) <- c("cSEMResults", "cSEMResults_default")
    
    ## Resample if requested:
    if(.resample_method != "none") {
      out <- resamplecSEMResults(
        .object               = out,
        .resample_method      = .resample_method,
        .resample_method2     = .resample_method2,
        .R                    = .R,
        .R2                   = .R2,
        .handle_inadmissibles = .handle_inadmissibles,
        .user_funs            = .user_funs,
        .eval_plan            = .eval_plan  
      )
    }
  }
  
  return(out)
}
