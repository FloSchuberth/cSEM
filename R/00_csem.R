#' Composite-based SEM
#'
#' Estimate linear, nonlinear, hierachical or multigroup structural equation
#' models using a composite-based approach. In \pkg{cSEM} 
#' any method or approach that involves linear compounts (scores/proxies/composites)
#' of observables (indicators/items/manifest variables) is defined as composite-based.
#' See the \href{https://m-e-rademaker.github.io/cSEM/articles/cSEM.html}{Get started} 
#' section of the \href{https://m-e-rademaker.github.io/cSEM/index.html}{cSEM website}
#' for details.
#'
#' `csem()` estimates linear, nonlinear, hierachical or multigroup structural 
#' equation models using a composite-based approach. 
#' 
#' \subsection{Data and model:}{
#' The `.data` and `.model` arguments are required. Data must be
#' in a `matrix` or a `data.frame` with column names matching
#' the indicator names used in the model description of the measurement model.
#' Alternatively, a `list` of data sets (matrices or data frames) may be provided
#' in which case estimation is repeated for each data set. 
#' Possible column types/classes of the data provided are: "logical", 
#' "numeric" ("double" or "integer"), "factor" ("ordered" and/or "unordered"),
#' "character", or a mix of several types. Character columns will be treated 
#' as (undordered) factors.
#'
#' To provide a model use the [lavaan model syntax][lavaan::model.syntax].
#' Note, however, that \pkg{cSEM} currently only supports the "standard" lavaan
#' model syntax (Types 1, 2, 3, and 7 as described on the help page). 
#' Therefore, specifying e.g. a threshold or scaling factors is ignored. 
#' Alternatively a standardized (possibly incomplete) [cSEMModel]-list may be supplied.
#' See [parseModel()] for details.
#' }
#'
#' \subsection{Weights and path coefficients:}{
#' By default weights are estimated using the partial least squares (path) 
#' algorithm (`"PLS-PM"`).
#' A range of alternative weightning algorithms may be supplied to 
#' `.approach_weights`. Currently, the following approaches are implemented 
#' \enumerate{
#' \item{(Default) Partial least squares path modeling (`"PLS-PM"`). The algorithm
#'    can be customized. See [calculateWeightsPLS()] for details.}
#' \item{Generalized structured component analysis (`"GSCA"`). The algorithm 
#'   can be customized. See [calculateWeightsGSCA()] for details.}
#' \item{Generalized canonical correlation analysis (*GCCA*), including 
#'   `"SUMCORR"`, `"MAXVAR"`, `"SSQCORR"`, `"MINVAR"`, `"GENVAR"`.}
#' \item{Principal component analysis (`"PCA"`)}
#' \item{Factor score regression using sum scores (`"unit"`), 
#'    regression (`"regression"`) or bartlett scores (`"bartlett"`)}
#' }
#' 
#' Its possible to supply starting values for the weightning algorithm 
#' via `.starting_values`. The argument accepts a named list of vectors where the
#' list names are the construct names whose indicator weights the user
#' wishes to set. The vectors must be named vectors of `"indicator_name" = value` 
#' pairs, where `value` is the starting weight. See the examples section below for details.
#'
#' Composite-indicator and composite-composite correlations are properly
#' disattenuated by default to yield consistent loadings, construct correlations, 
#' and path coefficients if any of the concepts in the model are modeled as a 
#' common factor. 
#' 
#' For *PLS-PM* disattenuation is done using *PLSc* \insertCite{Dijkstra2015}{cSEM}.
#' For *GSCA* disattenuation is done implicitly by using *GSCAm* \insertCite{Hwang2017}{cSEM}. 
#' Weights obtained by *GCCA*, *unit*, *regression*, *bartlett* or *PCA* are 
#' disattenuated using Croon's approach \insertCite{Croon2002}{cSEM}.
#' Disattenuation my be suppressed by setting `.disattenuate = FALSE`. 
#' Note, however, that quantities in this case are inconsistent 
#' estimates for their construct level counterparts if any of the constructs in 
#' the structural model are modeled as a common factor!
#' 
#' By default. path coefficients are estimated using OLS (`.approach_path = "OLS"`). 
#' For linear models, two-stage least squares (`"2SLS"`) is available , however, *only if* 
#' *instruments are internal*, i.e. part of the structural model. Future versions
#' will add support for external instruments. Instruments must be supplied to 
#' `.instruments` as a named list where the names
#' of the list elements are the names of the dependent constructs of the structural 
#' equations whose explanatory variables are believed to be endogenous. 
#' The list consists of vectors of names of instruments corresponding to each equation. 
#' Note that exogenous variables of a given equation **must** be supplied as 
#'instruments for themselves. See the examples section.
#'
#' If reliabilities are known they can be supplied as `"name" = value` pairs to 
#' `.reliabilities`, where `value` is a numeric value between 0 and 1. 
#' Currently, only supported for "PLS-PM".
#' }
#'
#' \subsection{Nonlinear models:}{
#' If the model is nonlinear `csem()` estimates a polynomial structural equation model
#' using a non-iterative method of moments approach described in
#' \insertCite{Dijkstra2014}{cSEM}. Nonlinear terms include interactions and
#' exponential terms. The latter is described in model syntax as an
#' "interaction with itself", e.g., `x_1^3 = x1.x1.x1`. Currently only exponential
#' terms up to a power of three (i.e. three-way interactions) are allowed.
#'
#' The current version of the package allows two kinds of estimation:
#' estimation of the reduced form equation (`.approach_nl = "replace"`) and 
#' sequential estimation (`.approach_nl = "sequential"`, the default). The latter does not 
#' allow for multivariate normality of all exogenous variables, i.e., 
#' the latent variables and the error terms.
#'
#' Distributional assumptions are kept to a minimum (an i.i.d. sample from a 
#' population with finite moments for the relevant order); for higher order models, 
#' that go beyond interaction, we work in this version with the assumption that
#' as far as the relevant moments are concerned certain combinations of 
#' measurement errors behave as if they were Gaussian.
#' For details see: \insertCite{Dijkstra2014;textual}{cSEM}.
#' }
#' 
#' \subsection{Second-order model}{
#'  Second-order models are specified using the operators `=~` and `<~`. These
#'  operators are usually used with indicators on their right-hand side. For 
#'  second-order models the right-hand side variables are constructs instead:
#' \preformatted{my_model <- "
#' # Structural model
#' SAT ~ QUAL
#' VAL ~ SAT
#'
#' # Measurement model
#' EXPE <~ expe1 + expe2
#'
#' SAT =~ sat1 + sat2
#' VAL =~ val1 + val2
#' IMAG <~ imag1 + imag2
#' 
#' # Second-order term (in this case a second-order common factor)
#' QUAL =~ IMAG + EXPE
#' "
#' }
#'}
#' Currently, two approaches are available: `"repeated_indicators"` and 
#' `"2stage"` (the default). 
#' \subsection{Multigroup analysis}{
#' To perform multigroup analysis provide either a list of data sets or one 
#' data set containing a group-identifyer-column whose column 
#' name must be provided to `.id`. Values of this column are interpreted as group 
#' identifiers and `csem()` will split the data by levels of that column and run
#' the estimation for each level separately. Note that the more levels
#' the group-identifyer-column has, the more estimation runs are required.
#' This can considerably slow down estimation, especially if resampling is
#' requested. For the latter it will generally be faster to use 
#' `.eval_plan = "multiprocess"`.
#' } 
#' \subsection{Inference:}{
#' Inference is done via resampling. See [resamplecSEMResults] for details.
#' }
#' 
#' @usage csem(
#'   .data                    = NULL,
#'   .model                   = NULL,
#'   .approach_2ndorder       = c("3stage", "repeated_indicators"),
#'   .approach_nl             = c("sequential", "replace"),
#'   .approach_paths          = c("OLS", "2SLS", "3SLS"),
#'   .approach_weights        = c("PLS-PM", "SUMCORR", "MAXVAR", "SSQCORR", 
#'                                "MINVAR", "GENVAR", "GSCA", "PCA", 
#'                                "unit", "bartlett", "regression"),
#'   .disattenuate            = TRUE,
#'   .id                      = NULL,
#'   .instruments             = NULL,
#'   .normality               = FALSE,
#'   .reliabilities           = NULL,
#'   .starting_values         = NULL,
#'   .resample_method         = c("none", "bootstrap", "jackknife"),
#'   .resample_method2        = c("none", "bootstrap", "jackknife"),
#'   .R                       = 499,
#'   .R2                      = 199,
#'   .handle_inadmissibles    = c("drop", "ignore", "replace"),
#'   .user_funs               = NULL,
#'   .eval_plan               = c("sequential", "multiprocess"),
#'   .seed                    = NULL,
#'   .sign_change_option      = c("none", "individual", "individual_reestimate", 
#'                                "construct_reestimate"),
#'   ...
#'   )
#'
#' @param .data A `data.frame` or a `matrix` of standardized or unstandarized 
#'   data (indicators/items/manifest variables). 
#'   Additionally, a `list` of data sets (data frames or matrices) is accepted in which 
#'   case estimation is repeated for each data set. Possible column types or classes 
#'   of the data provided are: "logical", "numeric" ("double" or "integer"), 
#'   "factor" ("ordered" and/or "unordered"), "character" (converted to factor),
#'   or a mix of several types.
#' @inheritParams csem_arguments
#' @param ... Further arguments to be passed down to lower level functions of `csem()`.
#'   See [args_csem_dotdotdot] for a complete list of available arguments.
#'
#' @return
#' An object of class `cSEMResults` with methods for all postestimation generics.
#' Note that, technically, a call to [csem()] results in an object with at least 
#' two class attributes. The first class attribute is always `cSEMResults`. 
#' The second is one of `cSEMResults_default`, `cSEMResults_multi`, or 
#' `cSEMResults_2ndorder` and depends on the estimated model and/or the type of 
#' data provided to the `.model` and `.data` arguments. The third class attribute
#' `cSEMResults_resampled` is only added if resampling was conducted. 
#' Technically, method dispatch for all postestimation 
#' functions is based on the second class attribute. For a details see the 
#' [cSEMResults helpfile ][cSEMResults].
#'
#' @references
#'   \insertAllCited{}
#'
#' @seealso [args_default], [cSEMArguments], [cSEMResults], [foreman] 
#' 
#' @example inst/examples/example_csem.R
#' 
#' @export
#' 

csem <- function(
  .data                  = NULL,
  .model                 = NULL,
  .approach_2ndorder     = c("3stage", "repeated_indicators"),
  .approach_nl           = c("sequential", "replace"),
  .approach_paths        = c("OLS", "2SLS", "3SLS"),
  .approach_weights      = c("PLS-PM", "SUMCORR", "MAXVAR", "SSQCORR", 
                             "MINVAR", "GENVAR","GSCA", "PCA",
                             "unit", "bartlett", "regression"),
  .disattenuate          = TRUE,
  .id                    = NULL,
  .instruments           = NULL,
  .normality             = FALSE,
  .reliabilities         = NULL,
  .starting_values       = NULL,
  .resample_method       = c("none", "bootstrap", "jackknife"),
  .resample_method2      = c("none", "bootstrap", "jackknife"),
  .R                     = 499,
  .R2                    = 199,
  .handle_inadmissibles  = c("drop", "ignore", "replace"),
  .user_funs             = NULL,
  .eval_plan             = c("sequential", "multiprocess"),
  .seed                  = NULL,
  .sign_change_option    = c("none", "individual", "individual_reestimate", 
                             "construct_reestimate"),
  ...
  ) {
  ## Match arguments
  .approach_nl          <- match.arg(.approach_nl)
  .approach_2ndorder    <- match.arg(.approach_2ndorder)
  .approach_paths       <- match.arg(.approach_paths)
  .approach_weights     <- match.arg(.approach_weights)
  .resample_method      <- match.arg(.resample_method)
  .resample_method2     <- match.arg(.resample_method2)
  .handle_inadmissibles <- match.arg(.handle_inadmissibles)
  .eval_plan            <- match.arg(.eval_plan)
  .sign_change_option   <- match.arg(.sign_change_option)
  
  ## Collect and handle arguments
  # Note: all.names = TRUE is neccessary for otherwise arguments with a leading
  #       dot will be omitted!
  args_used <- c(as.list(environment(), all.names = TRUE), list(...))
  args_used["..."] <- NULL
  args        <- handleArgs(args_used)
  args_needed <- args[intersect(names(args), names(as.list(formals(foreman))))]
  
  ## Parse model
  model_original <- parseModel(.model, .instruments = .instruments)
  
  ## Modify model if model contains second order constructs
  if(any(model_original$construct_order == "Second order")) {
    model_1stage <- convertModel(
      .csem_model        = model_original, 
      .approach_2ndorder = args$.approach_2ndorder,
      .stage             = "first")
    ## Update model
    model_1stage$construct_order <- model_original$construct_order
    args_needed[[".model"]] <- model_1stage
  } else {
    args_needed[[".model"]] <- model_original
  }
    
  ## Check data
  if(!any(class(.data) %in% c("data.frame", "matrix", "list"))) {
    stop2(
      "The following error occured in the `csem()` function:\n",
      "Data must be provided as a `matrix`, a `data.frame` or a `list`. ", 
      ".data has class: ", 
      paste0(class(.data), collapse = ", ")
      )
  }
  ## Select cases
  if(!is.null(.id) && !inherits(.data, "list")) {

    if(length(.id) != 1) {
      stop2(
        "The following error occured in the `csem()` function:\n",
        "`.id` must be a character string or an integer identifying one single column."
        )
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
      ## NOTE: 
      #  01.03.2019: Using do.call(foreman, args_needed) would be more elegant but is much 
      #              much much! slower (especially for larger data sets).
      #
      #  15.05.2019: Apparently do.call(foreman, args_needed) is not that bad
      #              after all. I did several comparisons using microbenchmark
      #              but there was no speed difference (anymore?!). When I compared and
      #              thought that do.call is slow I fixed other parts of 
      #              foreman as well...maybe they had been the real culprit.
      #              So we use do.call again, as it avoids redundancy
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
  if(inherits(.data, "list") | !is.null(.id)) {
    
    # Sometimes the original (unstandardized) pooled data set is required 
    # (e.g. for permutation and permutation based tests).
    # By convention the original, pooled dataset is therefore added to the first 
    # element of "out" (= results for the first group/dataset)! 
    # If ".data" was a list of data they are combined to on pooled dataset
    # If ".data" was originally pooled and subsequently split by ".id"
    # The original unsplit data is returned.
    
    out[[1]]$Information$Data_pooled <- if(inherits(.data, "list")) {
      data_pooled <- do.call(rbind, .data)
      data_pooled <- as.data.frame(data_pooled)
      data_pooled[, "id"] <- rep(names(out), times = sapply(.data, nrow))
      data_pooled
    } else {
      .data
    }
    
    ## Set class
    class(out) <- c("cSEMResults", "cSEMResults_multi")

    ## Second order (2/3 stage approach)
    if(any(model_original$construct_order == "Second order") && 
       args$.approach_2ndorder == "3stage") {
      
      ### Second step
      out <- lapply(out, function(x) {
        out2 <- calculate2ndOrder(model_original, x)
        x <- list("First_stage" = x, "Second_stage" = out2)
        
        ## Append original arguments needed as they are required by e.g. testOMF.
        # Since 
        args_needed[[".model"]] <- model_original
        x$Second_stage$Information$Arguments_original <- args_needed
        
        # Return x
        x
      })
      class(out) <- c("cSEMResults", "cSEMResults_multi", "cSEMResults_2ndorder")
    }
    
  } else if(any(model_original$construct_order == "Second order") && 
            args$.approach_2ndorder == "3stage") {
    
    ### Second step
    # Note: currently data supplied as a list or grouped data is not allowed
    out2 <- calculate2ndOrder(model_original, out)
    out <- list("First_stage" = out, "Second_stage" = out2)
    
    ## Append original arguments needed as they are required by e.g. testOMF.
    # Since 
    args_needed[[".model"]] <- model_original
    out$Second_stage$Information$Arguments_original <- args_needed
    
    class(out) <- c("cSEMResults", "cSEMResults_2ndorder")
    
  } else {
    
    class(out) <- c("cSEMResults", "cSEMResults_default")
    
  }
  
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
      .eval_plan            = .eval_plan,
      .seed                 = .seed,
      .sign_change_option   = .sign_change_option,
      ...
    )
  }
  
  return(out)
}
