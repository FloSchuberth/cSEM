#' Composite-based SEM
#'
#' Estimate linear, nonlinear, hierachical or multigroup structural equation
#' models using a composite-based approach. In \pkg{cSEM} 
#' any method or approach that involves linear compounts (scores/proxies/composites)
#' of observables (indicators/items/manifest variables) is defined as composite-based.
#' See the [the cSEM vignette](../doc/vignette-cSEM.html) for details.
#' 
#'
#' `csem()` estimates linear, nonlinear, hierachical or multigroup structural 
#' equation models using a composite-based approach. 
#' 
#' \subsection{Data and model:}{
#' The `.data` and `.model` arguments are required. Data must be
#' provided as either a `matrix` or a `data.frame` with column names matching
#' the indicator names used in the model description of the measurement model.
#' Alternatively, a list of matrices or `data.frame`'s may be provided
#' in which case estimation is repeated for each data set. 
#' The data provided via `.data` may contain **one** character column whose column name 
#' must be provided to `.id`. Values of this column are interpreted as group 
#' identifiers and `csem()` will split the data by levels of that column and run
#' the estimation for each level separately.
#'
#' To provide a model use the [lavaan model syntax][lavaan::model.syntax] 
#' with two notable extensions/changes. First: the "`<~`" operator in `cSEM` is
#' used to define a composite instead of a causal formative common factor. Second:
#' the "`.`" is used to indicate interactions between constructs as in e.g.,
#' `construct1.construct2`. Alternatively a standardized (possibly incomplete)
#' [cSEMModel]-list may be supplied.
#' }
#'
#' \subsection{Weights and path coefficients:}{
#' By default weights are estimated using the partial least squares (path) algorithm (`"PLS-PM"`).
#' A broad range of alternative weightning algorithms may be supplied to `.appraoch_weights`.
#' Currently the following approaches are implemented 
#' \enumerate{
#' \item{(Default) Partial least squares path modeling (`"PLS"`). The algorithm
#'    can be customized. See [calculateWeightsPLS()] for details.}
#' \item{Generalized structured component analysis (`"GSCA"`)}
#' \item{Generalized canoncial correlation analysis (*GCCA*), including 
#'   `"SUMCORR"`, `"MAXVAR"`, `"SSQCORR"`, `"MINVAR"`, `"GENVAR"`}
#' \item{Principal component analysis (`"PCA"`)}
#' \item{Factor score regression using sum scores (`"unit"`), 
#'    regression (`"regression"`) or bartlett scores (`"bartlett"`)}
#' }
#'
#' Composite-indicator and composite-composite correlations are properly
#' disattenuated by default to yield consistent loadings, construct correlations, 
#' and path coefficients if any of the constructs in the model are modeled as a 
#' common factor. 
#' 
#' For *PLS-PM* disattenuation is done using *PLSc* \insertCite{Dijkstra2015}{cSEM}.
#' For *GSCA* disattenuation is done implicitly by using *GSCAm*. Weights obtained
#' by GCCA, unit weights and fixed weights are disattenuated using a  
#' Disattenuation my be suppressed by setting 
#' `.disattenuate = FALSE`. Note, however, that quantities in this case are inconsistent 
#' estimates for their construct level counterparts if any of the constructs in the structural model is
#' modeled as a common factor!
#' }
#'
#' \subsection{Nonlinear models:}{
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
#' For details see: \insertCite{Dijkstra2014;textual}{cSEM}.
#' }
#' 
#' \subsection{Inference:}{
#' Inference is done via resampling. See [resamplecSEMResults] for details.
#' }
#' 
#' @usage csem(
#'   .data                        = NULL,
#'   .model                       = NULL,
#'   .approach_nl                 = c("sequential", "replace"),
#'   .approach_paths              = c("OLS", "2SLS"),
#'   .approach_weights            = c("PLS-PM", "SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR", 
#'                                    "GSCA", "PCA", "unit", "bartlett", "regression"),
#'   .disattenuate                = TRUE,
#'   .id                          = NULL,
#'   .normality                   = TRUE,
#'   .reliabilities               = NULL,
#'   .resample_method             = c("none", "bootstrap", "jackknife"),
#'   .resample_method2            = c("none", "bootstrap", "jackknife"),
#'   .R                           = 499,
#'   .R2                          = 199,
#'   .handle_inadmissibles        = c("drop", "ignore", "replace"),
#'   .user_funs                   = NULL,
#'   .eval_plan                   = c("sequential", "multiprocess"),
#'   .seed                        = NULL,
#'   ...
#'   )
#'
#' @param .data A `data.frame` or a `matrix` of standardized or unstandarized 
#'   data (indicators/items/manifest variables). 
#'   Additionally, a `list` of `data.frame`(s) or `matrice`(s) is accepted in which 
#'   case estimation is repeated for each data set. Possible column types or classes 
#'   of the data provided are: "logical", "numeric" ("double" or "integer"), 
#'   "factor" ("ordered" and/or "unordered") or a mix of several types. 
#'   The data may also include *one* character column whose column name must be 
#'   given to `.id`. This column is assumed to contain group identifiers used to split 
#'   the data into groups.
#' @inheritParams csem_arguments
#' @param ... Further arguments to be passed down to lower level functions of `csem()`.
#'   See [csem_dotdotdot] for a complete list of available arguments.
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
#' @examples
#' 
#' @export
#' 

csem <- function(
  .data                  = NULL,
  .model                 = NULL,
  .approach_2ndorder     = c("3stage", "repeated_indicators"),
  .approach_nl           = c("sequential", "replace"),
  .approach_paths        = c("OLS", "2SLS"),
  .approach_weights      = c("PLS-PM", "SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR",
                             "GSCA", "unit", "PCA", "bartlett", "regression"),
  .disattenuate          = TRUE,
  .id                    = NULL,
  .normality             = TRUE,
  .reliabilities         = NULL,
  .resample_method       = c("none", "bootstrap", "jackknife"),
  .resample_method2      = c("none", "bootstrap", "jackknife"),
  .R                     = 499,
  .R2                    = 199,
  .handle_inadmissibles  = c("drop", "ignore", "replace"),
  .user_funs             = NULL,
  .eval_plan             = c("sequential", "multiprocess"),
  .seed                  = NULL,
  .sign_change_option    = c("no",'individual','individual_reestimate','construct_reestimate'),
  .starting_values       = NULL,
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
  model_original <- parseModel(.model)
  
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
      ## NOTE: using do.call(foreman, args_needed) would be more elegant but is much 
      # much much! slower (especially for larger data sets). There is now a
      # terible redundancy because foreman is repeated 3 times. This should be
      # addressed at some point, however, it is much better than do.call in terms
      # of speed
      # out <- do.call(foreman, args_needed)
      foreman(
        .data                        = args_needed$.data,
        .model                       = args_needed$.model,
        .approach_cor_robust         = args_needed$.approach_cor_robust,
        .approach_nl                 = args_needed$.approach_nl,
        .approach_paths              = args_needed$.approach_paths,
        .approach_weights            = args_needed$.approach_weights,
        .conv_criterion              = args_needed$.conv_criterion,
        .disattenuate                = args_needed$.disattenuate,
        .dominant_indicators         = args_needed$.dominant_indicators,
        .estimate_structural         = args_needed$.estimate_structural,
        .id                          = args_needed$.id,
        .iter_max                    = args_needed$.iter_max,
        .normality                   = args_needed$.normality,
        .PLS_approach_cf             = args_needed$.PLS_approach_cf,
        .PLS_ignore_structural_model = args_needed$.PLS_ignore_structural_model,
        .PLS_modes                   = args_needed$.PLS_modes,
        .PLS_weight_scheme_inner     = args_needed$.PLS_weight_scheme_inner,
        .reliabilities               = args_needed$.reliabilities,
        .tolerance                   = args_needed$.tolerance
      )
    })
  } else if(any(class(.data) == "list")) {
    
    out <- lapply(.data, function(x) {
      
      args_needed[[".data"]] <- x
      ## NOTE: using do.call(foreman, args_needed) would be more elegant but is much 
      # much much! slower (especially for larger data sets). 
      # out <- do.call(foreman, args_needed)
      foreman(
        .data                        = args_needed$.data,
        .model                       = args_needed$.model,
        .approach_cor_robust         = args_needed$.approach_cor_robust,
        .approach_nl                 = args_needed$.approach_nl,
        .approach_paths              = args_needed$.approach_paths,
        .approach_weights            = args_needed$.approach_weights,
        .conv_criterion              = args_needed$.conv_criterion,
        .disattenuate                = args_needed$.disattenuate,
        .dominant_indicators         = args_needed$.dominant_indicators,
        .estimate_structural         = args_needed$.estimate_structural,
        .id                          = args_needed$.id,
        .iter_max                    = args_needed$.iter_max,
        .normality                   = args_needed$.normality,
        .PLS_approach_cf             = args_needed$.PLS_approach_cf,
        .PLS_ignore_structural_model = args_needed$.PLS_ignore_structural_model,
        .PLS_modes                   = args_needed$.PLS_modes,
        .PLS_weight_scheme_inner     = args_needed$.PLS_weight_scheme_inner,
        .reliabilities               = args_needed$.reliabilities,
        .tolerance                   = args_needed$.tolerance
      )
    })
    if(is.null(names(.data))) {
      names(out) <- paste0("Data_", 1:length(out))
    } else {
      names(out) <- names(.data)
    }
  } else {
    ## NOTE: using do.call(foreman, args_needed) would be more elegant but is much 
    # much much! slower (especially for larger data sets). 
    # out <- do.call(foreman, args_needed)
    out <- foreman(
      .data                        = args_needed$.data,
      .model                       = args_needed$.model,
      .approach_cor_robust         = args_needed$.approach_cor_robust,
      .approach_nl                 = args_needed$.approach_nl,
      .approach_paths              = args_needed$.approach_paths,
      .approach_weights            = args_needed$.approach_weights,
      .conv_criterion              = args_needed$.conv_criterion,
      .disattenuate                = args_needed$.disattenuate,
      .dominant_indicators         = args_needed$.dominant_indicators,
      .estimate_structural         = args_needed$.estimate_structural,
      .id                          = args_needed$.id,
      .iter_max                    = args_needed$.iter_max,
      .normality                   = args_needed$.normality,
      .PLS_approach_cf             = args_needed$.PLS_approach_cf,
      .PLS_ignore_structural_model = args_needed$.PLS_ignore_structural_model,
      .PLS_modes                   = args_needed$.PLS_modes,
      .PLS_weight_scheme_inner     = args_needed$.PLS_weight_scheme_inner,
      .reliabilities               = args_needed$.reliabilities,
      .tolerance                   = args_needed$.tolerance
    )
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
    # If ".data" was originally pooled and subsequently split by ".id"
    # The original unsplit data is returned.
    
    out[[1]]$Information$Data_pooled <- if(any(class(.data) == "list")) {
      data_pooled <- do.call(rbind, .data)
      data_pooled <- as.data.frame(data_pooled)
      data_pooled[, "id"] <- rep(names(out), times = sapply(.data, nrow))
      data_pooled
    } else {
      .data
    }

    class(out) <- c("cSEMResults", "cSEMResults_multi")
    
  } else if(any(model_original$construct_order == "Second order") && 
            args$.approach_2ndorder == "3stage") {
    
    ### Second step
    # Note: currently only data supplied as a list or grouped data is not allowed
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
    
    if(is.null(.seed)) {
      .seed <- sample(.Random.seed, 1)
    }
    
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
      .sign_change_option   = .sign_change_option
    )
  }
  
  return(out)
}
