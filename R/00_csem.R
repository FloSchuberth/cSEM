#' Composite-based SEM
#'
#'\lifecycle{stable}
#'
#' Estimate linear, nonlinear, hierarchical or multigroup structural equation
#' models using a composite-based approach. In \pkg{cSEM} 
#' any method or approach that involves linear compounds (scores/proxies/composites)
#' of observables (indicators/items/manifest variables) is defined as composite-based.
#' See the \href{https://m-e-rademaker.github.io/cSEM/articles/cSEM.html}{Get started} 
#' section of the \href{https://m-e-rademaker.github.io/cSEM/index.html}{cSEM website}
#' for a general introduction to composite-based SEM and \pkg{cSEM}.
#'
#' `csem()` estimates linear, nonlinear, hierarchical  or multigroup structural 
#' equation models using a composite-based approach. 
#' 
#' \subsection{Data and model:}{
#' The `.data` and `.model` arguments are required. `.data` must be given
#' a `matrix` or a `data.frame` with column names matching
#' the indicator names used in the model description. Alternatively, 
#' a `list` of data sets (matrices or data frames) may be provided
#' in which case estimation is repeated for each data set. 
#' Possible column types/classes of the data provided are: "`logical`", 
#' "`numeric`" ("`double`" or "`integer`"), "`factor`" ("`ordered`" and/or "`unordered`"),
#' "`character`", or a mix of several types. Character columns will be treated 
#' as (unordered) factors.
#' 
#' Depending on the type/class of the indicator data provided cSEM computes the indicator 
#' correlation matrix in different ways. See [calculateIndicatorCor()] for details.
#'
#' In the current version `.data` must not contain missing values. Future versions
#' are likely to handle missing values as well.
#' 
#' To provide a model use the [lavaan model syntax][lavaan::model.syntax].
#' Note, however, that \pkg{cSEM} currently only supports the "standard" lavaan
#' model syntax (Types 1, 2, 3, and 7 as described on the help page). 
#' Therefore, specifying e.g., a threshold or scaling factors is ignored. 
#' Alternatively, a standardized (possibly incomplete) [cSEMModel]-list may be supplied.
#' See [parseModel()] for details.
#' }
#'
#' \subsection{Weights and path coefficients:}{
#' By default weights are estimated using the partial least squares path modeling 
#' algorithm (`"PLS-PM"`).
#' A range of alternative weighting algorithms may be supplied to 
#' `.approach_weights`. Currently, the following approaches are implemented 
#' \enumerate{
#' \item{(Default) Partial least squares path modeling (`"PLS-PM"`). The algorithm
#'    can be customized. See [calculateWeightsPLS()] for details.}
#' \item{Generalized structured component analysis (`"GSCA"`) and generalized 
#'   structured component analysis with uniqueness terms (GSCAm). The algorithms
#'   can be customized. See [calculateWeightsGSCA()] and [calculateWeightsGSCAm()] for details.
#'   Note that GSCAm is called indirectly when the model contains constructs
#'   modeled as common factors only and `.disattenuate = TRUE`. See below.}
#' \item{Generalized canonical correlation analysis (*GCCA*), including 
#'   `"SUMCORR"`, `"MAXVAR"`, `"SSQCORR"`, `"MINVAR"`, `"GENVAR"`.}
#' \item{Principal component analysis (`"PCA"`)}
#' \item{Factor score regression using sum scores (`"unit"`), 
#'    regression (`"regression"`) or bartlett scores (`"bartlett"`)}
#' }
#' 
#' It is possible to supply starting values for the weighting algorithm 
#' via `.starting_values`. The argument accepts a named list of vectors where the
#' list names are the construct names whose indicator weights the user
#' wishes to set. The vectors must be named vectors of `"indicator_name" = value` 
#' pairs, where `value` is the starting weight. See the examples section below for details.
#'
#' Composite-indicator and composite-composite correlations are properly
#' disattenuated by default to yield consistent loadings, construct correlations, 
#' and path coefficients if any of the concepts are modeled as a 
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
#' By default path coefficients are estimated using ordinary least squares (`.approach_path = "OLS"`). 
#' For linear models, two-stage least squares (`"2SLS"`) is available, however, *only if* 
#' *instruments are internal*, i.e., part of the structural model. Future versions
#' will add support for external instruments if possible. Instruments must be supplied to 
#' `.instruments` as a named list where the names
#' of the list elements are the names of the dependent constructs of the structural 
#' equations whose explanatory variables are believed to be endogenous. 
#' The list consists of vectors of names of instruments corresponding to each equation. 
#' Note that exogenous variables of a given equation **must** be supplied as 
#'instruments for themselves.
#'
#' If reliabilities are known they can be supplied as `"name" = value` pairs to 
#' `.reliabilities`, where `value` is a numeric value between 0 and 1. 
#' Currently, only supported for "PLS-PM".
#' }
#'
#' \subsection{Nonlinear models:}{
#' If the model contains nonlinear terms `csem()` estimates a polynomial structural equation model
#' using a non-iterative method of moments approach described in
#' \insertCite{Dijkstra2014;textual}{cSEM}. Nonlinear terms include interactions and
#' exponential terms. The latter is described in model syntax as an
#' "interaction with itself", e.g., `xi^3 = xi.xi.xi`. Currently only exponential
#' terms up to a power of three (e.g., three-way interactions or cubic terms) are allowed:
#' \enumerate{
#' \item{- Single, e.g., `eta1`}
#' \item{- Quadratic, e.g., `eta1.eta1`}
#' \item{- Cubic, e.g., `eta1.eta1.eta1`}
#' \item{- Two-way interaction, e.g., `eta1.eta2`}
#' \item{- Three-way interaction, e.g., `eta1.eta2.eta3`}
#' \item{- Quadratic and two-way interaction, e.g., `eta1.eta1.eta3`}
#' }
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
#' \subsection{Models containing second-order constructs}{
#' Second-order constructs are specified using the operators `=~` and `<~`. These
#' operators are usually used with indicators on their right-hand side. For 
#' second-order constructs the right-hand side variables are constructs instead.
#' If c1, and c2 are constructs forming or measuring a higher-order
#' construct, a model would look like this:
#' \preformatted{my_model <- "
#' # Structural model
#' SAT  ~ QUAL
#' VAL  ~ SAT
#'
#' # Measurement/composite model
#' QUAL =~ qual1 + qual2
#' SAT  =~ sat1 + sat2
#' 
#' c1 =~ x11 + x12
#' c2 =~ x21 + x22
#' 
#' # Second-order construct (in this case a second-order composite build by common
#' # factors)
#' VAL <~ c1 + c2
#' "
#' }
#' Currently, two approaches are explicitly implemented: 
#' \itemize{
#' \item{(Default) `"2stage"`. The (disjoint) two-stage approach as proposed by \insertCite{Agarwal2000;textual}{cSEM}.
#' Note that by default a correction for attenuation is applied if common factors are 
#' involved in modeling second-order constructs. For instance, the three-stage approach
#'  proposed by \insertCite{VanRiel2017;textual}{cSEM} is applied in case of a second-order construct specified as a 
#'  composite of common factors. On the other hand, if no common factors are involved the two-stage approach 
#'  is applied as proposed by \insertCite{Schuberth2020;textual}{cSEM}.}
#' \item{`"mixed"`. The mixed repeated indicators/two-stage approach as proposed by \insertCite{Ringle2012;textual}{cSEM}.}
#' }
#'  
#' 
#' The repeated indicators approach as proposed by \insertCite{Joereskog1982b;textual}{cSEM}
#' and the extension proposed by \insertCite{Becker2012;textual}{cSEM} are 
#' not directly implemented as they simply require a respecification  of the model. 
#' In the above example the repeated indicators approach
#' would require to change the model and to append the repeated indicators to 
#' the data supplied to `.data`. Note that the indicators need to be renamed in this case as 
#' `csem()` does not allow for one indicator to be attached to multiple constructs.
#' \preformatted{my_model <- "
#' # Structural model
#' SAT  ~ QUAL
#' VAL  ~ SAT
#' 
#' VAL ~ c1 + c2
#'
#' # Measurement/composite model
#' QUAL =~ qual1 + qual2
#' SAT  =~ sat1 + sat2
#' VAL  =~ x11_temp + x12_temp + x21_temp + x22_temp
#' 
#' c1 =~ x11 + x12
#' c2 =~ x21 + x22
#' "
#' }
#' According to the extended approach indirect effects of `QUAL` on `VAL` via `c1` 
#' and `c2` would have to be specified as well.
#'}
#' 
#' \subsection{Multigroup analysis}{
#' To perform a multigroup analysis provide either a list of data sets or one 
#' data set containing a group-identifier-column whose column 
#' name must be provided to `.id`. Values of this column are taken as levels of a
#' factor and are interpreted as group 
#' identifiers. `csem()` will split the data by levels of that column and run
#' the estimation for each level separately. Note that the more levels
#' the group-identifier-column has, the more estimation runs are required.
#' This can considerably slow down estimation, especially if resampling is
#' requested. For the latter it will generally be faster to use 
#' `.eval_plan = "multisession"` or `.eval_plan = "multicore"`.
#' } 
#' \subsection{Inference:}{
#' Inference is done via resampling. See [resamplecSEMResults()] and [infer()] for details.
#' }
#' 
#' @usage csem(
#' .data                  = NULL,
#' .model                 = NULL,
#' .approach_2ndorder     = c("2stage", "mixed"),
#' .approach_cor_robust   = c("none", "mcd", "spearman"),
#' .approach_nl           = c("sequential", "replace"),
#' .approach_paths        = c("OLS", "2SLS"),
#' .approach_weights      = c("PLS-PM", "SUMCORR", "MAXVAR", "SSQCORR", 
#'                            "MINVAR", "GENVAR","GSCA", "PCA",
#'                            "unit", "bartlett", "regression"),
#' .conv_criterion        = c("diff_absolute", "diff_squared", "diff_relative"),
#' .disattenuate          = TRUE,
#' .dominant_indicators   = NULL,
#' .estimate_structural   = TRUE,
#' .id                    = NULL,
#' .instruments           = NULL,
#' .iter_max              = 100,
#' .normality             = FALSE,
#' .PLS_approach_cf       = c("dist_squared_euclid", "dist_euclid_weighted", 
#'                            "fisher_transformed", "mean_arithmetic",
#'                            "mean_geometric", "mean_harmonic",
#'                            "geo_of_harmonic"),
#' .PLS_ignore_structural_model = FALSE,
#' .PLS_modes                   = NULL,
#' .PLS_weight_scheme_inner     = c("path", "centroid", "factorial"),
#' .reliabilities         = NULL,
#' .starting_values       = NULL,
#' .resample_method       = c("none", "bootstrap", "jackknife"),
#' .resample_method2      = c("none", "bootstrap", "jackknife"),
#' .R                     = 499,
#' .R2                    = 199,
#' .handle_inadmissibles  = c("drop", "ignore", "replace"),
#' .user_funs             = NULL,
#' .eval_plan             = c("sequential", "multicore", "multisession"),
#' .seed                  = NULL,
#' .sign_change_option    = c("none", "individual", "individual_reestimate", 
#'                            "construct_reestimate"),
#' .tolerance             = 1e-05
#' )
#'
#' @param .data A `data.frame` or a `matrix` of standardized or unstandardized  
#'   data (indicators/items/manifest variables). 
#'   Additionally, a `list` of data sets (data frames or matrices) is accepted in which 
#'   case estimation is repeated for each data set. Possible column types or classes 
#'   of the data provided are: "`logical`", "`numeric`" ("`double`" or "`integer`"), 
#'   "`factor`" ("`ordered`" and/or "`unordered`"), "`character`" (will be converted to factor),
#'   or a mix of several types.
#' @inheritParams csem_arguments
#'
#' @return
#' An object of class `cSEMResults` with methods for all postestimation generics.
#' Technically, a call to [csem()] results in an object with at least 
#' two class attributes. The first class attribute is always `cSEMResults`. 
#' The second is one of `cSEMResults_default`, `cSEMResults_multi`, or 
#' `cSEMResults_2ndorder` and depends on the estimated model and/or the type of 
#' data provided to the `.model` and `.data` arguments. The third class attribute
#' `cSEMResults_resampled` is only added if resampling was conducted. 
#' For a details see the [cSEMResults helpfile ][cSEMResults].
#'
#' @section Postestimation:
#' \describe{
#' \item{[assess()]}{Assess results using common quality criteria, e.g., reliability,
#'   fit measures, HTMT, R2 etc.}
#' \item{[infer()]}{Calculate common inferential quantities, e.g., standard errors, 
#'   confidence intervals.}
#' \item{[predict()]}{Predict endogenous indicator scores and compute common prediction metrics.}
#' \item{[summarize()]}{Summarize the results. Mainly called for its side-effect the print method.}
#' \item{[verify()]}{Verify/Check admissibility of the estimates.}
#' }
#' 
#' Tests are performed using the test-family of functions. Currently the following
#' tests are implemented:
#' \describe{
#' \item{[testOMF()]}{Bootstrap-based test for overall model fit based on 
#'   \insertCite{Beran1985;textual}{cSEM}}
#' \item{[testMICOM()]}{Permutation-based test for measurement invariance of composites
#' proposed by \insertCite{Henseler2016;textual}{cSEM}}
#' \item{[testMGD()]}{Several (mainly) permutation-based tests for multi-group comparisons.}
#' \item{[testHausman()]}{Regression-based Hausman test to test for endogeneity.}
#' }
#' 
#' Other miscellaneous postestimation functions belong do the do-family of functions.
#' Currently three do functions are implemented:
#' \describe{
#' \item{[doIPMA()]}{Performs an importance-performance matrix analyis (IPMA).}
#' \item{[doNonlinearEffectsAnalysis()]}{Perform a nonlinear effects analysis as
#'   described in e.g.,
#'   \insertCite{Spiller2013;textual}{cSEM}}
#' \item{[doRedundancyAnalysis()]}{Perform a redundancy analysis (RA) as proposed by 
#' \insertCite{Hair2016;textual}{cSEM} with reference to \insertCite{Chin1998;textual}{cSEM}}
#' }
#' 
#' @references
#'   \insertAllCited{}
#'
#' @seealso [args_default()], [cSEMArguments], [cSEMResults], [foreman()], [resamplecSEMResults()],
#'   [assess()], [infer()], [predict()], [summarize()], [verify()], [testOMF()],
#'   [testMGD()], [testMICOM()], [testHausman()]
#' 
#' @example inst/examples/example_csem.R
#' 
#' @export

csem <- function(
  .data                  = NULL,
  .model                 = NULL,
  .approach_2ndorder     = c("2stage", "mixed"),
  .approach_cor_robust   = c("none", "mcd", "spearman"),
  .approach_nl           = c("sequential", "replace"),
  .approach_paths        = c("OLS", "2SLS"),
  .approach_weights      = c("PLS-PM", "SUMCORR", "MAXVAR", "SSQCORR", 
                             "MINVAR", "GENVAR","GSCA", "PCA",
                             "unit", "bartlett", "regression"),
  .conv_criterion        = c("diff_absolute", "diff_squared", "diff_relative"),
  .disattenuate          = TRUE,
  .dominant_indicators   = NULL,
  .estimate_structural   = TRUE,
  .id                    = NULL,
  .instruments           = NULL,
  .iter_max              = 100,
  .normality             = FALSE,
  .PLS_approach_cf       = c("dist_squared_euclid", "dist_euclid_weighted", 
                             "fisher_transformed", "mean_arithmetic",
                             "mean_geometric", "mean_harmonic",
                             "geo_of_harmonic"),
  .PLS_ignore_structural_model = FALSE,
  .PLS_modes                   = NULL,
  .PLS_weight_scheme_inner     = c("path", "centroid", "factorial"),
  .reliabilities         = NULL,
  .starting_values       = NULL,
  .resample_method       = c("none", "bootstrap", "jackknife"),
  .resample_method2      = c("none", "bootstrap", "jackknife"),
  .R                     = 499,
  .R2                    = 199,
  .handle_inadmissibles  = c("drop", "ignore", "replace"),
  .user_funs             = NULL,
  .eval_plan             = c("sequential", "multicore", "multisession"),
  .seed                  = NULL,
  .sign_change_option    = c("none", "individual", "individual_reestimate", 
                             "construct_reestimate"),
  .tolerance             = 1e-05
  ) {
  ## Match arguments
  .approach_2ndorder    <- match.arg(.approach_2ndorder)
  .approach_cor_robust  <- match.arg(.approach_cor_robust)
  .approach_nl          <- match.arg(.approach_nl)
  .approach_paths       <- match.arg(.approach_paths)
  .approach_weights     <- match.arg(.approach_weights)
  .conv_criterion       <- match.arg(.conv_criterion)
  .eval_plan            <- match.arg(.eval_plan)
  .handle_inadmissibles <- match.arg(.handle_inadmissibles)
  .PLS_approach_cf      <- match.arg(.PLS_approach_cf)
  .PLS_weight_scheme_inner <- match.arg(.PLS_weight_scheme_inner)
  .resample_method      <- match.arg(.resample_method)
  .resample_method2     <- match.arg(.resample_method2)
  .sign_change_option   <- match.arg(.sign_change_option)
  
  ## Collect and handle arguments
  # Note: all.names = TRUE is neccessary for otherwise arguments with a leading
  #       dot will be omitted!
  args_used <- c(as.list(environment(), all.names = TRUE))
  args        <- handleArgs(args_used)
  args_needed <- args[intersect(names(args), names(as.list(formals(foreman))))]
  
  ## Check if columnames of .data contain a "." and stop if so
  ## Check data
  if(!any(class(.data) %in% c("data.frame", "matrix", "list"))) {
    stop2(
      "The following error occured in the `csem()` function:\n",
      "Data must be provided as a `matrix`, a `data.frame` or a `list`. ", 
      ".data has class: ", 
      paste0(class(.data), collapse = ", ")
    )
  }
  if(inherits(.data,  "list")) {
    c_names <- unique(unlist(lapply(.data, colnames)))
  } else {
    c_names <- colnames(.data)
  }
  
  if(length(grep("\\.", c_names)) > 0) {
   stop2(
   "At least one variable name in your data set contain a `.` (dot).",
   " Dots are a reserved special character in cSEM. Please rename these variables in your data and the model description.") 
  }
  
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
    # If ".data" was a list of data they are combined to one pooled dataset
    # If ".data" was originally pooled and subsequently split by ".id"
    # The original unsplit data is returned.
    
    out[[1]]$Information$Data_pooled <- if(inherits(.data, "list")) {
      # Clean and order columns according to the first data set
      data_cleaned <- lapply(.data, function(x) {
        # Order data according to the ordering of the measurement model; delete
        # all columns that are not needed
        x <- x[, setdiff(colnames(model_original$measurement), model_original$vars_attached_to_2nd)]
        x
      })
      
      # Combine
      data_pooled <- do.call(rbind, data_cleaned)
      data_pooled <- as.data.frame(data_pooled)
      # data_pooled[, "id"] <- rep(names(out), times = sapply(.data, nrow))
      data_pooled
    } else {
      .data
    }
    ## Add second order approach to $Information
    if(any(model_original$construct_order == "Second order")) {
      out <- lapply(out, function(x){
        x$Information$Approach_2ndorder <- .approach_2ndorder
        x$Information$Model_original    <- model_original
        x
      })
    } else {
      out <- lapply(out, function(x){
        x$Information$Approach_2ndorder <- NA
        x
      })
    }
    ## Set class
    class(out) <- c("cSEMResults", "cSEMResults_multi")

    ## Second order (2/3 stage approach)
    if(any(model_original$construct_order == "Second order") && 
       args$.approach_2ndorder %in% c("2stage", "mixed")) {
      
      ### Second step
      out <- lapply(out, function(x) {
        calculate2ndStage(model_original, x, args_needed, .approach_2ndorder)
        })
      
      ## Set class
      class(out) <- c("cSEMResults", "cSEMResults_multi", "cSEMResults_2ndorder")
    }
    
  } else if(any(model_original$construct_order == "Second order") && 
            args$.approach_2ndorder %in% c("2stage", "mixed")) {
    
    ### Second step
    out <- calculate2ndStage(model_original, out, args_needed, .approach_2ndorder)
    
  } else {
    ## Add second order approach to $Information
    if(any(model_original$construct_order == "Second order")) {
      out$Information$Approach_2ndorder  <- .approach_2ndorder
      out$Information$Model_original     <- model_original

    } else {
      out$Information$Approach_2ndorder <- NA
    }
    
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
      .force                = FALSE,
      .seed                 = .seed,
      .sign_change_option   = .sign_change_option
    )
  }
  
  return(out)
}
