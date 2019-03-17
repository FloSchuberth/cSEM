#' cSEMArguments
#'
#' An alphabetical list of all arguments used by functions of the `cSEM` package,
#' including their description and defaults.
#' Mainly used for internal purposes (parameter inheritance). To list all arguments
#' and their defaults, use [args_default()]. To list only those arguments to be 
#' passed down to lower level functions via the `...` argument of any 
#' function having `...` as a formal argument, set the `.only_dots` argument 
#' of the [args_default()] function to `TRUE`.
#'
#' @param .alpha An integer or a numeric vector of significance levels. 
#'   Defaults to `0.05`.
#' @param .approach Character string. The Kettenring approach to use. One of 
#' "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*" or "*GENVAR*". Defaults to
#' "*SUMCORR*".
#' @param .approach_2ndorder Character string. Approach used for models containing
#'   second order constructs. One of: "*3stage*" or "*repeated_indicators*". 
#'   Defaults to "*3stage*".
#' @param .approach_cor_robust Character string. Approach used to obtain a robust 
#'   indicator correlation matrix. One of: "*none*" in which case nothing is done,
#'   "*spearman*" for the spearman correlation, or
#'   "*mcd*" via \code{\link[MASS:cov.rob]{MASS::cov.rob()}} for a robust correlation matrix. 
#'   Defaults to "*none*".
#' @param .approach_nl Character string. Approach used to estimate nonlinear
#'   structural relationships. One of: "*sequential*" or "*replace*".
#'   Defaults to "*sequential*".
#' @param .approach_paths Character string. Approach used to estimate the
#'   structural coefficients. One of: "*OLS*" or "*2SLS*" (not yet implemented).
#'   Defaults to "*OLS*".
#' @param .approach_weights Character string. Approach used to
#'   obtain composite weights. One of: "*PLS-PM*", "*SUMCORR*", "*MAXVAR*",
#'   "*SSQCORR*", "*MINVAR*", "*GENVAR*", "*GSCA*", "*PCA*", "*unit*", "*bartlett*", 
#'   or "*regression*".
#'   Defaults to "*PLS-PM*".
#' @param .args_used A list of function argument names to `fun()` whose value 
#'   was modified by the user.
#' @param .bias_corrected Logical. Should the standard and the tStat
#'   confidence intervall be bias-corrected using the bootstraped bias estimate? 
#'   If `TRUE` the confidence intervall for some estimated parameter `theta` 
#'   is centered at `2*theta - theta*_hat`,
#'   where `theta*_hat` is the average over all .R bootstrap estimates of `theta`.
#'   Defaults to `TRUE`
#' @param .C A (J x J) composite variance-covariance matrix.
#' @param .choices Logical. Should candidate values for the arguments be returned?
#'   Defaults to `FALSE`.
#' @param .ci A vector of character strings naming the confidence interval to compute.
#'   For possible choices see [infer()].
#' @param .conv_criterion Character string. The criterion to use for the convergence check.
#'   One of: "*diff_absolute*", "*diff_squared*", or "*diff_relative*". Defaults
#'   to "*diff_absolute*".
#' @param .csem_model A (possibly incomplete) [cSEMModel]-list.
#' @param .csem_resample A list resulting from a call to [resamplecSEMResults()].
#' @param .cv_folds Integer. The number of cross-validation folds to use. Setting
#'   `.cv_folds` to `N` (the number of observations) produces 
#'   leave-one-out cross-validation samples. Defaults to `10`.
#' @param .data A `data.frame` or a `matrix` of standardized or unstandarized data. 
#'   Possible column types or classes of the data provided are: logical, 
#'   numeric (double or integer), factor (ordered and unordered) 
#'   or a mix of several types. The data may also include
#'   *one* character column whose column name must be given to `.id`. 
#'   This column is assumed to contain group identifiers used to split 
#'   the data into groups.
#' @param .disattenuate Logical. If possible, should composite/proxy correlations 
#'   be disattenuated if the construct is modeled as a common factor? 
#'   Defaults to `TRUE`.
#' @param .dist Character string. The distribution to use for the critical value.
#'  One of *"t"* for Student's t-distribution or *"z"* for the standard normal distribution.
#'  Defaults to *"z"*.
#' @param .distance Character string. A distance measure. One of: "*geodesic*"
#'   or "*squared_euclidian*". Defaults to "*geodesic*".
#' @param .df Character string. The method for obtaining the degrees of freedom.
#'   Choices are "*type1*" and "*type2*". Defaults to "*type1*" .
#' @param .dominant_indicators A character vector of `"name" = "value"` pairs, 
#'   where `"value"` is a character string giving the name of the dominant indicator
#'   and `"name"` a character string of the corresponding construct name.
#'   Dominant indicators may be specified for a subset of the constructs. 
#' @param .E A (J x J) matrix of inner weights.
#' @param .estimate_structural Logical. Should the structural coefficients
#'   be estimated? Defaults to `TRUE`.
#' @param .eval_plan Character string. The evaluation plan to use. One of 
#'   "*sequential*" or "*multiprocess*". In the latter case 
#'   all available cores will be used. Defaults to "*sequential*".
#' @param .first_resample A list containing the `.R` resamples based on the original
#'   data obtained by resamplecSEMResults()
#' @param .second_resample A list containing `.R2` resamples for each of the `.R`
#'   resamples of the first run.
#' @param .H The (N x J) matrix of construct scores.
#' @param .handle_inadmissibles Character string. How should inadmissible results 
#'   be treated? One of "*drop*", "*ignore*", or "*replace*". If "*drop*", all
#'   replications/resamples yielding an inadmissible result will be dropped (
#'   the number of results returned will be less than .R). For "*ignore*" all results are returned 
#'   even if they are inadmissible (number of results returned = .R). For "*replace*"
#'   resampling continues until there are exactly .R admissible solutions. 
#'   Defaults to "*drop*".
#' @param .id Character string or integer. The name or position of the column of 
#'   `.data` used to split the data into groups.
#'    Defaults to `NULL`.
#' @param .iter_max Integer. The maximum number of iterations allowed.
#'   If `iter_max = 1` and `.approach_weights = "PLS-PM"` one-step weights are returned. 
#'   If the algorithm exceeds the specified number, weights of iteration step 
#'   `.iter_max - 1`  will be returned with a warning. Defaults to `100`.
#' @param .n Integer. The number of observations of the original data.
#' @param .matrix1 A `matrix` to compare.
#' @param .matrix2 A `matrix` to compare.
#' @param .matrices A list of at least two matrices.
#' @param .model A model in \code{\link[lavaan:model.syntax]{lavaan model syntax}}
#'   or a [cSEMModel]-list.
#' @param .modes A vector giving the mode for each construct in the form `"name" = "mode"`. 
#'   Only used internally. 
#' @param .normality Logical. Should joint normality be assumed in the nonlinear model? 
#'  Defaults to `TRUE`. Ignored if the model is linear.
#' @param .object An R object of class [cSEMResults] resulting from a call to [csem()].
#' @param .only_common_factors Logical. Should only common factors be included? 
#'   Defaults to `FALSE`.
#' @param .only_dots Logical. Should only arguments to be passed to lower level 
#'   functions via the  `...` argument of the `fun` function be returned. 
#'   Defaults to `FALSE`.
#' @param .P A (J x J) construct variance-covariance matrix (possibly disattenuated).
#' @param .parallel Logical. Use parallel computing. Defaults to `FALSE`. Note:
#'   requires the `doSNOW` and the `parallel` package to be installed.
#' @param .PLS_approach_cf Character string. Approach used to obtain the correction
#'   factors for PLSc. One of: "*dist_squared_euclid*", "*dist_euclid_weighted*",
#'   "*fisher_transformed*", "*mean_arithmetic*", "*mean_geometric*", "*mean_harmonic*",
#'   "*geo_of_harmonic*". Defaults to "*dist_squared_euclid*". 
#'   Ignored if `.disattenuate = FALSE` or if `.approach_weights` is not PLS-PM.
#' @param .PLS_ignore_structural_model Logical. Should the structural model be ignored
#'   when calculating the inner weights of the PLS-PM algorithm? Defaults to `FALSE`.
#'   Ignored if `.approach_weights` is not PLS-PM.
#' @param .PLS_modes Either a named list specifying the mode that should be used for
#'   each construct in the form `"name" = "mode"`, a single character
#'   string giving the mode that should be used for all constructs, or `NULL`.
#'   Possible choices for `"mode"` are: "*modeA*", "*modeB*", "*modeBNNLS*", "*unit*", a single number (weight) or 
#'   a vector of fixed weights of the same length as there are indicators for the
#'   construct given by `"name"`. If only a single number is provided this is identical to
#'   using unit weights, as weights are rescaled such that the related composite 
#'   has unit variance.  Defaults to `NULL`.
#'   If `NULL` the appropriate mode according to the type
#'   of construct used is choosen. Ignored if `.approach_weight` is not PLS-PM.  
#' @param .PLS_weight_scheme_inner Character string. The inner weighting scheme
#'   used in PLS-PM. One of: "*centroid*", "*factorial*", or "*path*".
#'   Defaults to "*path*". Ignored if `.approach_weight` is not PLS-PM.
#' @param .probs A vector of probabilities.
#' @param .quantity Character string. Which statistic should be returned?
#'   One of (TODO) 
#'   Defaults to (TODO).
#' @param .Q A vector of composite-construct correlations with element names equal to
#'   the names of the J construct names used in the measurement model. Note 
#'   Q^2 is also called the reliability coefficient.
#' @param .reliabilities A character vector of `"name" = value` pairs, 
#'   where `value` is a number between 0 and 1 and `"name"` a character string
#'   of the corresponding construct name, or `NULL`. Reliabilities
#'   may be given for a subset of the constructs. Defaults to `NULL` in which case
#'   reliabilities are estimated by `csem()`. Currently, only supported for
#'   `.approach_weights = "PLS-PM"`.
#' @param .resample_method Character string. The resampling method to use. One of: 
#'  "*bootstrap*" or "*jackknife*". Defaults to "*bootstrap*".
#' @param .resample_method2 Character string. The resampling method to use when resampling
#'   from a resample. One of: "*none*", "*bootstrap*" or "*jackknife*". For 
#'   "*bootstrap*" the number of draws is provided via `.R2`. Currently, 
#'   resampling from each resample is only required for the studentized confidence
#'   intervall computed by the [infer()] function. Defaults to "*none*".  
#' @param `.resample_object` An R object of class `cSEMResults_resampled`
#'   obtained from [resamplecSEMResults()] or by setting `.resample_method = "bootstrap"`
#'   or `"jackknife"` when calling [csem()].
#' @param .R Integer. The number of bootstrap replications. Defaults to `499`.
#' @param .R2 Integer. The number of bootstrap replications to use when 
#'   resampling from a resample. Defaults to `199`.
#' @param .S The (K x K) empirical indicator correlation matrix.
#' @param .saturated Logical. Should a saturated structural model be used? Defaults to `FALSE`.
#' @param .seed Integer. The random seed to use. Defaults to `NULL` in which
#'   case an arbitrary seed is choosen.
#' @param .stage Character string. The stage the model is need for.
#'   One of "*first*" or "*second*". Defaults to "*first*".
#' @param .terms A vector of construct names to be classified.
#' @param .tolerance Double. The tolerance criterion for convergence. 
#'   Defaults to `1e-05`.
#' @param .type_vcv Character string. Indicates which model-implied correlation matrix is calcuted
#'  One of "*indicator*" or "*construct*". Defaults to "*indicator*".   
#' @param .verbose Logical. Should information be printed to the console? Defaults
#'   to `TRUE`.
#' @param .user_funs A function or a (named) list of functions to apply to every
#'   resample. Takes `.object` as an input (e.g., `myFun <- function(.object) {...}`).
#'   Output should preferably be a (named)
#'   vector but matrices are also accepted. However, the output will be 
#'   vectorized (columnwise) in this case. See the examples section for details.
#' @param .vector1 A vector of numeric values.
#' @param .vector2 A vector of numeric values.
#' @param .W A (J x K) matrix of weights.
#' @param .W_new A (J x K) matrix of weights.
#' @param .W_old A (J x K) matrix of weights.
#' @param .which_fun Character string. The `...` argument names and values of which 
#'   function should be returned? One of: `"csem"` or `"cca"`. Defaults to `"csem"`. 
#'   Currently ignored if `.only_dots = FALSE`. 
#' @param .X A matrix of processed data (scaled, cleaned and ordered).
#' @param .X_cleaned A data.frame of processed data (cleaned and ordered). Note: `X_cleaned`
#'   may not be scaled!
#'
#' @name csem_arguments
#' @aliases cSEMArguments
#' 
NULL

#' Show argument defaults or candidates
#'
#' Show all arguments used by package functions including default or candidate
#' values. For argument descriptions, see: [csem_arguments].
#'
#' By default `args_default()`returns a list of default values by argument name.
#' If the list of accepted candidate values is required instead, use `.choices = TRUE`.
#' 
#' If `.only_dots = TRUE` arguments passed to lower level 
#' functions via the `...` argument of `which.fun` are returned.
#'
#' @usage args_default(
#'   .choices   = FALSE,
#'   .only_dots = FALSE,
#'   .which_fun = "csem"
#' )
#' 
#' @inheritParams  csem_arguments
#' 
#' @return A named list of argument names and defaults or accepted candidates.
#'
#' @examples
#'
#' args_default() # List all possible arguments and their defaults.
#' args_default(.only_dots = TRUE) # list only those accepted by the `...` argument
#'                                 # of the csem function.
#' args_default(.only_dots = TRUE, .which_fun = "csem") # the same
#' 
#' ## Show accepted candidates:
#' args_default(.choices = TRUE)
#'
#' @seealso [handleArgs()], [csem_arguments], [cca()], [csem()], [foreman()]
#'
#' @export

args_default <- function(
  .choices   = FALSE, 
  .only_dots = FALSE, 
  .which_fun = "csem"
  ) {
  
  args <- list(
    .alpha                   = 0.05,
    .approach                = c("SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR"),
    .approach_paths          = c("OLS", "2SLS"),
    .approach_weights        = c("PLS-PM", "SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR",
                                 "GSCA", "PCA", "unit", "bartlett", "regression"), 
    .arguments               = NULL,
    .bias_corrected          = TRUE,
    .C                       = NULL,
    .choices                 = FALSE,
    .ci                      = c("CI_standard_z", "CI_standard_t", "CI_percentile", 
                                 "CI_basic", "CI_bc", "CI_bca", "CI_t_intervall"),
    .csem_model              = NULL,
    .csem_resample           = NULL,
    .cv_folds                = 10,
    .data                    = NULL,
    .dist                    = c("z", "t"),
    .distance                = c("geodesic", "squared_euclidian"),
    .df                      = c("type1", "type2"),
    .E                       = NULL,
    .eval_plan               = c("sequential", "multiprocess"),
    .first_resample          = NULL,
    .handle_inadmissibles    = c("drop", "ignore", "replace"),
    .H                       = NULL,
    .id                      = NULL,
    .listMatrices            = NULL, 
    .matrix1                 = NULL,
    .matrix2                 = NULL,
    .matrices                = NULL,
    .model                   = NULL,
    .modes                   = NULL,
    .only_common_factors     = TRUE,
    .object                  = NULL,
    .only_dots               = FALSE,
    .P                       = NULL,
    .parallel                = FALSE,
    .probs                   = NULL,
    .quantity                = c("all", "mean", "sd", "bias", "CI_standard_z", "CI_standard_t",
                                 "CI_percentile", "CI_basic", "CI_bc", "CI_bca", "CI_t_intervall"),
    .Q                       = NULL,
    .R                       = 499,
    .R2                      = 199,
    .resample_method         = c("none", "bootstrap", "jackknife"),
    .resample_method2        = c("none", "bootstrap", "jackknife"),
    .resample_object         = NULL,
    .S                       = NULL,
    .saturated               = FALSE,
    .second_resample         = NULL,
    .seed                    = NULL,
    .terms                   = NULL,
    .type_vcv                = c("indicator", "construct"),
    .user_funs               = NULL,
    .verbose                 = TRUE,
    .W                       = NULL,
    .which_fun               = c("csem", "cca"),
    .x                       = NULL,
    .X                       = NULL,
    .X_cleaned               = NULL,
    .y                       = NULL
  )
  
  args_dotdotdot_csem <- list(
    # Arguments passed to convertModel()
    .approach_2ndorder       = c("3stage", "repeated_indicators"),
    .stage                   = c("first", "second"),
    # Arguments passed to calculateWeightsPLS
    .iter_max                = 100,
    .PLS_modes               = NULL,
    .tolerance               = 1e-05,
    
    # Arguments passed to calculateInnerWeightsPLS
    .PLS_ignore_structural_model = FALSE,
    .PLS_weight_scheme_inner     = c("path", "centroid", "factorial"),
    
    # Arguments passed to calculateCorrectionFactors
    .PLS_approach_cf         = c("dist_squared_euclid", "dist_euclid_weighted", 
                                 "fisher_transformed", "mean_arithmetic",
                                 "mean_geometric", "mean_harmonic",
                                 "geo_of_harmonic"),
    # Arguments passed to checkConvergence
    .conv_criterion          = c("diff_absolute", "diff_squared", "diff_relative"),
    
    # Arguements passed to estimatorPathOLS
    .approach_nl             = c("sequential", "replace"),
    .normality               = TRUE,
    
    #  Arguments passed to foreman
    .approach_cor_robust     = c("none", "mcd", "spearman"),
    .disattenuate            = TRUE,
    .dominant_indicators     = NULL,
    .estimate_structural     = TRUE,
    .reliabilities           = NULL
  )
  
  args_dotdotdot_cca <- args_dotdotdot_csem
  args_dotdotdot_cca[c(".approach_nl", ".normality", ".estimate_structural")] <- NULL

  if(!.choices) {
    args <- lapply(args, function(x) eval(x)[1])
    args_dotdotdot_csem <- lapply(args_dotdotdot_csem, function(x) eval(x)[1])
    args_dotdotdot_cca <- lapply(args_dotdotdot_cca, function(x) eval(x)[1])
  }
  
  if(.only_dots) {
    cat("Possible `...` argument names and ", 
        ifelse(.choices, "possible candidates","default values"), " for function: `", .which_fun, "()` \n\n", 
        sep = "", "See `?cSEMArguments` for their description.\n\n")
    switch (.which_fun,
      "csem" = return(args_dotdotdot_csem),
      "cca"  = return(args_dotdotdot_cca),
      stop("Function `", .which_fun, "()` does not have `...` arguments",
           call. = FALSE)
    )
  } else {
    
    args_all <- c(args, args_dotdotdot_csem)
    args_all <- args_all[sort(names(args_all))]
    
    return(args_all)
  }

}

#' Internal: Handle arguments
#'
#' Internal helper function to handle arguments passed to any function within `cSEM`.
#'
#' @param .args_used A list of argument names and user picked values.
#'
#' @return The [args_default] list, with default values
#' changed to the values given by the user.
#' @keywords internal

handleArgs <- function(.args_used) {
  
  args_default  <- args_default()
  
  args_default_names  <- names(args_default)
  args_used_names <- names(.args_used)
  
  ## Are all argument names used valid?
  args_diff <- setdiff(args_used_names, args_default_names)
  
  if (length(args_diff) > 0) {
    stop("The following argument(s) are not known to cSEM: ",
         paste0("`", args_diff, '`', collapse = ", "), "\n",
         "Please type 'args_default()' to find all valid arguments and their defaults. Or check ?args_default.",
         call. = FALSE)
  }
  
  ## Are arguments valid
  # Note: for now, I will only check if string arguments are valid. In the
  # future we could refine and check if e.g. numeric values are within a reasonable
  # range or larger than zero if required etc.
  # choices_logical <- Filter(function(x) any(is.logical(x)), args_default(.choices = TRUE))
  # choices_numeric <- Filter(function(x) any(is.numeric(x)), args_default(.choices = TRUE))
  choices_character <- Filter(function(x) any(is.character(x)), args_default(.choices = TRUE))

  character_args <- intersect(names(choices_character), args_used_names)
  x <- Map(function(x, y) x %in% y, 
      x = .args_used[character_args], 
      y = choices_character[character_args]
      )

  lapply(seq_along(x), function(i) {
    if(isFALSE(x[[i]])) {
      n <- names(x[i])
      a <- args_default(.choices = TRUE)[[n]]
      stop(paste0("`", .args_used[n], "` is not a valid choice for ", "`", 
                  names(.args_used[n]),"`.\n"), 
           "Choices are: ", paste0("`", a[-length(a)],"`", collapse = ", "), 
           " or " , paste0("`", a[length(a)], "`"), call. = FALSE)
    }

  })
  ## Replace all arguments that were changed or explicitly given and keep
  #  the default values for the others
  
  for(i in args_used_names) {
    args_default[i] <- .args_used[i]
  }
  
  return(args_default)
}
