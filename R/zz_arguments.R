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
#' @param .data A `data.frame` or a `matrix` containing the raw data. Possible
#'   data column types or classes are: logical, numeric (double or integer), factor 
#'   (ordered and unordered) or a mix of several types.
#' @param .model A model in \code{\link[lavaan:model.syntax]{lavaan model syntax}}
#'   or a [cSEMModel]-list.
#' @param .alpha An integer or a numeric vector of significance levels. 
#'   Defaults to `0.05`.
#' @param .approach Character string. The Kettenring approach to use. One of 
#' "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*" or "*GENVAR*". Defaults to
#' "*SUMCORR*".
#' @param .approach_2ndorder Character string. Approach used for models containing
#'   second order constructs. One of: "*2stage*" or "*repeated_indicators*". 
#'   Defaults to "*2stage*".
#' @param .approach_cor_robust Character string. Approach used to obtain a robust 
#'   indicator correlation matrix. One of: "*none*" in which case nothing is done,
#'   "*theil-sen*" or (TODO)
#'   Defaults to "*none*".
#' @param .approach_nl Character string. Approach used to estimate nonlinear
#'   structural relationships. One of: "*sequential*" or "*replace*".
#'   Defaults to "*sequential*".
#' @param .approach_paths Character string. Approach used to estimate the
#'   structural coefficients. One of: "*OLS*" or "*2SLS*".
#'   Defaults to "*OLS*".
#' @param .approach_weights Character string. Approach used to
#'   obtain composite weights. One of: "*PLS*", "*SUMCORR*", "*MAXVAR*",
#'   "*SSQCORR*", "*MINVAR*", "*GENVAR*", "*GSCA*", "*fixed*", or "*unit*".
#'   Defaults to "*PLS*".
#' @param .args_used A list of function argument names to `fun()` whose value 
#'   was modified by the user.
#' @param .C A (J x J) composite variance-covariance matrix.
#' @param .choices Logical. Should candidate values for the arguments be returned?
#'   Defaults to `FALSE`.
#' @param .conv_criterion Character string. The criterion to use for the convergence check.
#'   One of: "*diff_absolute*", "*diff_squared*", or "*diff_relative*". Defaults
#'   to "*diff_absolute*".
#' @param .csem_model A (possibly incomplete) [cSEMModel]-list.
#' @param .disattenuate Logical. If possible, should composite correlations be disattenuated
#'   if the construct is modeled as a common factor? Defaults to `TRUE`.
#' @param .distance Character string. A distance measure. One of: "*geodesic*"
#'   or "*squared_euclidian*". Defaults to "*geodesic*".
#' @param .dominant_indicators A character vector of `"name" = "value"` pairs, 
#'   where `"value"` is a character string giving the name of the dominant indicator
#'   and `"name"` a character string of the corresponding construct name.
#'   Dominant indicators may be specified for a subset of the constructs. 
#' @param .drop_inadmissibles Logical. Should inadmissible solutions be dropped? 
#'   Defaults to `TRUE`.
#' @param .E A (J x J) matrix of inner weights.
#' @param .estimate_structural Logical. Should the structural coefficients
#'   be estimated? Defaults to `TRUE`.
#' @param .H The (N x J) matrix of construct scores.
#' @param .id Character string. The name of the column of `.data` used to split
#'   the data into groups. Defaults to `NULL`.
#' @param .iter_max Integer. The maximum number of iterations allowed.
#'   If `iter_max = 1` and `.approach_weights = "PLS"` one-step weights are returned. 
#'   If the algorithm exceeds the specified number, weights of iteration step 
#'   `.iter_max - 1`  will be returned with a warning. Defaults to `100`.
#' @param .matrix1 A `matrix` to compare.
#' @param .matrix2 A `matrix` to compare.
#' @param .matrices A list of at least two matrices.
#' @param .modes A vector specifying the mode that should be used for
#'   each construct in the form `"name" = "mode"`, where `"name"` refers to the
#'   construct name and `"mode"`` is one of *"ModeA"* or *"ModeB"*.
#' @param .normality Logical. Should joint normality be assumed in the nonlinear model?
#'  For details see: \insertCite{Dijkstra2014;textual}{cSEM}. 
#'  Defaults to `TRUE`. Ignored if the model is linear.
#' @param .object An R object of class `cSEM<class>` with corresponding method.
#' @param .only_dots Logical. Should only arguments to be passed to lower level 
#'   functions via the  `...` argument of the `fun` function be returned. 
#'   Defaults to `FALSE`.
#' @param .P A (J x J) construct variance-covariance matrix (possibly disattenuated).
#' @param .parallel Logical. Use parallel computing. Defaults to `FALSE`. Note:
#'   requires the `doSNOW` and the `parallel` package to be installed.
#' @param .PLS_approach_cf Character string. Approach used to obtain the correction
#'   factors for PLSc. One of: "*dist_squared_euclid*", "*dist_euclid_weighted*",
#'   "*fisher_transformed*", "*mean_arithmetic*", "*mean_geometric*", "*mean_harmonic*",
#'   "*geo_of_harmonic*". Defaults to "*dist_euclid*". 
#'   Ignored if `.disattenuate = FALSE` or if `.approach_weights` is not PLS.
#' @param .PLS_ignore_structural_model Logical. Should the structural model be ignored
#'   when calculating the inner weights of the PLS algorithm? Defaults to `FALSE`.
#'   Ignored if `.approach_weights` is not PLS.
#' @param .PLS_modes Either a named vector specifying the mode that should be used for
#'   each construct in the form `"name" = "mode"`, a single character
#'   string giving the mode that should be used for all constructs, or `NULL`.
#'   Possible choices are: "*ModeA*" or "*ModeB*". Defaults to `NULL`.
#'   If `NULL`, `csem()` will choose the appropriate mode according to the type
#'   of construct used. Ignored if `.approach_weight` is not PLS.  
#' @param .PLS_weight_scheme_inner Character string. The inner weighting scheme
#'   used in PLS. One of: "*centroid*", "*factorial*", or "*path*".
#'   Defaults to "*path*". Ignored if `.approach_weight` is not PLS.
#' @param .Q A vector of composite-construct correlations with element names equal to
#'   the names of the J construct names used in the measurement model. Note 
#'   Q^2 is also called the reliability coefficient.
#' @param .reliabilities A character vector of `"name" = value` pairs, 
#'   where `value` is a number between 0 and 1 and `"name"` a character string
#'   of the corresponding construct name, or `NULL`. Reliabilities
#'   may be given for a subset of the constructs. Defaults to `NULL` in which case
#'   reliabilities are estimated by `csem()`.
#' @param .runs Integer. How many runs should be performed? Defaults to `499`.
#' @param .S The (K x K) empirical indicator correlation matrix.
#' @param .saturated Logical. Should a saturated structural model be used? Defaults to `FALSE`.
#' @param .show_progress Logical. Show progress bar. Defaults to `TRUE`.
#' @param .stage Character string. The stage the model is need for.
#'   One of "*first*" or "*second*". Defaults to "*first*".
#' @param .terms A vector of construct names to be classified.
#' @param .tolerance Double. The tolerance criterion for convergence. 
#'   Defaults to `1e-05`.
#' @param .verbose Logical. Should information be printed to the console? Defaults
#'   to `TRUE`.
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
#' args.default(.choices = TRUE)
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
    .data                    = NULL,
    .model                   = NULL,
    .alpha                   = 0.05,
    .approach                = c("SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR"),
    .approach_paths          = c("OLS", "2SLS"),
    .approach_weights        = c("PLS", "SUMCORR", "MAXVAR", "SSQCORR", "MINVAR", "GENVAR",
                                 "GSCA", "fixed", "unit"), 
    .C                       = NULL,
    .choices                 = FALSE,
    .csem_model              = NULL,
    .distance                = c("geodesic", "squared_euclidian"),
    .drop_inadmissibles      = TRUE,
    .E                       = NULL,
    .H                       = NULL,
    .id                      = NULL,
    .matrix1                 = NULL,
    .matrix2                 = NULL,
    .matrices                = NULL,
    .modes                   = NULL,
    .object                  = NULL,
    .only_dots               = FALSE,
    .P                       = NULL,
    .parallel                = FALSE,
    .runs                    = 499,
    .Q                       = NULL,
    .S                       = NULL,
    .saturated               = FALSE,
    .show_progress           = TRUE,
    .terms                   = NULL,
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
    .approach_2ndorder       = c("2stage", "repeated_indicators"),
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
    .approach_cor_robust     = c("none", "theil-sen"),
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
