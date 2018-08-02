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
#' @param .data A `data.frame` or a `matrix` containing the raw data. 
#' @param .model A model in \code{\link[lavaan:model.syntax]{lavaan model syntax}}
#'   or a [cSEMModel]-list.
#' @param .alpha An integer or a numeric vector of significance levels. 
#'   Defaults to `0.05`.
#' @param .approach Character string. The Kettenring approach to use. One of 
#' "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*" or "*GENVAR*". Defaults to
#' "*SUMCORR*".
#' @param .approach_nl Character string. Approach used to estimate nonlinear
#'   structural relationships. One of: "*sequential*" or "*replace*".
#'   Defaults to "*sequential*".
#' @param .approach_paths Character string. Approach used to estimate the
#'   structural coefficients. One of: "*OLS*", "*2SLS*", or "*3SLS*".
#'   Defaults to "*OLS*".
#' @param .approach_weights Character string. Approach used to
#'   obtain composite weights. One of: "*PLS*", "*SUMCORR*", "*MAXVAR*",
#'   "*SSQCORR*", "*MINVAR*", "*GENVAR*", "*GSCA*", "*fixed*", or "*unit*".
#'   Defaults to "*PLS*".
#' @param .args_used A list of function argument names to `foo()` whose value 
#'   was modified by the user. Usually captured by `as.list(match_call())[-1]` 
#'   within the body of `foo()`.
#' @param .C A (J x J) composite variance-covariance matrix.
#' @param .csem_model A (possibly incomplete) [cSEMModel]-list.
#' @param .disattenuate Logical. If possible, should composite correlations be disattenuated
#'   if the construct is modelled as a common factor? Defaults to `TRUE`.
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
#' @param .id Charachter string. The name of the column of `.data` used to split
#'   the date into groups. 
#' @param .iter_max Integer. The maximum number of iterations allowed.
#'   If `iter_max = 1` and `.approach_weights = "PLS"` one-step weights are returned. 
#'   If the algorithm exceeds the specified number, weights of iteration step 
#'   `.iter_max - 1`  will be returned with a warning. Defaults to `100`.
#' @param .matrix1 A `matrix` to compare.
#' @param .matrix2 A `matrix` to compare.
#' @param .modes A named vector giving the mode used to obtain new outer weights.
#' @param .normality Logical. Should joint normality be assumed in the nonlinear model?
#'  For details see: \insertCite{Dijkstra2014;textual}{cSEM}. 
#'  Defaults to `TRUE`. Ignored if the model is linear.
#' @param .object An R object of class `cSEM<class>` with corresponding method.
#' @param .only_dots Logical. Should only arguments to be passed to lower level 
#'   functions via the  `...` argument of the [csem()] or [cca()] function be returned. 
#'   Defaults to `FALSE`.
#' @param .PLS_approach_cf Character string. Approach used to obtain the correction
#'   factors for PLSc. One of: "*dist_euclid*", "*dist_euclid_weighted*",
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
#'   Defaults to "*centroid*". Ignored if `.approach_weight` is not PLS.
#' @param .Q A vector of composite-construct correlations with element names equal to
#'   the names of the J construct names used in the measurement model. Note 
#'   Q^2 is also called the reliability coefficient.
#' @param .reliabilities A character vector of `"name" = value` pairs, 
#'   where `value` is a number between 0 and 1 and `"name"` a character string
#'   of the corresponding construct name, or `NULL`. Reliabilities
#'   may be given for a subset of the constructs. Defaults to `NULL` in which case
#'   reliabilites are estimated by `csem()`.
#' @param .runs Integer. How many runs should be performed? Defaults to `499`.
#' @param .S The (K x K) empirical indicator correlation matrix.
#' @param .terms A vector of construct names to be classified.
#' @param .tolerance Double. The tolerance criterion for convergence. 
#'   Defaults to `1e-05`.
#' @param .verbose Logical. Should information be printed to the console? Defaults
#'   to `TRUE`.
#' @param .W A list with elements `$W`, `$E`, `Modes`, `Conv_status`, and 
#'   `Iterations`.
#' @param ... Further arguments to be passed down to other functions.
#' See [args_default()] with `.only_dots = TRUE` for a complete list 
#' of possible arguments and their defaults.
#'
#' @name csem_arguments
#' @aliases cSEMArguments csem_parameters
NULL

#' List all arguments and their defaults
#'
#' The function `args_default()` lists all possible arguments and their defaults.
#' For a description of the arguments see [csem_arguments].
#'
#' This is the complete list of arguments, their current defaults, and their
#' possible alternatives. Each default may be changed by passing a `arg_name = value`
#' pair to the function that uses the argument. Changes to the default arguments
#' are checked and processed by the [handleArgs()] function.
#' 
#' If `.only_dots = TRUE` all possible arguments to be passed to lower level functions via the 
#' `...` argument of the [csem()] or [cca()] function are returned.
#'
#' @usage args_default(.only_dots = FALSE)
#' 
#' @param csem_arguments
#' 
#' @return A named list of argument names and their defaults.
#'
#' @examples
#'
#' args_default() # List all possible arguments and their defaults.
#' args_default(.only_dots = TRUE) # list only those accepted by `...`.
#'
#' @seealso [handleArgs()], [csem_arguments], [cca()], [csem()], [foreman()]
#'
#' @export

args_default <- function(.only_dots = FALSE) {
  
  args <- list(
    .data                    = NULL,
    .model                   = NULL,
    .alpha                   = 0.05,
    .approach                = "SUMCORR",
    .approach_nl             = "sequential",
    .approach_paths          = "OLS",
    .approach_weights        = "PLS", 
    .C                       = NULL,
    .csem_model              = NULL,
    .disattenuate            = TRUE,
    .drop_inadmissibles      = TRUE,
    .E                       = NULL,
    .estimate_structural     = TRUE,
    .H                       = NULL,
    .id                      = NULL,
    .matrix1                 = NULL,
    .matrix2                 = NULL,
    .mode                    = NULL,
    .object                  = NULL,
    .only_dots               = FALSE,
    .PLS_modes               = NULL,
    .runs                    = 499,
    .Q                       = NULL,
    .reliabilities           = NULL,
    .S                       = NULL,
    .terms                   = NULL,
    .verbose                 = TRUE,
    .W                       = NULL
  )
  
  args_dotdotdot <- list(
    # Arguments passed to calculateWeightsPLS
    .iter_max                = 100,
    .tolerance               = 1e-05,
    
    # Arguments passed to calculateInnerWeightsPLS
    .PLS_ignore_structural_model = FALSE,
    .PLS_weight_scheme_inner     = "centroid",
    
    # Arguments passed to calculateCorrectionFactors
    .PLS_approach_cf         = "dist_euclid",
    
    # Arguements passed to estimatorPathOLS
    .normality               = TRUE,
    
    #  Arguments passed to foreman
    .dominant_indicators     = NULL
    
  )

  if(.only_dots) {
    return(args_dotdotdot)
  } else {
    
    args_all <- c(args, args_dotdotdot)
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
#'
handleArgs <- function(.args_used) {
  
  x <- args_default()
  
  args_default_names  <- names(x)
  args_used_names <- names(.args_used)
  
  ## Are all argument names used valid?
  args_diff <- setdiff(args_used_names, args_default_names)
  
  if (length(args_diff) > 0) {
    stop("The following argument(s) are not known to cSEM: ",
         paste0("`", args_diff, '`', collapse = ", "), "\n",
         "Please type 'args_default()' to find all valid arguments and their defaults. Or check ?args_default.",
         call. = FALSE)
  }
  
  ## Which arguments need to be changed?
  args_intersect <- intersect(args_used_names, args_default_names)
  
  ## Replace all arguments that were changed or explicitly given and keep
  #  the default values for the others
  
  for(i in args_intersect) {
    x[i] <- .args_used[i]
  }
  
  return(x)
}
