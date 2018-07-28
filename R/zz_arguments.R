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
#' @param .model A model in \href{http://lavaan.ugent.be/tutorial/syntax1.html}{lavaan model syntax}
#'   or a [cSEMModel] object.
#' @param .alpha An integer or a numeric vector of significance levels. 
#'   Defaults to `0.05`.
#' @param .approach Character string. The Kettenring approach to use. One of 
#' "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*" or "*GENVAR*". Defaults to
#' "*SUMCORR*".
#' @param .approach_nl Character string. The approach used to handle non-linear
#'   structural relationships. Possible choices are: "*sequential*", or "*replace*".
#'   Defaults to "*sequential*".
#' @param .approach_paths Character string. Approach used to estimate the
#'   structural coefficients. Possible choices are: "*OLS*", "*2SLS*", or "*3SLS*".
#'   Defaults to "*OLS*".
#' @param .approach_weights Character string. The name of the approach used to
#'   obtain composite weights. Possible choices are: "*PLS*", "*SUMCORR*", "*MAXVAR*",
#'   "*SSQCORR*", "*MINVAR*", "*GENVAR*", "*GSCA*", "*fixed*", or "*unit*".
#'   Defaults to "*PLS*".
#' @param .args_used A list of argument names and user picked values. Usually captured by 
#'   `as.list(match_call())[-1]`.
#' @param .C A (J x J) proxy variance-covariance matrix.
#' @param .csem_model A (possibly incomplete) [cSEMModel] list.
#' @param .disattenuate Logical. Should proxy correlations be disattenuated
#'   if the construct is modeled as a common factor? Defaults to `TRUE`.
#' @param .dominant_indicators A character vector of `name = value` pairs, where `value` is 
#'   a character string giving the name of the dominant indicator and `name` 
#'   the corresponding construct name. Dominant indicators may be specified for 
#'   a subset of the constructs. 
#' @param .dropInadmissibles Logical. Should the inadmissible solution be dropped? Defaults to `TRUE`.
#' @param .E A (J x J) matrix of inner weights.
#' @param .estimate_structural Logical. Should the structural (path) coefficients
#'   be estimated? Defaults to `TRUE`.
#' @param .group_var  Character string. The name of the column used to split the data into groups.
#' @param .H The (N x J) matrix of proxy values.
#' @param .id Charachter string. The name of the column of `.data` used to split
#'   the date into groups. 
#' @param .ignore_structural_model Logical. Should the structural (path) model be ignored
#'   when calculating the inner weights of the PLS algorithm? Defaults to `FALSE`.
#' @param .iter_max Integer. The maximum number of iterations of the PLS algorithm.
#'   If `iter_max = 1` one-step weights are returned. If the algorithm exceeds
#'   the specified number, weights of iteration step `.iter_max - 1`  will be returned
#'   with a warning. Defaults to `100`. The argument is ignored if
#'   `.approach_weight` is not PLS.
#' @param .matrix1 A matrix. Defaults to `NULL`.
#' @param .matrix2 A matrix. Defaults to `NULL`. 
#' @param .modes A named vector giving the mode used to obtain new outer weights.
#' @param .normality Logical. Should joint normality be assumed if the model
#'   contains non-linear terms. For details see: \insertCite{Dijkstra2014;textual}{cSEM}. 
#'   Defaults to `TRUE`.
#' @param .object An R object of class [cSEMResults].
#' @param .only_dots Logical. Should only those arguments to be passed to lower level functions via the 
#' `...` argument of the [csem] or [cca] function be returned. Defaults to `FALSE`.
#' @param .permutations Integer. The number permutations used for step 2 and 3 of the test.
#' Defaults to `100`. Note however that the number should typically be a lot higher.
#' @param .approach_cf Character string. The approach to obtain the correction
#'   factor for PLSc. Possible choices are: "*dist_euclid*", "*dist_euclid_weighted*",
#'   "*fisher_transformed*", "*mean_arithmetic*", "*mean_geometric*", "*mean_harmonic*",
#'   "*geo_of_harmonic*". Defaults to "*dist_euclid*". Ignored if `.disattenuate = FALSE`.
#' @param .PLS_weight_scheme_inner Character string. The inner weighting scheme
#'   used in PLS. Possible choices are: "*centroid*", "*factorial*", or "*path*".
#'   Defaults to "*centroid*". The argument is ignored if `.approach_weight` is not PLS.
#' @param .PLS_mode Either a named vector specifying the mode that should be used for
#'   each construct in the form `name = "mode"`, a single character
#'   string giving the mode that should be used for all constructs or `NULL`.
#'   Possible choices are: "*ModeA*" or "*ModeB*". Defaults to `NULL`.
#'   If `NULL`, `cSEM` will choose the appropriate mode according to the type
#'   of construct used. The argument is ignored if `.approach_weight` is not PLS.
#' @param .Q A vector of proxy-construct correlations with element names equal to
#'   the names of the J construct names used in the measurement model. Note 
#'   Q^2 is also called the reliability coefficient.
#' @param .reliabilities A vector of `name = value` pairs of reliability coefficients. 
#'   Element names are equal to the names of the J construct names. Reliabilities
#'   may be given for a subset of the constructs. Defaults to `NULL` in which case
#'   reliabilites are estimated by `cSEM`.
#' @param .runs Integer. How many runs should be performed? Defaults to `499`.
#' @param .S The (scaled) empirical (K x K) indicator covariance/correlation matrix,
#'   where K is the number of indicators.
#' @param .terms A vector of construct names to be classified.
#' @param .tolerance Double. The tolerance criterion for convergence of the PLS
#'   algorithm. Defaults to `1e-05`.
#'   The argument is ignored if `.approach_weight` is not PLS.
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
#' are checked and processed by the [handleArgs] function.
#' 
#' If `.only_dots = TRUE` all possible arguments to be passed to lower level functions via the 
#' `...` argument of the [csem] or [cca] function are returned.
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
#' @seealso [handleArgs], [csem_arguments], [cca], [csem], [foreman]
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
    .dropInadmissibles       = TRUE,
    .E                       = NULL,
    .estimate_structural     = TRUE,
    .group_var               = NULL,
    .H                       = NULL,
    .id                      = NULL,
    .matrix1                 = NULL,
    .matrix2                 = NULL,
    .mode                    = NULL,
    .object                  = NULL,
    .only_dots               = FALSE,
    .permutations            = 100,
    .PLS_mode                = NULL,
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
    .ignore_structural_model = FALSE,
    .PLS_weight_scheme_inner = "centroid",
    
    # Arguments passed to calculateCorrectionFactors
    .approach_cf             = "dist_euclid",
    
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
