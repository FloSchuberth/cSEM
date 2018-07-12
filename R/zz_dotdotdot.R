#' Show possible "..." arguments and their defaults
#'
#' The function `args_default` lists all possible arguments to be passed to lower
#' level functions via the `...` argument of the [csem] or [cca] function.
#'
#' This is the complete list of arguments, their current defaults, and their
#' possible alternatives. Each default may be changed by passing a `arg_name = value`
#' pair to [csem] or [cca]. Arguments are processed by the internal [handleArgs]
#' function and, if valid, passed down to their corresponding function
#' via the [workhorse] function.
#'
#' **Arguments accepted by [csem]**
#' \describe{
#'   \item{`.approach_cf`}{*Default: __"dist_euclid"__*. A character string
#'   of the name of the approach used to obtain the correction factor for PLSc.
#'
#'   Currently seven methods are available (see: \insertCite{Dijkstra2013}{cSEM}
#'   for details on the different approaches):
#'
#'   \itemize{
#'     \item "dist_euclid"`
#'     \item "dist_euclid_weighted"
#'     \item "fisher_transformed"
#'     \item "mean_geometric"
#'     \item "mean_harmonic"
#'     \item "mean_arithmetic"
#'     \item "geo_of_harmonic" (not yet implemented)
#'   }
#'
#'   The argument is passed to: [calculateCorrectionFactors]}
#'   \item{`.ignore_structural_model`}{*Default: FALSE*. Should the structural (path)
#'   model be ignored when calculating the inner weights of the PLS algorithm?
#'
#'   The argument is passed to: [calculateInnerWeightsPLS]}
#'   \item{`.iter_max`}{*Default: __100__*. An integer indicating the maximum number of
#'   iterations of the PLS algorithm. If `iter_max = 1` one-step weights
#'   are returned. If the algorithm exceeds the specified number, weights of the
#'   previous iteration step will be returned with a warning.
#'
#'   The argument is passed to: [calculateWeightsPLS]}
#'   \item{`.tolerance`}{*Default: __1e-06__*. A numeric value specifying the tolerance
#'   criterion for convergence of the PLS algorithm.
#'
#'   The argument is passed to: [calculateWeightsPLS]}
#' }
#'
#' @usage args_default()
#'
#' @return A named list of argument names and their defaults.
#'
#' @examples
#' ## List all possible arguments to be passed to lower level functions via
#'    the ... argument of the csem or cca function.
#'
#' args_default()
#'
#' @seealso [handleArgs], [cca], [csem], [workhorse]
#'
#' @references
#' \insertRef{Dijkstra2013}{cSEM}
#'
#' @export

args_default <- function() {

  args_ls <- list(
    # Arguments passed to calculateWeightsPLS
    .iter_max                 = formals(workhorse)$.iter_max,
    .tolerance                = formals(workhorse)$.tolerance,


    # Arguments passed to calculateInnerWeightsPLS
    .ignore_structural_model  = formals(workhorse)$.ignore_structural_model,
    .PLS_weight_scheme_inner  = eval(formals(workhorse)$.PLS_weight_scheme_inner)[1],

    # Arguments passed to calculateCorrectionFactors
    .approach_cf              = eval(formals(workhorse)$.approach_cf)[1],

    # Arguements passed to estimator
    .normal                   = formals(workhorse)$.normal,
    
    #  Argument only used in workhorse
    .dominantIndicators = formals(workhorse)$.dominantIndicators
    
    
  )

  return(args_ls)
}

#' Internal: Handle "..." arguments
#'
#' Internal helper function to handle arguments passed via the `...` argument
#' to [csem] or [cca].
#'
#' @param .dotdotdot A named list of argument names and defaults usually captured by `list(...)`.
#'
#' @return The function returns the [args_default] list, with defaults
#' changed if necessary.
#'

handleArgs <- function(.dotdotdot) {

  x <- args_default()

  args_default_names <- names(x)
  args_dots_names    <- names(.dotdotdot)

  ## Are all argument names valid?
  args_diff <- setdiff(args_dots_names, args_default_names)
  if (length(args_diff) > 0) {
    stop("The following argument(s) are not known to cSEM: ",
         paste('"', args_diff, '"', sep = "", collapse = ", "), "\n",
         "Please type 'args_default()' to find all valid arguments and their defaults.\nOr check ?args_default.",
         call. = FALSE)
  }

  ## Which arguments need to be changed?
  args_intersect <- intersect(args_dots_names, args_default_names)

  ## Replace all arguments that were changed or explicitly given and keep
  #  the default values for the others

  for(i in args_intersect) {
    x[i] <- .dotdotdot[i]
  }

  return(x)
}
