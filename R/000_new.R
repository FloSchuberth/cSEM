### Use this file to add new functions whose name you have not yet decided on
# or when it is unclear where it belongs.


# Calculates the Heterotrait-Monotrait factor correlations, see Henseler et al. (2015)


HTMT = function(.object){
  S=.object$Estimates$Indicator_VCV
  
  # HTMT only for common factors
  cf_names=names(.object$Information$Construct_types[.object$Information$Construct_types=="Common factor"])
  
  stop('Not yet implemented') 
  
}


#' `cSEMResults` method for `summary()`
#'
#' The [cSEMResults] method for the generic function [summary()]. 
#' 
#' Computes a summary of the results obtained from running [csem], [cca], or
#' [foreman].
#'
#' @usage summary(.object, .what = NULL)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#' 
#' @inherit csem_summary return
#'
#' @export
#'
summary.cSEMResults <- function(.object, .what = NULL) {
  
  ## Structure loadings output
  temp <- x$Estimates$Loading_estimates
  names_loadings <- paste0(rep(rownames(temp), times = apply(temp, 1, function(x) sum(x != 0))),
                           " =~ ", colnames(temp))
  loading_estimates <- data.frame("Loading" = names_loadings,
                                  "Estimate" = unlist(t(temp)[t(temp) != 0 ]),
                                  stringsAsFactors = FALSE)
  
  ## Structure weights output
  temp <- x$Estimates$Weight_estimates
  names_weights <- paste0(rep(rownames(temp), times = apply(temp, 1, function(x) sum(x != 0))),
                          " -- ", colnames(temp))
  weight_estimates <- data.frame("Weights" = names_weights,
                                 "Estimate" = unlist(t(temp)[t(temp) != 0 ]),
                                 stringsAsFactors = FALSE)
  
  ## Create summary list
  summary_out <- list(
    "Construct_types"        = x$Meta_information$Construct_types,
    "Correction_factors"     = x$Estimates$Correction_factors,
    "Loading_estimates"      = loading_estimates,
    "Names_endogenous_var"   = x$Meta_information$modellist$vars_endo,
    "Number_of_observations" = x$Meta_information$Number_of_observations,
    "Path_estimates"         = x$Estimates$Path_estimates,
    "Path_estimator"         = x$Meta_information$Path_approach,
    "PLS_modes"              = x$Meta_information$PLS_Modes,
    "PLS_weight_scheme_inner"= x$Meta_information$PLS_Inner_Weightning_scheme,
    "Weight_estimates"       = weight_estimates,
    "Weight_estimator"       = x$Meta_information$Weight_approach
  )
  
  ## Set class for printing and return
  class(summary_out) <- "cSEMResultssummary"
  return(summary_out)
}

#' Effects
#'
#' Compute direct, indirect and total effects.
#'
#' @usage effects(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @examples
#' \dontrun{
#' # still to implement
#' }
#'
#' @export
#'
effects.cSEMResults <- function(.object) {
  # Implementation is inspired by the matrixpls package licensed under GPL-3
  
  ## Endogenous (lhs) variables
  vars_endo <- .object$Information$Model$vars_endo
  
  ## Matrix of direct effects:
  direct <- .object$Estimates$Path_estimates
  
  ## Matrix of total total effects: B = direct
  # Note: eta = B x eta + zeta
  #       (I - B)*eta = zeta
  
  B_star <- diag(nrow(direct)) - direct
  total <- solve(B_star) - diag(nrow(direct))
  
  ## Matrix of indirect effects:
  indirect <- total - direct
  
  list(direct = direct[vars_endo, , drop = FALSE], 
       indirect = indirect[vars_endo, , drop = FALSE], 
       total = total[vars_endo, , drop = FALSE])
}


## this is some test change