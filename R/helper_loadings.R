#' Calculate factor and composite loadings
#'
#' Calculates factor loadings (for constructs modeled as common factors) and
#' composite loadings (for constructs modeled as composites).
#'
#' Here more details
#'
#' @inheritParams csem_arguments
#'
#' @return The (J x K) matrix of loadings (disattenuated if requested).
#'
calculateLoadings <- function(
  .S                  = NULL,
  .W                  = NULL,
  .csem_model         = NULL,
  .disattenuate       = NULL,
  .modes              = NULL,
  .correction_factors = NULL
) {

  x <- .W %*% .S * .csem_model$measurement

  if(.disattenuate) {

    ## Get names of constructs modeled as composites
    names_c  <- .csem_model$construct_type[.csem_model$construct_type$Type == "Composite", ]$Name
    ## Get names of constructs modeled as common factors
    names_cf <- setdiff(rownames(.csem_model$structural), names_c)
    ## Get names of the common factors whose weights where estimated with "ModeA"
    names_modeA <- intersect(names(.modes[.modes == "ModeA"]), names_cf)
    ## Get names of the common factors whose weights where estimated with "ModeB"
    names_modeB <- intersect(names(.modes[.modes == "ModeB"]), names_cf)

    ##
    if(length(names_modeA) > 0) {

      x_modesA <- t(t(.W[names_modeA,, drop = FALSE]) %*%
                        diag(.correction_factors[names_modeA]))

      if(length(names_modeB) > 0) {
        stop("Variable(s): ", paste0("`", names_modeB, "`", collapse = ", "), 
             " where estimated using Mode B\n",
             "Currently correction for attenuation for constructs modeled as",
             " common factors is only possible for Mode A.",
             call. = FALSE)
      }

      x[names_modeA, ] <- x_modesA
    } else {

      stop("Variable(s): ", paste0("`", names_modeB, "`", collapse = ", "), 
           " where estimated using Mode B\n",
           "Currently correction for attenuation for constructs modeled as",
           " common factors is only possible for Mode A.",
           call. = FALSE)
    }
  }
  return(x)
}
