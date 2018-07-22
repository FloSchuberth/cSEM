#' Calculate factor and composite loadings and cross-loadings
#'
#' Calculates factor loadings (for constructs modeled as common factors) and
#' composite loadings (for constructs modeled as composites).
#'
#' Here more details
#'
#' @inheritParams csem_arguments
#'
#' @return The (J x K) matrix of loadings and "cross loadings" (attenuated if required).
#'
calculateLoadings <- function(
  .S                  = NULL,
  .W                  = NULL,
  .Q                  = NULL,
  .csem_model         = NULL,
  .disattenuate       = NULL,
  .modes              = NULL
  # .correction_factors = NULL
) {

  ## Matrix of loadings and cross-loadings
  Lambda <- .W %*% .S

  if(.disattenuate) {
    
    ## Get names of constructs modeled as composites
    names_c  <- names(.csem_model$construct_type[.csem_model$construct_type == "Composite"])
    ## Get names of constructs modeled as common factors
    names_cf <- setdiff(rownames(.csem_model$structural), names_c)
    ## Get names of the common factors whose weights were estimated with "ModeA"
    names_modeA <- intersect(names(.modes[.modes == "ModeA"]), names_cf)
    ## Get names of the common factors whose weights were estimated with "ModeB"
    names_modeB <- intersect(names(.modes[.modes == "ModeB"]), names_cf)
    
    if(length(names_cf) > 0) {
      if(length(names_modeA) > 0) {
        ## Disattenuate loadings and cross-loadings ------------------------------
        for(i in names_modeA) {
          temp  <- .W[i, ] # becomes a vector!
          temp1 <- temp[which(temp != 0)]
          temp2 <- temp[which(temp == 0)]
          
          Lambda[i, names(temp1)] <- .Q[i] * temp1 / c(t(temp1) %*% temp1)
          Lambda[i, names(temp2)] <- Lambda[i, names(temp2)] / .Q[i]
          
        }
      }
      if(length(names_modeB) > 0) {
        stop("Variable(s): ", paste0("`", names_modeB, "`", collapse = ", "), 
             " were estimated using Mode B\n",
             "Currently correction for attenuation for constructs modeled as",
             " common factors is only possible for Mode A.",
             call. = FALSE)
      }
    } else {# all constructs are modeled as composites
      return(Lambda)
    }
  } # END .disattentuate == TRUE
  return(Lambda)
}
