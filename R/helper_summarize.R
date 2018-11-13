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
#' @keywords internal
#'

effects <- function(.object) {
  # Implementation is inspired by the matrixpls package licensed under GPL-3
  
  m         <- .object$Information$Model$structural
  ## Endogenous (lhs) variables
  vars_endo <- rownames(m)[rowSums(m) != 0]
  
  ## Matrix of direct effects:
  direct <- .object$Estimates$Path_estimates
  
  ## Matrix of total total effects: B = direct
  # Note: eta = B x eta + zeta
  #       (I - B)*eta = zeta
  
  B_star <- diag(nrow(direct)) - direct
  total <- solve(B_star) - diag(nrow(direct))
  
  ## Matrix of indirect effects:
  indirect <- total - direct
  
  # MAtrix containing the variance accounted for (VAR)
  VAR <- indirect/total
  
  list(direct = direct[vars_endo, , drop = FALSE], 
       indirect = indirect[vars_endo, , drop = FALSE], 
       total = total[vars_endo, , drop = FALSE],
       var = VAR[vars_endo, , drop = FALSE])
}