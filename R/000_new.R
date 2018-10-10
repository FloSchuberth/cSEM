### Use this file to add new functions whose name you have not yet decided on
# or when it is unclear where it belongs.


# Calculates the average variance extracted (AVE)
# The current calculation requires standardized loadings
AVE=function(.object,.only_common_factors=TRUE){
  construct_names=names(.object$Information$Model$construct_type)
  
  # Extract loadings
  L=.object$Estimates$Loading_estimates

  AVEs=sapply(construct_names, function(x){
    lam=c(L[x,L[x,]!=0])
    ave=sum(lam^2)/(sum(lam^2)+sum(1-lam^2))
    ave
    
  })
  
  names(AVEs)=construct_names
  
  # By default, for composites the CR is set to 1
  if(.only_common_factors){
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite"])
    AVEs[co_names]=NULL
  }
  
  return(CRs)
  
  
}


# Calculates the composite reliability, see Raykov (1997).
# The current calculation requires standardized loadings, which is not an issue yet.
# In case of hierarchical models, we need to apply this function on both stages
CR=function(.object, .only_common_factors=TRUE){
  construct_names=names(.object$Information$Model$construct_type)
  
  # Extract loadings
  L=.object$Estimates$Loading_estimates
  
  
  # Calculate CR for all constructs
  CRs=sapply(construct_names, function(x){
    lam=c(L[x,L[x,]!=0])
    cr=sum(lam)^2/(sum(lam)^2+sum(1-lam^2))
    cr
    
  })
  
  names(CRs)=construct_names
  
  # By default, for composites the CR is set to 1
  if(.only_common_factors){
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite"])
    CRs[co_names]=1
  }
  
  return(CRs)
}


# Calculates the Heterotrait-Monotrait factor correlations, see Henseler et al. (2015)
HTMT = function(.object){
  S=.object$Estimates$Indicator_VCV
  
  # HTMT only for common factors
  cf_names=names(.object$Information$Construct_types[.object$Information$Construct_types=="Common factor"])
  
  stop('Not yet implemented') 
  
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
