### Use this file to add new functions whose name you have not yet decided on
# or when it is unclear where it belongs.


AVE <- function(.object) {
  UseMethod("AVE")
}


#' AVE
#'
#' Computes the average variance extracted \insertCite{Raykov1997}{cSEM} based on the standardized loadings.
#'
#' @usage AVE(
#'  .object              = args_default()$.object,
#'  .only_common_factors = args_default()$.only_common_factors
#' )
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
#' @references 
#' \insertAllCited{}
#'
#'
#' @export
#'
AVE.cSEMResults_default=function(.object=args_default()$.object,
             .only_common_factors=args_default()$.only_common_factors){
  
  construct_names=names(.object$Information$Model$construct_type)
  
  # Extract loadings
  L=.object$Estimates$Loading_estimates

  AVEs=sapply(construct_names, function(x){
    lam=c(L[x,L[x,]!=0])
    ave=sum(lam^2)/(sum(lam^2)+sum(1-lam^2))
    ave
    
  })
  
  names(AVEs)=construct_names
  
  # By default, for composites the AVEs is hidden
  if(.only_common_factors){
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite",drop=F])
    AVEs=AVEs[setdiff(construct_names,co_names),drop=F]
  }
  
  return(AVEs) 
}

AVE.cSEMResults_multi <- function(.object) {
  
  lapply(.object, AVE.cSEMResults_default)
}







CR <- function(.object) {
  UseMethod("CR")
}
#' CR
#'
#' Computes the composite reliability \insertCite{Raykov1997}{cSEM} based on standardized loading.
#'
#' @usage CR(
#'  .object              = args_default()$.object,
#'  .only_common_factors = args_default()$.only_common_factors
#' )
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
#' @references 
#' 
#' \insertAllCited{}
#'
#' @export
#'
CR.cSEMResults_default=function(.object=args_default()$.object,
            .only_common_factors=args_default()$.only_common_factors){
  
  
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
  
  # By default, for composites the CR is set to 1, otherwise they are returned
  if(.only_common_factors){
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite",drop=F])
    CRs[intersect(construct_names,co_names),drop=F]=1
  }
  
  return(CRs)
}


CR.cSEMResults_multi <- function(.object) {
  
  lapply(.object, CR.cSEMResults_default)
}

#' Cronbach_alpha
#'
#' Computes Cronbach's alpha \insertCite{Cronbach1951}{cSEM} based on the correlation matrix 
#'
#' @usage Cronbach_alpha(
#'  .object              = args_default()$.object,
#'  .only_common_factors = args_default()$.only_common_factors
#' )
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
#' @references 
#' \insertAllCited{}
#'
#' @export
#'
#'
#'

Cronbach_alpha.cSEMResults_default=function(.object=args_default()$.object,
            .only_common_factors=args_default()$.only_common_factors){
  construct_names=names(.object$Information$Model$construct_type)
  
  # calculate Cronbach's alpah for all constructs
  alphas=sapply(construct_names,function(x){
    relevant_indicators=colnames(.object$Information$Model$measurement[x,.object$Information$Model$measurement[x,]!=0,drop=FALSE])
    S_relevant=.object$Estimates$Indicator_VCV[relevant_indicators,relevant_indicators]
    alpha_temp=psych::alpha(S_relevant,delete=FALSE,na.rm=FALSE)
    alpha_temp$total$raw_alpha
    })
  
  names(alphas)=construct_names
  
  # By default, for composites the CR is set to 1
  if(.only_common_factors){
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite",drop=F])
    alphas[intersect(construct_names,co_names),drop=F]=NULL
  }
  
  return(alphas)
}

#' Fornell_Larcker
#'
#' Computes the Fornell-Larcker criterion \insertCite{Fornell1981}{cSEM}.
#'
#' @usage Fornell_Larcker(
#'  .object              = args_default()$.object,
#'  .only_common_factors = args_default()$.only_common_factors
#' )
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @examples
#' \dontrun{
#' # still to implement
#' }
#'@references
#'
#' \insertAllCited{}
#'
#' @export
#'
Fornell_Larcker=function(.object=args_default()$.object,
                                 .only_common_factors=args_default()$.only_common_factors){

  # Calculate the average variance extracted
  average_variance_extracted=AVE(.object,.only_common_factors=FALSE)

  .object

  construct_names=names(.object$Information$Model$construct_type)

  construct_vcv_sq=.object$Estimates$Construct_VCV^2
  diag(construct_vcv_sq)=average_variance_extracted[colnames(construct_vcv_sq)]
stop('Not implemented yet')

}



#' Calculate effect size
#'
#' (TODO)
#'
#' @usage calculateEffectSize(
#'  .object              = args_default()$.object,
#' )
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
calculateEffectSize <- function(.object) {
  
  # Get relevenat quantities
  H <- .object$Estimates$Construct_scores
  Q <- sqrt(.object$Estimates$Construct_reliabilities)
  P <- .object$Estimates$Construct_VCV
  csem_model  <- .object$Information$Model
  normality   <- .object$Information$Arguments$.normality
  approach_nl <- .object$Information$Arguments$.approach_nl
  
  s <- csem_model$structural

  vars_endo <- rownames(s)[rowSums(s) != 0]
  outer_out <- lapply(vars_endo, function(x) {
    
    # get colnames
    indep_vars <- colnames(s[x , s[x, ] != 0, drop = FALSE])
    
    inner_out <- lapply(indep_vars, function(i) {
      # update csem_model
      xx <- csem_model
      
      xx$structural[x, i] <- 0 
      
      out <- cSEM:::estimatePathOLS(
        .H = H,
        .Q = Q,
        .P = P,
        .csem_model = xx,
        .normality = normality,
        .approach_nl = approach_nl
      )
      # calculate 
      r2_excluded <- out$R2[x]
      r2_included <- .object$Estimates$R2[x]
      
      effect_size <- (r2_included - r2_excluded)/(1 - r2_included)
      list("r2_ex" = r2_excluded, "r2_in" = r2_included, "eff_size" = effect_size)
    })
    names(inner_out) <- indep_vars
    inner_out
  })
  names(outer_out) <- vars_endo
  outer_out
}


GoF.default=function(.object){
  
  
  m         <- .object$Information$Model$structural
  # Endogenous constructs
  vars_endo <- rownames(m)[rowSums(m) != 0]

  # Collect names of common factors
  ct=.object$Information$Model$construct_type
  cf=names(ct[ct == 'Common factor'])
  
  
  # Select loadings of the common factors
  L=.object$Estimates$Loading_estimates[cf,]
    L=L[L!=0]
  
    # Calculate GoF. 
    # It is defined as the mean of the R^2s of the structural model times
    # the variance in the indicators that is explained by the construct.
    # For the latter, only constructs modeled as common factors are considered
    # as they explain their indicators in contrast to the composite.
    sqrt(mean(L^2) * mean(.object$Estimates$R2))

}
