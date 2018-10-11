### Use this file to add new functions whose name you have not yet decided on
# or when it is unclear where it belongs.


#' AVE
#'
#' Compute the average variance extracted (AVE) based on the standardized loadings.
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
#' 
#' Fornell, C. & Larcker, D. F. (1981). Evaluating structural equation models with unobservable variables and measurement error. \emph{Journal of Marketing Research}, 18, 39–50.
#'
#'
#' @export
#'
AVE=function(.object=args_default()$.object,
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
  
  # By default, for composites the CR is set to 1
  if(.only_common_factors){
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite"])
    AVEs[co_names]=NULL
  }
  
  return(AVEs) 
}



#' CR
#'
#' Compute the composite reliability (CR) based on standardized loading, see Raykov (1997)
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
#' Raykov, T. (1997). Estimation of composite reliability for congeneric measures. \emph{Applied Psychological Measurement}, 21(2), 173–184.
#'
#' @export
#'
CR=function(.object=args_default()$.object,
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
  
  # By default, for composites the CR is set to 1
  if(.only_common_factors){
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite"])
    CRs[co_names]=1
  }
  
  return(CRs)
}

#' Cronbach_alpha
#'
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
#' 
#' Cronbach, L. J.(1951). Coefficinet alpha and the internal structure of tests. \emph{Psychometrika}, 16(3), 297–334.
#'
#' @export
#'
Cronbach_alpha=function(.object=args_default()$.object,
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
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite"])
    alphas[co_names]=NULL
  }
  
  return(alphas)
}

#' Fornell_Larcker
#'
#' Computes the Fornell-Larcker criterion
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
#' Fornell, C. & Larcker, D. F. (1981). Evaluating structural equation models with unobservable variables and measurement error. \emph{Journal of Marketing Research}, 18, 39–50.
#'
#' @export
#'
Fornell_Larcker=function(.object=args_default()$.object,
                                 .only_common_factors=args_default()$.only_common_factors){
  
  # Calculate the average variance extracted
  average_variance_extracted=AVE(.object,.only_common_factors=FALSE)
  
  .object
  
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
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite"])
    alphas[co_names]=NULL
  }
  
  return(alphas)
}

