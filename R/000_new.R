### Use this file to add new functions whose name you have not yet decided on
# or when it is unclear where it belongs.

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

AVE <- function(.object, 
                .only_common_factors=args_default()$.only_common_factors) {
  UseMethod("AVE")
}

#' @describeIn AVE (TODO)
#' @export
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

#' @describeIn AVE (TODO)
#' @export
AVE.cSEMResults_multi <- function(.object=args_default()$.object,
                                  .only_common_factors=args_default()$.only_common_factors) {
  
  lapply(.object, AVE.cSEMResults_default)
}


#' @describeIn AVE (TODO)
#' @export
AVE.cSEMResults_2ndorder <- function(.object=args_default()$.object,
                                     .only_common_factors=args_default()$.only_common_factors) {
  
  stop('Not implemented yet.')
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

CR <- function(.object=args_default()$.object,
               .only_common_factors=args_default()$.only_common_factors) {
  UseMethod("CR")
}

#' @describeIn CR (TODO)
#' @export
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
  
  # By default, for composites the CR is omitted, otherwise they are returned
  if(.only_common_factors){
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite",drop=F])
    CRs=CRs[setdiff(construct_names,co_names),drop=F]
  }
  
  return(CRs)
}

#' @describeIn CR (TODO)
#' @export
CR.cSEMResults_multi <- function(.object=args_default()$.object,
                                 .only_common_factors=args_default()$.only_common_factors) {
  
  lapply(.object, CR.cSEMResults_default)
}


#' @describeIn CR (TODO)
#' @export
CR.cSEMResults_2ndorder <- function(.object=args_default()$.object,
                                    .only_common_factors=args_default()$.only_common_factors) {
  
  stop('Not implemented yet')
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
Cronbach_alpha <- function(.object=args_default()$.object,
                           .only_common_factors=args_default()$.only_common_factors,
                           .alpha=args_default()$.alpha) {
  UseMethod("Cronbach_alpha")
}


#' @describeIn Cronbach_alpha (TODO)
#' @export
Cronbach_alpha.cSEMResults_default=function(.object=args_default()$.object,
            .only_common_factors=args_default()$.only_common_factors,
            .alpha=args_default()$.alpha){
  
  construct_names=names(.object$Information$Model$construct_type)
  
  # calculate Cronbach's alpah for all constructs
  alphas=sapply(construct_names,function(x){
    relevant_indicators=colnames(.object$Information$Model$measurement[x,.object$Information$Model$measurement[x,]!=0,drop=FALSE])
    S_relevant=.object$Estimates$Indicator_VCV[relevant_indicators,relevant_indicators]
    alpha_temp=psych::alpha(S_relevant,delete=FALSE,na.rm=FALSE)
    alphaHat=alpha_temp$total$raw_alpha
    
    # # Calculating confidence intervals based on Trinchera et al. (2018)
    # 
    # # Can we use standardized indicators for that purpose?
    # X = .object$Information$Data[,relevant_indicators]
    # sumscores=apply(X,1,sum)
    # # Variance and higher-order moments of the of the sum scores
    # varss=var(sumscores)
    # sigma4.S <- varss^2
    # sigma6.S <- varss^3
    # sigma8.S <- varss^4
    # 
    # NrIndicators=length(relevant_indicators)
    # nObs = nrow(X)
    # 
    # # alpha (equals the one from the psych package)
    # NrIndicators/(NrIndicators - 1)*(1 - sum(diag(cov(X)))/varss)
    # 
    # a <- matrix(,NrIndicators,NrIndicators)
    # for (pp in 1:NrIndicators){
    #   for (ll in 1:NrIndicators){
    #     a [pp,ll] = -diag(cov(X))[pp]*diag(cov(X))[ll]+(1/nObs)*sum(((X[,pp]
    #                                                    *mean(X[,pp]))^2)*((X[,ll]-mean(X[,ll]))^2))
    #   }
    # }
    # 
    # A <- sum(a)
    # b <- matrix(,1,NrIndicators)
    # for (pp in 1:NrIndicators){
    #   b [pp] <-
    #   -varss*diag(cov(X))[pp]+(1/nObs)*sum(((sumscores*mean(sumscores))^2)
    #                                       *((X[,pp]-mean(X[,pp]))^2))
    # }
    # B <- sum(diag(var(X)))*sum(b)
    # 
    # CC <- sum(diag(var(X)))^2*(-sigma4.S+(1/nObs)*sum((sumscores - mean(sumscores))^4))
    # teta.alpha <- (NrIndicators^2/(NrIndicators - 1)^2)*(A/sigma4.S-(2*B)/sigma6.S+CC/sigma8.S)
    # seAlpha <-sqrt(teta.alpha/nObs)
    # zvalue <- qnorm((1-.alpha/2), mean = 0, sd = 1)
    # up <- (alphaHat + zvalue*seAlpha)
    # low <- (alphaHat - zvalue*seAlpha)
    # return(list(point.estimate = alphaHat, seAlpha = seAlpha,
    # CI = list(low.bound = low, up.bound = up)))    
    
    alphaHat
    
    })# End of sapply
  
  names(alphas)=construct_names
  
  # By default, Cronbach's alpha is omitted for composite, otherwise it is returned 
  if(.only_common_factors){
    co_names=names(.object$Information$Model$construct_type[.object$Information$Model$construct_type=="Composite",drop=F])
    alphas=alphas[setdiff(construct_names,co_names),drop=F]
  }
  
  
  # Implement closed-form CIs see Trinchera et al. (2018)
  # Rcode adopted from Trinchera et al. (2018)
    # MyAlpha.CI <− function (X, sig = 0.05){
    # p <− ncol(X)
    # n <− nrow(X)
    # S <−apply(X,1,sum)
    # sigma.S <− var(S)
    # alphaHat <− (p/(p − 1))*(1−sum(diag(cov(X)))/sigma.S)
    # sigma4.S <− sigma.S^2
    # sigma6.S <− sigma.S^3
    # sigma8.S <− sigma.S^4
    # a <− matrix(,p,p)
    # for (pp in 1:p){
    #   for (ll in 1:p){
    #     a [pp,ll] <−
    #     −diag(cov(X))[pp]*diag(cov(X))[ll]+(1/n)*sum(((X[,pp]
    #                                                    *mean(X[,pp]))^2)*((X[,ll]−mean(X[,ll]))^2))
    #   }
    # }
    # A <− sum(a)
    # b <− matrix(,1,p)
    # for (pp in 1:p){
    #   b [pp] <−
    #   −sigma.S*diag(cov(X))[pp]+(1/n)*sum(((S*mean(S))^2)
    #                                       *((X[,pp]−mean(X[,pp]))^2))
    # }
    # B <− sum(diag(var(X)))*sum(b)
    # CC <− sum(diag(var(X)))^2*(−sigma4.S+(1/n)*sum((S
    #                                                 −mean(S))^4))
    # teta.alpha <− (p^2/(p − 1)^2)*(A/sigma4.S−(2*B)/sigma6.
    #                                S+CC/sigma8.S)
    # seAlpha<−sqrt(teta.alpha/n)
    # zvalue <− qnorm((1−sig/2), mean = 0, sd = 1)
    # up <− (alphaHat + zvalue*seAlpha)
    # low <− (alphaHat − zvalue*seAlpha)
    # return(list(point.estimate = alphaHat, seAlpha = seAlpha,
                # CI = list(low.bound = low, up.bound = up)))
  
  
  return(alphas)
}

#' @describeIn Cronbach_alpha (TODO)
#' @export
Cronbach_alpha.cSEMResults_multi <- function(.object=args_default()$.object,
                                             .only_common_factors=args_default()$.only_common_factors,
                                             .alpha=args_default()$.alpha) {
  
  lapply(.object, Cronbach_alpha.cSEMResults_default)
}

#' @describeIn Cronbach_alpha (TODO)
#' @export
Cronbach_alpha.cSEMResults_2ndorder <- function(.object=args_default()$.object,
                                                .only_common_factors=args_default()$.only_common_factors) {
  
  stop('Not implemented yet.')
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

calculateEffectSize <- function(.object=args_default()$.object) {
  UseMethod("calculateEffectSize")
}


#' @describeIn calculateEffectSize (TODO)
#' @export
calculateEffectSize.cSEMResults_default <- function(.object=args_default()$.object) {
  
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
      model_temp <- csem_model
      
      model_temp$structural[x, i] <- 0 
      
      out <- cSEM:::estimatePathOLS(
        .H = H,
        .Q = Q,
        .P = P,
        .csem_model = model_temp,
        .normality = normality,
        .approach_nl = approach_nl
      )
      # calculate 
      r2_excluded <- out$R2[x]
      r2_included <- .object$Estimates$R2[x]
      
      effect_size <- (r2_included - r2_excluded)/(1 - r2_included)
      # list("r2_ex" = r2_excluded, "r2_in" = r2_included, "eff_size" = effect_size)
      list("effect_size" = effect_size)
    })
    names(inner_out) <- indep_vars
    inner_out
  })
  names(outer_out) <- vars_endo
  outer_out
}

#' @describeIn calculateEffectSize (TODO)
#' @export
calculateEffectSize.cSEMResults_multi <- function(.object=args_default()$.object) {
  
  lapply(.object, calculateEffectSize.cSEMResults_default)
}

#' @describeIn calculateEffectSize (TODO)
#' @export
calculateEffectSize.cSEMResults_2ndorder <- function(.object=args_default()$.object) {
  
  stop('Not implemented yet.')
}






#' Goodness of fit
#'
#' (TODO)
#'
#' @usage GoF(
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
#'
GoF <- function(.object=args_default()$.object) {
  UseMethod("GoF")
}


#' @describeIn GoF (TODO)
#' @export
GoF.cSEMResults_default=function(.object){
  
  
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

#' @describeIn GoF (TODO)
#' @export
GoF.cSEMResults_multi <- function(.object=args_default()$.object) {
  
  lapply(.object, GoF.cSEMResults_default)
}

#' @describeIn GoF (TODO)
#' @export
GoF.cSEMResults_2ndorder <- function(.object=args_default()$.object) {
  
  stop('Not implemented yet.')
}



#' Variance Inflation Factor for weight estimates obtained by Mode B
#'
#' (TODO)
#'
#' @usage calculateVIFPLS(
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
#'
calculateVIFPLS <- function(.object=args_default()$.object) {
  UseMethod("calculateVIFPLS")
}


#' @describeIn calculateVIFPLS (TODO)
#' @export
calculateVIFPLS.cSEMResults_default=function(.object=args_default()$.object){
  
  if(.object$Information$Arguments$.approach_weights!='PLS-PM'){
    stop('VIF is only calculate when weights are obtained by PLS-PM') }
  
  # Calculation of the VIF makes only sense when Mode B is used 
  # Extract composite names that are build by Mode B 
  ModesUsed=.object$Information$Weight_info$Modes
  ModeB=names(which(ModesUsed=='ModeB'))
  
  
  measurement=.object$Information$Model$measurement
  VIF=list()
  for(blocks in ModeB){
    all = names(which(measurement[blocks, ] == 1))
    if (length(all) != 1) {# if a none single-indicator block
      VIF[[blocks]]=lapply(1:length(all), function(x) {
        indep = all[-x]
        dep = all[x]
        # paste0(dep,'~',paste0(indep,sep= '+ '))
        form = paste0(dep, '~', paste0(indep, collapse = '+ '))
        
        res = lm(form, data = as.data.frame(.object$Information$Data))
        R2 = cor(res$fitted.values, .object$Information$Data[, dep]) ^ 2
        names(R2) = dep
        1/(1-R2)
      })
      
    } else{#if single-indicator block
      VIF[[blocks]]=NaN
    }
  }
  return(VIF)
}

#' @describeIn calculateVIFPLS (TODO)
#' @export
calculateVIFPLS.cSEMResults_multi <- function(.object=args_default()$.object) {
  
  lapply(.object, calculateVIFPLS.cSEMResults_default)
}

#' @describeIn calculateVIFPLS (TODO)
#' @export
calculateVIFPLS.cSEMResults_2ndorder <- function(.object=args_default()$.object) {
  
  stop('Not implemented yet.')
}
  