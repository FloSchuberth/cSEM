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
  data = .object$Information$Data
  
  # calculate Cronbach's alpha & CIs for all constructs
  #  Calculation of the CIs are based on Trinchera et al. (2018)
  # The code for the calculation of the CIs is addpated from the paper Trinchera et al. (2018) 
    alphas=lapply(construct_names,function(x){
    indicators_relevant=colnames(.object$Information$Model$measurement[x,.object$Information$Model$measurement[x,]!=0,drop=FALSE])
    nInd=length(indicators_relevant)
    data_relevant=data[,indicators_relevant]
    nObs=nrow(data_relevant)
    S_relevant=.object$Estimates$Indicator_VCV[indicators_relevant,indicators_relevant]
    
    # sum scores
    SS=Matrix::rowSums(data_relevant)
    sig2SS = var(SS)
    sig4SS = var(SS)^2  
    sig6SS = var(SS)^3
    sig8SS = var(SS)^4
    # Calculation of alpha
    alphaHat=nInd/(nInd - 1)*(1 - sum(diag(cov(data_relevant)))/sig2SS)
    
    # alpha_temp=psych::alpha(S_relevant,delete=FALSE,na.rm=FALSE)
    # alphaHat=alpha_temp$total$raw_alpha
    
    a <- matrix(,nInd,nInd)
    for (pp in 1:nInd){
      for (ll in 1:nInd){
        a [pp,ll] = -diag(S_relevant)[pp]*diag(S_relevant)[ll]+(1/nObs)*sum(((data_relevant[,pp]
                                                       *mean(data_relevant[,pp]))^2)*((data_relevant[,ll]-mean(data_relevant[,ll]))^2))
      }
    }
    
    A <- sum(a)
    b <- matrix(,1,nInd)
    for (pp in 1:nInd){
      b [pp] <- -sig2SS*diag(S_relevant)[pp]+(1/nObs)*sum(((SS*mean(SS))^2)
                                          *((data_relevant[,pp]-mean(data_relevant[,pp]))^2))
    }
    B <- sum(diag(S_relevant))*sum(b)
    CC <- sum(diag(S_relevant))^2*(-sig4SS+(1/nObs)*sum((SS - mean(SS))^4))
    teta.alpha <- (nInd^2/(nInd - 1)^2)*(A/sig4SS-(2*B)/sig6SS+CC/sig8SS)
    seAlpha <-sqrt(teta.alpha/nObs)
    zvalue <- qnorm((1-.alpha/2), mean = 0, sd = 1)
    up <- alphaHat + zvalue*seAlpha
    low <- alphaHat - zvalue*seAlpha
    # return(list(point.estimate = alphaHat, seAlpha = seAlpha,
    # CI = list(low.bound = low, up.bound = up)))
    
    list(point.estimate = alphaHat, seAlpha = seAlpha,
         CI = list(low.bound = low, up.bound = up))
    
    
    
    
    # # Calculating confidence intervals based on Trinchera et al. (2018)
    # 
    # # Can we use standardized indicators for that purpose?
    # X = .object$Information$Data[,indicators_relevant]
    # sumscores=apply(X,1,sum)
    # # Variance and higher-order moments of the of the sum scores
    # varss=var(sumscores)
    # sigma4.S <- varss^2
    # sigma6.S <- varss^3
    # sigma8.S <- varss^4
    # 
    # NrIndicators=length(indicators_relevant)
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
    
    # alphaHat
    
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

# Empirical Bayes correction for loadings suggested by Dijkstra (2018) can be extended to other parameter as well 
# I recommend to extend it to construct correlations instead of the path coefficients and restimate them as we now the bounds of such correlations
# However, this requires that we bootstrap the SEs for the construct correlations.
quasiEmpiricalBayesCorrection <- function(.object,.method=c('median','mean')){
  # Loadings=.object$Estimates$Loading_estimates
  # Loadings[1,2]=2
  # Indmatrix=which(abs(Loadings)>1,arr.ind = F)

  
  L=.object$Estimates$Loading_estimates[.object$Information$Model$measurement!=0]
  

    
    sig_hat=apply(.object$Estimates$Estimates_resample$Estimates1$Loading_estimates$Resampled,2,sd)
    
    # adjustedLoading=c()
    if(.method=='median'){
    for(i in which(abs(L)>1)){
      A=pnorm((-1-L[i])/sig_hat[i])
      B=pnorm((1-L[i])/sig_hat[i])
      # Overwrite the old loadings
      L[i]=L[i] + sig_hat[i] * qnorm(mean(c(A,B)),0,1)
      }
    }

    if(.method=='mean'){
      for(i in which(abs(L)>1)){
        A=pnorm((-1-L[i])/sig_hat[i])
        B=pnorm((1-L[i])/sig_hat[i])
        # Overwrite the old loadings
        L[i]=L[i] + sig_hat[i]/(B-A) *( dnorm((-1-L[i])/sig_hat[i],0,1) - dnorm((1-L[i])/sig_hat[i],0,1))
      }
     }
    
    # Overwrite the old loadings
    .object$Estimates$Loading_estimates[.object$Information$Model$measurement!=0]=L
    
    .object
  }
 

#' Implementation of PLS predict adopted from Shmueli et al. (2016) Table 1
#' 
#' 
predictPLS=function(.object, testDataset){
  
  # Check whether the test dataset contains the same inidcators as the train dataset
  if(sum(!(colnames(.object$Information$Model$measurement)%in% colnames(testDataset)))!=0){
    stop('The same indicators as in the original estimation need to be provided') #Perhaps we do not need to be tha strict as we only need the exogenous indicators
  }
  
  # throw error if model is non-linear. See danks et al how this can be addressed.
  
  # Perform check of the provided dataset
  # Needs to be implemented
  
  # Order the provided dataset
  testData=testDataset[,colnames(.object$Information$Model$measurement)]
  
  # save descriptives of the original unscaled train dataset
  mean_train <- colMeans(.object$Information$Data)
  sd_train <- apply(.object$Information$Data,2,sd)
  W_train <- .object$Estimates$Weight_estimates
  Loadings_train <- .object$Estimates$Loading_estimates
  
  # Path coefficients estimates based on the trainig dataset
  pathtrain <- .object$Estimates$Path_estimates

  # Identifiy exogenous construct in the structural model
  Cons_exo <- rownames(pathtrain)[rowSums(pathtrain)==0]
  Cons_endo <- rownames(pathtrain)[rowSums(pathtrain)!=0]
  
  # Path coefficients of exogenous and endogenous constructs
  B_train      <- .object$Estimates$Path_estimates[Cons_endo, Cons_endo, drop = FALSE]
  Gamma_train  <- .object$Estimates$Path_estimates[Cons_endo, Cons_exo, drop = FALSE]
  
  # Predict scores for the exogenous constructs
  exogscores <- testData%*%t(W_train[exog,,drop = FALSE])
  
  
  # calculate predictions of the endogenous constructs
  endoscores <- exogscores%*%t(Gamma_train) %*% solve(diag(nrow(B_train)) - t(B_train))
  
  
  xhat <- endoscores %*% loadingstrain[endo,,drop = FALSE]
  
  # Denormalize predictions
  xhatrescale= sapply(colnames(xhat),function(x){
    xhat[,x]*sd_train[x]+mean_train[x]
  } )
  
  return(xhatrescale)
  
}





 
  
  
#' Implementation of REBUS-PLS
#' Code obtained by Laura Trinchera.
#' 
REBUS_cluster <- function(.object){
  # if (class(pls)!="plspm") 
    # stop("An object of class 'plspm' was expected")
  # if (any(pls$model[[4]]!="A"))# checking reflective modes DOES it also work with PLSc?
  #   stop("REBUS only works for reflective modes")
  # DOES it also work with PLSc? I chekc whether all constructs are modeled as common factor
  if(any(.object$Information$Model$construct_type == "Common factor")){
    stop("REBUS only works for reflective modes")
  }
  # Is not needed as in cSEM data is always scaled
  # if (!pls$model[[5]])# checking scaled data
  #   stop("REBUS only works with scaled='TRUE'")
  
  DM = .object$Information$Data
  Y.lvs <- .object$Estimates$Construct_scores# recovering LV scores from pls
  loads <- .object$Estimates$Loading_estimates# recovering loadings from pls
  PaCo <- .object$Estimates$Path_estimates# recovering path coeffs from pls
  
  Construct_names <- as.list(rownames(.object$Information$Model$measurement))
  indicator_blocks <- lapply( Construct_names, function(x) {
    names(.object$Information$Model$measurement[x,][.object$Information$Model$measurement[x,]!=0])
  })
  
  IDM <- .object$Information$Model$structural
  
  endo <- rowSums(IDM)
  endo[endo!=0] <- 1  # indicator of endogenous LVs
  
  out.res <- DM# matrix for storing outer resids
  inn.res <- Y.lvs[,endo==1]# matrix for storing inner resids
  # computation of outer residuals
  for (j in 1:length(Construct_names)){
    # q <- which(blocklist==j) 
    X.hat <- Y.lvs[,Construct_names[[j]]] %*% loads[Construct_names[[j]],indicator_blocks[[j]], drop=FALSE]
    out.res[,indicator_blocks[[j]]] <- DM[,indicator_blocks[[j]]] - X.hat[,indicator_blocks[[j]]]# outer residuals
  }
  
  # computationi of inner residuals
  endo_names = rownames(.object$Information$Model$structural[rowSums(.object$Information$Model$structural)!=0,,drop=FALSE])
  predicted=Y.lvs%*%t(PaCo)
  inn.res=Y.lvs[,endo_names,drop=F]- predicted[,endo_names,drop=FALSE]
  
  # if (sum(endo)!=1)# more than 1 endogenous LV
  #   Y.hat <- Y.lvs %*% t(PaCo[endo==1,])        
  # if (sum(endo)==1)# only 1 endogenous LV
  #   Y.hat <- Y.lvs %*% PaCo[endo==1,]        
  # inn.res <- Y.lvs[,endo==1] - Y.hat# inner residuals
  
  # ====================== cluster analysis =====================
  # hierarchical cluster analysis with Ward method using function "hcluster"
  res <- cbind(out.res, inn.res)    
  res.clus <- amap::hcluster(res, method="euclidean", diag=FALSE, upper=FALSE,
                       link="ward", members=NULL, nbproc=2, doubleprecision=TRUE)
  # plot of the dendrogram
  plot(res.clus, main=c("REBUS", "Cluster Dendrogram of Outer and Inner Residuals"),
       hang=-1, cex.main=.9, cex.axis=.5, xlab="Hierarchical Clustering", 
       sub="Ward method", labels=FALSE)
  return(res.clus)
}  
  
REBUS_iterative <-
  function(.object, hclus.res, nk, stop.crit=0.005, iter.max=100)
  {
    # ========================== it.reb function ==========================
    # Performs iterative steps of Response-Based Unit Segmentation  
    # in Partial Least Squares Path Modeling (PLS-PM)        
    # =========================== ARGUMENTS ==============================
    # pls: object of class "plspm"
    # hclus.res: object of class "hclust" obtained from function "res.clus"
    # nk: an integer lasrger than 1 indicating the number of classes 
    # stop.crit: stop criterion number (must be between 0 and 1)
    # iter.max: maximum number of iterations (must be an integer)
    
    # ==================== Checking function arguments ===================
    # if (class(pls)!="plspm") 
    #   stop("argument 'pls' must be an object of class 'plspm'")
    # if (any(pls$model[[4]]!="A"))# checking reflective modes
    #   stop("REBUS only works for reflective modes")
    if(any(.object$Information$Model$construct_type == "Common factor")){
      stop("REBUS only works for reflective modes")
    }
    # if (!pls$model[[5]])# checking scaled data
    #   stop("REBUS only works with scaled='TRUE'")
    if (missing(hclus.res))
      stop("argument 'hclus.res' is missing")
    if (class(hclus.res)!="hclust")
      stop("argument 'hclus.res' must be an object of class 'hclust'")
    if (missing(nk))
      stop("argument 'nk' (number of classes) is missing")
    if (mode(nk)!="numeric" || length(nk)!=1 || 
        nk<=1 || (nk%%1)!=0)
      stop("Invalid number of classes 'nk'. Must be an integer larger than 1")
    if (mode(stop.crit)!="numeric" || length(stop.crit)!=1 || 
        stop.crit<0 || stop.crit>=1)
    {
      warning("Invalid stop criterion 'stop.crit'. Deafult value 0.005 is used")
      stop.crit <- 0.005
    }
    if (mode(iter.max)!="numeric" || length(iter.max)!=1 || 
        iter.max<=1 || (iter.max%%1)!=0)
    {
      warning("Invalid number of maximum iterations 'iter.max'. Deafult value 100 is used")
      iter.max <- 100
    }
    
    # ========================== INPUTS SETTING ==========================
    # IDM <- pls$model[[1]]# Inner Design Matrix
    IDM <- .object$Information$Model$structural# Inner Design Matrix
    lvs.names <- as.list(rownames(.object$Information$Model$measurement))
    Construct_names <- as.list(rownames(.object$Information$Model$measurement))
    blocks <- lapply( Construct_names, function(x) {
      names(.object$Information$Model$measurement[x,][.object$Information$Model$measurement[x,]!=0])
    })
    # blocks <- pls$model[[2]]# cardinality of blocks
    scheme <- pls$model[[3]]# inner weighting scheme
    modes <- pls$model[[4]]# measurement modes
    scaled <- pls$model[[5]]# type of scaling
    plsr <- pls$model[[7]]# pls-regression
    if (plsr) 
      warning("path coefficients will be calculated with OLS regression")
    plsr <- FALSE
    # DM <- pls$data# original data matrix
    DM <- .object$Information$Data# original data matrix     WE USE THE STANDARDIZED ONE
    lvs <- length(Construct_names)# number of LVs
    # lvs.names <- rownames(IDM)# names of LVs
    mvs <- length(colnames(.object$Information$Model$measurement))
    # mvs <- sum(blocks)# number of MVs
    mvs.names <-colnames(.object$Information$Model$measurement)
    # mvs.names <- colnames(DM)
    
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
      blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist) 
    N <- nrow(DM)# number of observations
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1  # indicator of endog LVs
    n.end <- sum(endo)# number of enfod LVs
    # apply the selected scaling
    if(scaled) X=scale(DM) 
    # initial partition
    ini.part <- cutree(hclus.res, nk)# cutting dendrogram in 'nk' clusters
    nclus <- nlevels(factor(ini.part))  # number of clusters
    
    # =============== creating objects to store results ======================
    w.locals <- as.list(1:nclus)# outer.weights
    Y.locals <- as.list(1:nclus)# std latent variables 
    loads.locals <- as.list(1:nclus)# loadings
    path.locals <- as.list(1:nclus)# path coefficients
    R2.locals <- as.list(1:nclus)# R2
    comu.locals <- as.list(1:nclus)# mvs communalities
    outres.locals <- as.list(1:nclus)# communality residuals
    innres.locals <- as.list(1:nclus)# structural residuals
    CM.locals <- matrix(NA,N,nclus)# CM distance of each unit from each model
    
    # ====================== iterative process =====================
    old.clas <- ini.part
    iter.ch <- 0
    repeat 
    {
      # define MV matrix for each initial class
      split.DM <- as.list(1:nclus)
      split.X <- as.list(1:nclus)
      for (k in 1:nclus){
        split.DM[[k]] <- DM[old.clas==k,]       
      }
      # local models computation
      for (k in 1:nclus)
      {   
        mean.k <- apply(split.DM[[k]],2,mean)# local mean
        sd.k <- apply(split.DM[[k]],2,sd)# local std.dev
        # spliting data matrix for each class
        split.X[[k]] <- scale(split.DM[[k]], center=mean.k, scale=sd.k)
        
        # Here we could use the cSEM function to calculate the weights and path coef per class
        
        # calculating outer weights for each class
        out.ws  <- pls.weights(split.X[[k]], IDM, blocks, modes, scheme)
        w.locals[[k]] <- out.ws[[2]]
        # calculating LV scores for each class
        Y.k <- split.X[[k]] %*% out.ws[[2]]
        # calculating path coefficients for each class
        pathmod <- pls.paths(IDM, Y.k, plsr)
        path.locals[[k]] <- pathmod[[2]]
        R2.locals[[k]] <- pathmod[[3]][endo==1]
        # calculating loadings and communalities for each class
        loadcomu <- pls.loads(split.X[[k]], Y.k, blocks)   
        loads.locals[[k]] <- loadcomu[[1]]
        comu.locals[[k]] <- loadcomu[[2]]
        # latent variables for each unit in each local model
        X.k <- scale(DM, center=mean.k, scale=sd.k)
        Y.locals[[k]] <- X.k %*% w.locals[[k]]
        # computation of communality residuals
        out.res <- DM
        for (j in 1:lvs)
        {
          q <- which(blocklist==j) 
          X.hat <- Y.locals[[k]][,j] %*% t(loads.locals[[k]][q])
          out.res[,q] <- (X.k[,q] - X.hat)^2# outer residuals
        }
        outres.locals[[k]] <- out.res
        # computation of inner residuals
        if (sum(endo)!=1)
          Y.hat <- Y.locals[[k]] %*% t(path.locals[[k]][endo==1,])   
        if (sum(endo)==1)
          Y.hat <- Y.locals[[k]] %*% path.locals[[k]][endo==1,]        
        innres.locals[[k]] <- (Y.locals[[k]][,endo==1] - Y.hat)^2
        # computation super normalized residual of outer model
        res.num1 <- outres.locals[[k]] %*% diag(1/comu.locals[[k]],mvs,mvs)
        supres.outer <- rowSums(res.num1) / (sum(rowSums(res.num1))/(N-2))
        # computation of super normalized residual of inner model
        res.num2 <- innres.locals[[k]] %*% diag(1/R2.locals[[k]],n.end,n.end)
        supres.inner <- rowSums(res.num2) / (sum(rowSums(res.num2))/(N-2)) 
        # computation of the CM distance
        CM.locals[,k] <- sqrt(supres.outer * supres.inner)
      }
      # allocating the units to their closest class
      new.clas <- old.clas
      for (i in 1:N)
        new.clas[i] <- which(CM.locals[i,]==min(CM.locals[i,]))[1]
      # checking convergence
      dif.clas <- new.clas - old.clas 
      unit.change <- length(which(dif.clas!=0))
      iter.ch <- iter.ch + 1
      # rate of unit change
      if (unit.change/N < stop.crit || iter.ch==iter.max)
        break
      old.clas <- new.clas
      if (any(table(new.clas)<=5)) 
        stop("Too few units: a class with less than 6 units was detected") 
    }
    
    # ==================== computation of final local models =======================
    DM.cl.rebus <- as.list(1:nclus)    # data matrices for each class
    for (k in 1:nclus)
      DM.cl.rebus[[k]] <- DM[new.clas==k,]
    DM.rebus <- cbind(DM, rebus.class=new.clas)# data matrix with units membership
    redun.locals <- as.list(1:nclus)
    # table of path coefficients for local models
    path.labs <- NULL
    for (j in 1:lvs)
      for (i in j:lvs)
        if (IDM[i,j]==1) 
          path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))
    reb.labs <- paste(rep("Class",nclus),1:nclus,sep=".")
    path.rebus <- matrix(NA,length(path.labs),nclus)
    loads.rebus <- matrix(NA,mvs,nclus)
    qual.rebus <- matrix(NA,sum(lvs+2*n.end+1),nclus)
    for (k in 1:nclus)
    {               
      out.ws <- pls.weights(scale(DM.cl.rebus[[k]]), IDM, blocks, modes, scheme)
      Y.k <- scale(DM.cl.rebus[[k]]) %*% out.ws[[2]]
      pathmod <- pls.paths(IDM, Y.k, plsr)
      path.locals[[k]] <- pathmod[[2]]
      R2.locals[[k]] <- pathmod[[3]][endo==1]
      loadcomu <- pls.loads(DM.cl.rebus[[k]], Y.k, blocks)    
      loads.locals[[k]] <- loadcomu[[1]]
      comu.locals[[k]] <- loadcomu[[2]]
      redun <- rep(0,mvs)
      aux <- 0
      for (j in 1:lvs)
        if (endo[j]==1)
        {
          aux <- aux + 1
          redun[blocklist==j] <- comu.locals[[k]][blocklist==j] * R2.locals[[k]][aux]
        }
      redun.locals[[k]] <- redun
      
      path.rebus[,k] <- as.vector(path.locals[[k]][IDM==1])# path coeffs for local models        
      loads.rebus[,k] <- loads.locals[[k]]# loadings for local models  
      # table of quality indexes
      comu.aveg <- rep(NA,lvs) 
      redun.aveg <- rep(NA,n.end)
      R2.aux <- rep(NA,n.end)
      aux <- 0
      for (j in 1:lvs)
      {
        comu.aveg[j] <- mean(comu.locals[[k]][which(blocklist==j)]) 
        if (endo[j]==1) 
        {
          aux <- aux + 1      
          redun.aveg[aux] <- mean(redun.locals[[k]][which(blocklist==j)])
          R2.aux[aux] <- R2.locals[[k]][aux]
        }
      }
      qual.rebus[1:lvs,k] <- comu.aveg
      qual.rebus[(lvs+1):(lvs+n.end),k] <- redun.aveg
      qual.rebus[(lvs+n.end+1):(lvs+2*n.end),k] <- R2.aux
      qual.rebus[(lvs+2*n.end+1),k] <- sqrt(mean(comu.aveg)*mean(R2.aux))
    }
    gqi <- round(GQI(pls, new.clas), 4)
    dimnames(path.rebus) <- list(path.labs, reb.labs)
    dimnames(loads.rebus) <- list(mvs.names, reb.labs)
    v1 <- paste(rep("Com",lvs),lvs.names,sep=".")
    v2 <- paste(rep("Red",n.end),lvs.names[endo==1],sep=".")
    v3 <- paste(rep("R2",n.end),lvs.names[endo==1],sep=".")
    dimnames(qual.rebus) <- list(c(v1,v2,v3,"GoF"), reb.labs)
    qual.rebus <- round(qual.rebus,4)
    aux <- list(lvs,n.end,unit.change/N,stop.crit,iter.max,iter.ch,gqi)
    res <- list(path.coef=path.rebus, loadings=loads.rebus, quality=qual.rebus,
                segments=new.clas, origdata.clas=DM.rebus, aux=aux)
    class(res) <- "rebus"
    return(res)
  }
