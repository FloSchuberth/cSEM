### Use this file to add new functions whose name you have not yet decided on
# or when it is unclear where it belongs.

# Empirical Bayes correction for loadings suggested by Dijkstra (2018) can be 
# extended to other parameter as well .
# I recommend to extend it to construct correlations instead of the path 
# coefficients and restimate them as we now the bounds of such correlations
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
 
predict=function(.object, testDataset){
  
  # Check whether the test dataset contains the same inidcators as the train dataset
  if(sum(!(colnames(.object$Information$Model$measurement)%in% colnames(testDataset)))!=0){
    stop('The same indicators as in the original estimation need to be provided') #Perhaps we do not need to be tha strict as we only need the exogenous indicators
  }
  
  # Return error if model is non-linear. See danks et al how this can be addressed.
  if(.object$Information$Model$model_type != 'Linear'){
    stop2('Currenlty, predictPLS works only for linear models.')
  }
  
  if(.object$Information$Type_of_indicator_correlation != 'Pearson'){
    stop2('Currently, predictPLS works only in combination with Pearson correlation.')
  }
  
  # Perform check of the provided dataset
  # Needs to be implemented
  
  # Order the provided dataset
  testData <- testDataset[,colnames(.object$Information$Model$measurement)]
  
  
  # save descriptives of the original unscaled train dataset
  trainData <- .object$Information$Arguments$.data
  mean_train <- colMeans(trainData)
  sd_train <- apply(trainData,2,sd)
  
  # Scale the test dataset with the descriptives of taining dataset
  testDatascale <- sapply(colnames(testData),function(x){
    (testData[,x]-mean_train[x])/sd_train[x]
  })
  
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
  exogscores <- testDatascale%*%t(W_train[Cons_exo,,drop = FALSE])
  
  
  # calculate predictions of the endogenous constructs
  endoscores <- exogscores%*%t(Gamma_train) %*% solve(diag(nrow(B_train)) - t(B_train))
  
  
  xhat <- endoscores %*% Loadings_train[Cons_endo,,drop = FALSE]
  
  # Denormalize predictions
  xhatrescale= sapply(colnames(xhat),function(x){
    xhat[,x]*sd_train[x]+mean_train[x]
  } )
  
  # Return the predicted values of the indicators connected to endogenous constructs
  Ind_endo=colnames(.object$Information$Model$measurement)[colSums(.object$Information$Model$measurement[Cons_endo,])!=0]

  # Calculate the difference between original and predicted values
  residuals = testData[,Ind_endo] - xhatrescale[,Ind_endo]
  
  list(Pred_Indicators_endo = xhatrescale[,Ind_endo],
       Residuals_Indicators_endo = residuals, 
       Pred_Construct_endo = endoscores)
  
  
  }


 
#' Output for floodlight analysis
#'
#' Calculate the the effect of an independent variable depending on different values
#' of its moderator variable to perform a floodlight analysis \insertCite{Spiller2013;textual}{cSEM}.
#' 
#'
#' @usage effect_moderator_two_way=function(.object=args_default()$.object,
#'                                          .steps = seq(-1.5,1.5,0.01),
#'                                          .alpha = args_default()$.alpha, 
#'                                          .independent = NULL,
#'                                          .moderator = NULL,
#'                                          .dependent = NULL )
#'
#' @inheritParams csem_arguments
#'
#' @references
#'   \insertAllCited{}
#'   
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
effect_moderator_two_way=function(.object=args_default()$.object,
                                  .steps = seq(-1.5,1.5,0.01),
                                  .alpha = args_default()$.alpha, 
                                  .independent = NULL,
                                  .moderator = NULL,
                                  .dependent = NULL ){
  
  # Check whether .object is of class cSEMResults_resampled
  if(!("cSEMResults_resampled" %in% class(.object))){
    stop2('.object needs to be of class cSEMResults_resampled')
  }
  
  # Check whether it is a non-linear model
  if("Nonlinear" != .object$Information$Model$model_type){
    stop2("Model type must be nonlinear.")
  }
  
  # Works only for one significance level, i.e., no vector of significances is allowed
  if(length(.alpha)!=1){
    stop2('Only one significant level is allowed, not a vector')
  }
  
  # Check whether dependent, independent, and moderator variable are provided
  if(is.null(.dependent) | is.null(.independent) | is.null(.moderator)){
    stop2("Not all variables (dependent, independent, and/or moderator) are supplied.")
  }
  
  structural = .object$Information$Model$structural
  
  # Character string containing the names of the dependent and independent variables
  dependent_vars=rownames(structural[rowSums(structural)!=0,,drop=FALSE])
  independent_vars = colnames(structural[,colSums(structural[.dependent,,drop = FALSE])!=0,drop = FALSE])
  
  # Check whether dependent variable is valid
  if(!(.dependent %in% dependent_vars)){
    stop2("Defined dependent variables is not a dependent variable in the specified model.")
  }
  
  
  # Check whether supplied variables are also used in the model
  if(!all(c(.independent, .moderator) %in% colnames(structural))){
    stop2("Independent and moderator variable do not occur together in the equation of the dependent variable.")
  }
  
  # Add stop if the variables are included in a higher-order moderation, e.g., cubic term 
  
  # Possible names of the interaction term. One could adjust the arguement that sth lik x.z needs to be provided to .moderator
  possible_names=c(paste(.independent,.moderator, sep= '.'),paste(.moderator,.independent, sep= '.'))
  
  # Name of the interaction term 
  name_interaction = possible_names[possible_names %in% independent_vars]
  
  if(length(name_interaction) != 1){
    stop2("The defined interaction term does not exist in the model or is not part of the equation of the dependent variable.")
  }
  
  
  # Effect names
  name_mod_effect = paste(.dependent, name_interaction, sep=' ~ ')
  name_single_effect = paste(.dependent, .independent, sep=' ~ ')
  

  dataplot_temp = lapply(.steps, function(x){
    
    # bootstrap effect of x on y at level step, i.e., z=.steps
    effect_boot=.object$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,name_single_effect]+
      .object$Estimates$Estimates_resample$Estimates1$Path_estimates$Resampled[,name_mod_effect]*x 
    
    # Value of the originally estimated effect at level .steps
    effect_at_steps=.object$Estimates$Estimates_resample$Estimates1$Path_estimates$Original[name_single_effect]+
      .object$Estimates$Estimates_resample$Estimates1$Path_estimates$Original[name_mod_effect]*x
    
    out=c(lb=quantile(effect_boot,.alpha/2),
      ub=quantile(effect_boot,1-.alpha/2),
      de=effect_at_steps,
      value=x)
   
  })
  
# Prepare and return output
  
    res=list(out = do.call(rbind,dataplot_temp), 
             Information = c(alpha=.alpha, independent = .independent, moderator= .moderator, dependent = .dependent))
    class(res) = c("Two_Way_Effect", class(res))
    return(res)
} 


# plot method for effect_moderator_two_way
plot.effect_moderator_two_way = function(.TWobject){
  
  require(ggplot2)
  plot1=ggplot(as.data.frame(.TWobject$out),aes(x=.TWobject$out[,4],y=.TWobject$out[,3]))+
    geom_line()+
    geom_ribbon(aes(ymin=.TWobject$out[,1], ymax=.TWobject$out[,2]),alpha=0.2)+
    labs(x=paste('Level of ',.TWobject$Information['moderator']) , 
    y=paste('Effect of', .TWobject$Information['independent'], 'on \n', .TWobject$Information['dependent']))+
    theme_bw()
    # scale_x_continuous(breaks=seq(-3,3,0.5))+
    # theme(panel.grid.minor = element_blank())
  # geom_point(x=-0.224,y=0,size=2)
  plot1
  
  
}
  

