### Use this file to add new functions whose name you have not yet decided on
# or when it is unclear where it belongs.

# Empirical Bayes correction for loadings suggested by Dijkstra (2018) can be 
# extended to other parameter as well .
# I recommend to extend it to construct correlations instead of the path 
# coefficients and restimate them as we know the bounds of such correlations
# However, this requires that we bootstrap the SEs for the construct correlations.

## Things to discuss
# 1. Correction only for one single data set not each bootstrap data set, as this
#    does not make much sense because we can always set .handle_inadmissibles = "replace"
# 2. Not clear how to tackle inadmissible results due to not semi-positive 
#    construct correlation matrix. not semi-positve definite does not mean
#    that correlations are larger/smaller than 1 /-1...
# 2. Not clear how to handle reliabilities. 
#    - if loadings are corrected: reliabilities change --> construct correlations
#      change --> path coefficients change
#    - if construct correlations are corrected, path coefficients change but for 
#      nonlinear models we would need the reliabilities 
correctInadmissibles <- function(
  .object = NULL, 
  .method = c("median","mean")
  ){
  
  .object <- res1
  
  if(inherits(.object, "cSEMResults_2ndorder")) {

    stop("Not implemented yet")
    
  } else {
    x2 <- .object$Estimates
    sd <- cSEM:::SdResample(.object$Estimates$Estimates_resample$Estimates1,
                            .resample_method = "bootstrap") 
  }

  if(verify(.object)["2"]) {
    ### Loadings larger than 1 in absolute value
    L    <- x2$Loading_estimates[.object$Information$Model$measurement != 0]
    L_sd <- sd$Loading_estimates
    
    for(i in which(abs(L) > 1)) {
      A <- pnorm((-1 - L[i]) / L_sd[i])
      B <- pnorm((1 - L[i]) / L_sd[i]) 
      
      ## Replace VCV elements
      if(.method == "median"){
        ## Overwrite the old loadings
        L[i] <- L[i] + L_sd[i] * qnorm(mean(c(A , B)))
      }
      
      if(.method == "mean"){
        ## Overwrite the old loadings
        L[i] <- L[i] + L_sd[i] / (B - A) * 
          (dnorm((-1 - L[i]) / L_sd[i]) - 
             dnorm((1 - L[i]) / L_sd[i]))
      }
    }
    # Overwrite original construct VCV
    .object$Estimates$Loading_estimates[.object$Information$Model$measurement != 0] <- L
    
    # Recompute reliabilities based on the new loadings
    L  <- .object$Estimates$Loading_estimates
    W  <- .object$Estimates$Weight_estimates
    Q2 <- .object$Estimates$Reliabilities
    
    for(j in rownames(L)) {
      Q2[j] <- c(W[j, ] %*% L[j, ])^2
    }
    
    .object$Estimates$Reliabilities <- Q2
  }
  
  if(verify(.object)["3"]) {
    ### Construct VCV not positive semi-definite (because construct correlations 
    ### are larger than 1)
    VCV <- x2$Construct_VCV
    VCV_sd <- matrix(sd$User_fun, nrow = nrow(VCV), ncol = ncol(VCV),
                     dimnames = list(rownames(VCV), colnames(VCV)))
    
    for(j in rownames(VCV)) {
      for(i in colnames(VCV)) {
        if(abs(VCV[j, i]) > 1) {
          A <- pnorm((-1 - VCV[j, i]) / VCV_sd[j, i])
          B <- pnorm((1 - VCV[j, i]) / VCV_sd[j, i]) 
          
          ## Replace VCV elements
          if(.method == "median"){
            ## Overwrite the old loadings
            VCV[j, i] <- VCV[j, i] + VCV_sd[j, i] * qnorm(mean(c(A , B)))
          }
          
          if(.method == "mean"){
            ## Overwrite the old loadings
            VCV[j, i] <- VCV[j, i] + VCV_sd[j, i] / (B - A) * 
              (dnorm((-1 - VCV[j, i]) / VCV_sd[j, i]) - 
                 dnorm((1 - VCV[j, i]) / VCV_sd[j, i]))
          }
        }
      }
    }
    # Overwrite original construct VCV
    .object$Estimates$Construct_VCV <- VCV
    
    # Reestimate path model
    P <- .object$Estimates$Construct_VCV
    
    path <- cSEM:::estimatePath(
      .approach_nl    = .object$Information$Arguments$.approach_nl,
      .approach_paths = .object$Information$Arguments$.approach_paths,
      .csem_model     = .object$Information$Arguments$.model,
      .H              = .object$Estimates$Construct_scores,
      .normality      = .object$Information$Arguments$.normality,
      .P              = .object$Estimates$Construct_VCV,
      .Q              = .object$Estimates$Reliabilities # only for nonlinear models
    )
    
    .object$Estimates$Path_estimates <- path
  }
  
  ## Return
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