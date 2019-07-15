### Use this file to add new functions whose name you have not yet decided on
# or when it is unclear where it belongs.


 
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