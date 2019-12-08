#' Predict indicator scores
#'
#' Predict the indicator scores of a endogenous constructs.
#' 
#' Details and argument description
#'
#' @return An object of class cSEMPredict.
#'   
#' @usage predict(
#'  .object               = NULL,
#'  .cv_folds             = 10,
#'  .handle_inadmissibles = c("stop", "ignore", "set_NA"),
#'  .r                    = 10,
#'  .seed                 = NULL
#'  )
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @example inst/examples/example_predict.R
#' 
#' @export

predict <- function(
  .object               = NULL, 
  .cv_folds             = 10,
  .handle_inadmissibles = c("stop", "ignore", "set_NA"),
  .r                    = 10, 
  .seed                 = NULL
  ) {
  
  .handle_inadmissibles <- match.arg(.handle_inadmissibles)
  ## Errors and warnings -------------------------------------------------------
  # Stop if nonlinear. See Danks et al. (?) for how this can be addressed.
  if(.object$Information$Model$model_type != 'Linear'){
    stop2('Currently, `predict()` works only for linear models.')
  }
  
  # Stop if indicator correlation is not Bravais-Pearson
  if(.object$Information$Type_of_indicator_correlation != 'Pearson'){
    stop2('Currently, `predict()` works only in combination with Pearson correlation.')
  }
  
  # Stop if indicator correlation is not Bravais-Pearson
  if(inherits(.object, "cSEMResults_2ndorder")) {
    stop2('Currently, `predict()` is not implemented for models containing higher-order constructs.')
  }
  
  ## Get relevant quantities
  args <- .object$Information$Arguments
  
  out_all <- list()
  for(i in 1:.r) {
    
    ## Draw cross-validation samples
    dat <- resampleData(
      .object          = .object, 
      .resample_method = "cross-validation", 
      .cv_folds        = .cv_folds,
      .R               = 1,
      .seed            = NULL
    )[[1]]
    
    out_cv <- list() 
    for(j in 1:length(dat)) {
      
      X_train    <- as.matrix(do.call(rbind, dat[-j]))
      X_test     <- dat[[j]]
      
      mean_train      <- colMeans(X_train)
      sd_train        <- matrixStats::colSds(as.matrix(X_train))
      names(sd_train) <- names(mean_train)
      
      # Scale the test data set with the descriptives of the training data set
      # Reason: estimated parameters are based on a standardized data set, i.e., in
      #         a standardized metric. Observations of the test data must have the
      #         same scale as these estimates. 
      X_test_scaled <- sapply(1:ncol(X_test), function(x){
        (X_test[, x] - mean_train[x]) / sd_train[x]
      })
      
      # Keep rownames to be able to find individual observations
      rownames(X_test_scaled) <- rownames(X_test)
      
      # Estimate using dat and original arguments
      args$.data <- X_train
      Est        <- do.call(foreman, args) 
      
      # Identify exogenous construct in the structural model
      cons_exo  <- Est$Information$Model$cons_exo
      cons_endo <- Est$Information$Model$cons_endo
      
      # Which indicators are connected to endogenous constructs?
      endo_indicators <- colnames(Est$Information$Model$measurement)[colSums(Est$Information$Model$measurement[cons_endo, ]) != 0]
      # Which indicators are connected to exogenous constructs?
      exo_indicators <- setdiff(colnames(Est$Information$Model$measurement), endo_indicators)
      
      W_train        <- Est$Estimates$Weight_estimates
      loadings_train <- Est$Estimates$Loading_estimates
      path_train     <- Est$Estimates$Path_estimates
      
      # Path coefficients of exogenous and endogenous constructs
      B_train      <- path_train[cons_endo, cons_endo, drop = FALSE]
      Gamma_train  <- path_train[cons_endo, cons_exo, drop = FALSE]
      
      # Check status
      status_code <- sum(unlist(verify(Est)))
      
      ## Compute predictions based on path and measurement model ("target prediction")
      # Compute predictions if status is ok or inadmissibles should be ignored
      if(status_code == 0 | (status_code != 0 & .handle_inadmissibles == "ignore")) {
        
        ## Predict scores for the exogenous constructs (validity prediction)
        eta_hat_exo  <- X_test_scaled %*% t(W_train[cons_exo, ,drop = FALSE])
        
        # Predict scores for the endogenous constructs (structural prediction)
        eta_hat_endo <- eta_hat_exo %*% t(Gamma_train) %*% t(solve(diag(nrow(B_train)) - B_train))
        
        # Predict scores for indicators of endogenous constructs (communality prediction)
        X_hat <- eta_hat_endo %*% loadings_train[cons_endo, , drop = FALSE]
        
        # Denormalize predictions
        X_hat_rescaled <- sapply(colnames(X_hat), function(x) {
          mean_train[x] + X_hat[, x] * sd_train[x]
        })
        
        # Select only endogenous indicators
        X_hat_rescaled <- X_hat_rescaled[, endo_indicators]
        
        # Calculate the difference between original and predicted values
        residuals_target <- X_test[, endo_indicators] - X_hat_rescaled[, endo_indicators]
      } else if(status_code != 0 & .handle_inadmissibles == "set_NA"){
        X_hat_rescaled  <- residuals_target <- X_test[, endo_indicators] 
        X_hat_rescaled[] <- NA
        residuals_target[] <- NA
      } else {
        stop2("Estimation based on one of the cross-validation folds yielded an inadmissible results.\n",
              " Consider setting handle_inadmissibles = 'ignore'.")
      }
      ## Compute naiv mean-based predictions
      residuals_mb <- t(t(X_test[, endo_indicators]) - mean_train[endo_indicators])
      
      ## Compute naiv predictions based on a linear model that explains each
      ## endogenous indicator by all exogenous indicators
      beta_exo <- solve(t(X_train[, exo_indicators]) %*% 
                          X_train[, exo_indicators]) %*% 
        t(X_train[, exo_indicators]) %*% X_train[, endo_indicators, drop = FALSE]
      
      X_hat_lm <- as.matrix(X_test[, exo_indicators]) %*% beta_exo
      residuals_lm <- X_test[, endo_indicators] - X_hat_lm
      
      ## Output
      out_cv[[j]] <- list(
        "Predictions_target" = X_hat_rescaled,
        "Residuals_target"   = residuals_target,
        "Residuals_mb"       = residuals_mb,
        "Residuals_lm"       = residuals_lm
      )
    } # END for j in 1:length(dat)  
    
    out_temp <- lapply(purrr::transpose(out_cv), function(x) {
      x <- do.call(rbind, x)
      x <- x[order(as.numeric(rownames(x))), ]
      x
    })
    
    out_all[[i]] <- out_temp
  }
  
  # Compute average prediction over all .r runs that are not NA
  out_temp <- lapply(purrr::transpose(out_all), function(x) {
    
    a <- apply(abind::abind(x, along = 3), 1:2, function(y) sum(y, na.rm = TRUE))
    b <- Reduce("+", lapply(x, function(y) !is.na(y)))

    a / b
  })
  
  ## Compute prediction metrics ------------------------------------------------
  
  mae_target  <- apply(out_temp$Residuals_target, 2, function(x) mean(abs(x - mean(x))))
  mae_lm      <- apply(out_temp$Residuals_lm, 2, function(x) mean(abs(x - mean(x))))
  rmse_target <- apply(out_temp$Residuals_target, 2, function(x) sqrt(mean((x - mean(x))^2)))
  rmse_lm     <- apply(out_temp$Residuals_lm, 2, function(x) sqrt(mean((x - mean(x))^2)))
  
  q2_predict  <- c()
  for(i in colnames(out_temp$Residuals_target)) {
    
    q2_predict[i] <- 1- sum((out_temp$Residuals_target[, i] - mean(out_temp$Residuals_target[, i]))^2) /
      sum((out_temp$Residuals_lm[, i] - mean(out_temp$Residuals_lm[, i]))^2)
  } 
  
  ## Create data fram
  df_metrics <- data.frame(
    "Name"        = endo_indicators,
    "MAE_target"  = mae_target,
    "RMSE_target" = rmse_target,
    "MAE_lm"      = mae_lm,
    "RMSE_lm"     = rmse_lm,
    "Q2_predict"  = q2_predict,
    stringsAsFactors = FALSE
  )
  rownames(df_metrics) <- NULL
  
  out <- list(
    "Actual_target"      = .object$Information$Data[, endo_indicators],
    "Predictions_target" = out_temp$Predictions_target,
    "Residuals_target"   = out_temp$Residuals_target,
    "Residuals_lm"       = out_temp$Residuals_target,
    "Prediction_metrics" = df_metrics,
    "Information"        = list(
      "Number_of_observations" = nrow(.object$Information$Data),
      "Number_of_observations_training" = nrow(X_train),
      "Number_of_observations_test" = nrow(X_test),
      "Number_of_folds"        = .cv_folds,
      "Number_of_repetitions"  = .r,
      "Handle_inadmissibles"   = .handle_inadmissibles
    )
  )
  
  ## Return
  class(out) <- "cSEMPredict"
  out
}