#' Predict indicator scores
#'
#' Predict the indicator scores of endogenous constructs.
#' 
#' Predict uses the procedure introduced by \insertCite{Shmueli2016;textual}{cSEM} in the context of
#' PLS (commonly called: "PLSPredict" \insertCite{Shmueli2019}{cSEM}). 
#' Predict uses k-fold cross-validation to randomly 
#' split the data into training and test data and subsequently predicts the 
#' relevant values in the test data based on the model parameter estimates obtained 
#' using the training data. The number of cross-validation folds is 10 by default but
#' may be changed using the `.cv_folds` argument.
#' By default, the procedure is repeated `.r = 10` times to avoid irregularities
#' due to a particular split. See \insertCite{Shmueli2019;textual}{cSEM} for 
#' details.
#' 
#' Alternatively, users may supply a matrix of `.test_data` with the same column names
#' as those in the data used to obtain `.object` (the training data). 
#' In this case, arguments `.cv_folds` and `.r` are
#' ignored and predict uses the estimated coefficients from `.object` to
#' predict the values in the columns of `.test_data`.
#' 
#' In \insertCite{Shmueli2016;textual}{cSEM} PLS-based predictions for indicator `i`
#' are compared to the predictions based on a multiple regression of indicator `i`
#' on all available exogenous indicators (`.benchmark = "lm"`) and 
#' a simple mean-based prediction summarized in the Q2_predict metric.
#' `predict()` is more general in that is allows users to compare the predictions
#' based on a so-called target model/specificiation to predictions based on an
#' alternative benchmark. Available benchmarks include predictions
#' based on a linear model, PLS-PM weights, unit weights (i.e. sum scores), 
#' GSCA weights, PCA weights, and MAXVAR weights.
#' 
#' Each estimation run is checked for admissibility using [verify()]. If the 
#' estimation yields inadmissible results, `predict()` stops with an error (`"stop"`).
#' Users may choose to `"ignore"` inadmissible results or to simply set predictions
#' to `NA` (`"set_NA"`) for the particular run that failed. 
#'
#' @return An object of class `cSEMPredict` with print and plot methods.
#'   Technically, `cSEMPredict` is a 
#'   named list containing the following list elements:
#'
#' \describe{
#'   \item{`$Actual`}{A matrix of the actual values/indicator scores of the endogenous constructs.}
#'   \item{`$Prediction_target`}{A matrix of the predicted indicator scores of the endogenous constructs 
#'     based on the target model. Target refers to procedure used to estimate 
#'     the parameters in `.object`.}
#'   \item{`$Residuals_target`}{A matrix of the residual indicator scores of the endogenous constructs 
#'     based on the target model.}
#'   \item{`$Residuals_benchmark`}{A matrix of the residual indicator scores 
#'     of the endogenous constructs based on a model estimated by the procedure
#'     given to `.benchmark`.}
#'   \item{`$Prediction_metrics`}{A data frame containing the predictions metrics
#'     MAE, RMSE, and Q2_predict.}
#'   \item{`$Information`}{A list with elements
#'     `Target`, `Benchmark`,
#'     `Number_of_observations_training`, `Number_of_observations_test`, `Number_of_folds`,
#'     `Number_of_repetitions`, and `Handle_inadmissibles`.}
#'     }
#'   
#' @usage predict(
#'  .object               = NULL,
#'  .benchmark            = c("lm", "unit", "PLS-PM", "GSCA", "PCA", "MAXVAR"),
#'  .cv_folds             = 10,
#'  .handle_inadmissibles = c("stop", "ignore", "set_NA"),
#'  .r                    = 10,
#'  .test_data            = NULL
#'  )
#'
#' @inheritParams csem_arguments
#' @param .handle_inadmissibles Character string. How should inadmissible results 
#'   be treated? One of "*stop*", "*ignore*", or "*set_NA*". If "*stop*", [predict()] 
#'   will stop immediatly if estimation yields an inadmissible result.
#'   For "*ignore*" all results are returned even if all or some of the estimates
#'   yielded inadmissible results. 
#'   For "*set_NA*" predictions based on inadmissible parameter estimates are
#'   set to `NA`. Defaults to "*stop*"
#'
#' @seealso [csem], [cSEMResults]
#' 
#' @references
#'   \insertAllCited{}
#'
#' @example inst/examples/example_predict.R
#' 
#' @export

predict <- function(
  .object               = NULL, 
  .benchmark            = c("lm", "unit", "PLS-PM", "GSCA", "PCA", "MAXVAR"),
  .cv_folds             = 10,
  .handle_inadmissibles = c("stop", "ignore", "set_NA"),
  .r                    = 10,
  .test_data            = NULL
  ) {
  
  .benchmark            <- match.arg(.benchmark)
  .handle_inadmissibles <- match.arg(.handle_inadmissibles)
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, predict, 
                  .benchmark = .benchmark, 
                  .cv_folds = .cv_folds,
                  .handle_inadmissibles = .handle_inadmissibles,
                  .r = .r,
                  .test_data = .test_data
                  )
    
    class(out) <- c("cSEMPredict", "cSEMPredict_multi")
    return(out)
  } else {
    
    ## Errors and warnings -------------------------------------------------------
    # Stop if second order
    if(inherits(.object, "cSEMResults_2ndorder")) {
      stop2('Currently, `predict()` is not implemented for models containing higher-order constructs.')
    }
    
    # Stop if second order
    if(all(.object$Information$Model$structural == 0)) {
      stop2("`predict()` requires a structural model.")
    }
    
    # Stop if nonlinear. See Danks et al. (?) for how this can be addressed.
    if(.object$Information$Model$model_type != 'Linear'){
      stop2('Currently, `predict()` works only for linear models.')
    }
    
    # Stop if indicator correlation is not Bravais-Pearson
    if(.object$Information$Type_of_indicator_correlation != 'Pearson'){
      stop2('Currently, `predict()` works only in combination with Pearson correlation.')
    }
    
    ## Get arguments and relevant indicators
    #  Note: It is possible that the original data set used to obtain .object
    #        contains variables that have not been used in the model. These need
    #        to be deleted. Thats why we take the column names of .object$Information$Data.
    args <- .object$Information$Arguments
    indicators <- colnames(.object$Information$Data) # the data used for the estimation 
    # (standardized and clean) with variables
    # ordered accoriding to model$measurement.  
    
    ## Is the benchmark the same as what was used to obtain .object
    if(.benchmark == args$.approach_weights) {
      warning2(
        "The following warning occured in the `predict()` function:\n",
        "Original estimation is based on the same approach as the benchmark approach.",
        " Target and benchmark predicitons are identical."
      )
    }
    
    if(args$.disattenuate & .benchmark %in% c("unit", "GSCA", "MAXVAR") & 
       any(.object$Information$Model$construct_type == "Composite")) {
      args$.disattenuate <- FALSE
      warning2(
        "The following warning occured in the `predict()` function:\n",
        "Disattenuation only applicable if all constructs are modeled as common factors.",
        " Results based on benchmark = `", .benchmark, "` are not disattenuated."
      )
    }
    
    ## Has .test_data been supplied?
    if(!is.null(.test_data)) {
      .r <- 1
      .cv_folds <- NA
      
      dat_train <- args$.data[, indicators]
      
      if(length(setdiff(colnames(.test_data), colnames(dat_train))) > 0) {
        stop2("The following error occured in the `predict()` function:\n",
              "Some variable names in the test data are not part of the training data.")
      }
      # Warn if .test_data doesnt have row names
      if(is.null(rownames(.test_data))) {
        warning2(
          "The following warning occured in the `predict()` function:\n",
          "The test data does not have row names to identify observations. "
        )
      }
      dat_test  <- .test_data[, indicators]
      
      dat <- list("test" = dat_test, "train" = dat_train)
    }
    
    out_all <- list()
    for(i in 1:.r) {
      
      if(is.null(.test_data)) {
        ## Draw cross-validation samples
        dat <- resampleData(
          .object          = .object, 
          .resample_method = "cross-validation", 
          .cv_folds        = .cv_folds,
          .R               = 1,
          .seed            = NULL
        )[[1]]
        
        ## Clean data
        dat <- lapply(dat, processData, .model = .object$Information$Model)
        
        ii <- length(dat)
      } else {
        ii <- 1
      }
      
      out_cv <- list() 
      for(j in 1:ii) {
        
        X_train    <- as.matrix(do.call(rbind, dat[-j]))[, indicators]
        X_test     <- dat[[j]][, indicators]
        
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
        colnames(X_test_scaled) <- colnames(X_test)
        rownames(X_test_scaled) <- rownames(X_test)
        
        ### Estimate target and benchmark using training data and original arguments
        args$.data <- X_train
        
        args_target <- args_benchmark <- args
        
        if(.benchmark %in% c("unit", "PLS-PM", "GSCA", "PCA", "MAXVAR")) {
          args_benchmark$.approach_weights <- .benchmark
          kk <- 2
        } else {
          kk <- 1
        }
        
        ## Run for target and benchmark 
        args_list <- list(args_target, args_benchmark)
        results <- list()
        for(k in 1:kk) {
          
          Est        <- do.call(foreman, args_list[[k]]) 
          
          # Identify exogenous construct in the structural model
          cons_exo  <- Est$Information$Model$cons_exo
          cons_endo <- Est$Information$Model$cons_endo
          
          # Which indicators are connected to endogenous constructs?
          endo_indicators <- colnames(Est$Information$Model$measurement)[colSums(Est$Information$Model$measurement[cons_endo, , drop = FALSE]) != 0]
          # Which indicators are connected to exogenous constructs?
          exo_indicators <- colnames(Est$Information$Model$measurement)[colSums(Est$Information$Model$measurement[cons_exo, , drop = FALSE]) != 0]
          
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
          results[[k]] <- list(X_hat_rescaled, residuals_target)
        }
        
        if(.benchmark %in% c("unit", "PLS-PM", "GSCA", "PCA", "MAXVAR")) {
          predictions_benchmark <- results[[2]][[1]]
          residuals_benchmark   <- results[[2]][[2]]
        } else if(.benchmark == "lm") {
          ## Compute naiv predictions based on a linear model that explains each
          ## endogenous indicator by all exogenous indicators
          beta_exo <- solve(t(X_train[, exo_indicators]) %*% 
                              X_train[, exo_indicators]) %*% 
            t(X_train[, exo_indicators]) %*% X_train[, endo_indicators, drop = FALSE]
          
          predictions_benchmark <- as.matrix(X_test[, exo_indicators]) %*% beta_exo
          residuals_benchmark   <- X_test[, endo_indicators] - predictions_benchmark
        }
        ## Compute naiv mean-based predictions and residuals
        residuals_mb   <- t(t(X_test[, endo_indicators]) - mean_train[endo_indicators])
        
        ## Output
        out_cv[[j]] <- list(
          "Predictions_target"    = results[[1]][[1]],
          "Residuals_target"      = results[[1]][[2]],
          "Predictions_benchmark" = predictions_benchmark,
          "Residuals_benchmark"   = residuals_benchmark,
          "Residuals_mb"          = residuals_mb
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
    
    mae_target    <- apply(out_temp$Residuals_target, 2, function(x) mean(abs(x - mean(x))))
    mae_benchmark <- apply(out_temp$Residuals_benchmark, 2, function(x) mean(abs(x - mean(x))))
    rmse_target   <- apply(out_temp$Residuals_target, 2, function(x) sqrt(mean((x - mean(x))^2)))
    rmse_benchmark<- apply(out_temp$Residuals_benchmark, 2, function(x) sqrt(mean((x - mean(x))^2)))
    
    q2_predict  <- c()
    for(i in colnames(out_temp$Residuals_target)) {
      
      q2_predict[i] <- 1- sum((out_temp$Residuals_target[, i] - mean(out_temp$Residuals_target[, i]))^2) /
        sum((out_temp$Residuals_mb[, i] - mean(out_temp$Residuals_mb[, i]))^2)
    }
    
    ## Create data fram
    df_metrics <- data.frame(
      "Name"           = endo_indicators,
      "MAE_target"     = mae_target,
      "MAE_benchmark" = mae_benchmark,
      "RMSE_target"    = rmse_target,
      "RMSE_benchmark" = rmse_benchmark,
      "Q2_predict"     = q2_predict,
      stringsAsFactors = FALSE
    )
    rownames(df_metrics) <- NULL
    
    out <- list(
      "Actual"      = if(is.null(.test_data)) {
        .object$Information$Arguments$.data[, endo_indicators]
      } else {
        .test_data[, endo_indicators]
      },
      "Predictions_target"  = out_temp$Predictions_target,
      "Residuals_target"    = out_temp$Residuals_target,
      "Residuals_benchmark" = out_temp$Residuals_benchmark,
      "Prediction_metrics"  = df_metrics,
      "Information"         = list(
        "Target"                 = .object$Information$Arguments$.approach_weights,
        "Benchmark"              = .benchmark,
        "Handle_inadmissibles"   = .handle_inadmissibles,
        "Number_of_observations_training" = nrow(X_train),
        "Number_of_observations_test" = nrow(X_test),
        "Number_of_folds"        = .cv_folds,
        "Number_of_repetitions"  = .r
      )
    )
    
    ## Return
    class(out) <- "cSEMPredict"
    out 
  }
}