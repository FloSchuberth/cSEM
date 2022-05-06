#' Predict indicator scores
#'
#'\lifecycle{maturing}
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
#' By default, the procedure is not repeated (`.r = 1`). You may choose to repeat
#' cross-validation by setting a higher `.r` to be sure not to have a particular 
#' (unfortunate) split. See \insertCite{Shmueli2019;textual}{cSEM} for 
#' details. Typically `.r = 1` should be sufficient though.
#' 
#' Alternatively, users may supply a matrix or a data frame of `.test_data` with 
#' the same column names as those in the data used to obtain `.object` (the training data). 
#' In this case, arguments `.cv_folds` and `.r` are
#' ignored and predict uses the estimated coefficients from `.object` to
#' predict the values in the columns of `.test_data`.
#' 
#' In \insertCite{Shmueli2016;textual}{cSEM} PLS-based predictions for indicator `i`
#' are compared to the predictions based on a multiple regression of indicator `i`
#' on all available exogenous indicators (`.benchmark = "lm"`) and 
#' a simple mean-based prediction summarized in the Q2_predict metric.
#' `predict()` is more general in that is allows users to compare the predictions
#' based on a so-called target model/specification to predictions based on an
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
#'     MAE, RMSE, and Q2_predict. In case of categorical indicators, the misclassification error rate
#'     is also included. Please note that the misclassification error rate can only be obtained for categorical
#'     indicators. In case of continuous indicators, the misclassification error rate is set to the
#'     MAE value.}
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
#'  .r                    = 1,
#'  .test_data            = NULL,
#'  .approach_score_target= c("mean", "median", "mode"),
#'  .sim_points           = 100,
#'  .disattenuate         = TRUE,
#'  .treat_as_continuous  = TRUE,
#'  .approach_score_benchmark = c("mean", "median", "mode", "round")
#'  )
#'
#' @inheritParams csem_arguments
#' @param .handle_inadmissibles Character string. How should inadmissible results 
#'   be treated? One of "*stop*", "*ignore*", or "*set_NA*". If "*stop*", [predict()] 
#'   will stop immediately if estimation yields an inadmissible result.
#'   For "*ignore*" all results are returned even if all or some of the estimates
#'   yielded inadmissible results. 
#'   For "*set_NA*" predictions based on inadmissible parameter estimates are
#'   set to `NA`. Defaults to "*stop*"
#' @param .disattenuate Logical. Should the benchmark predictions be based on 
#'   disattenuated parameter estimates? Defaults to `TRUE`.
#'
#' @seealso [csem], [cSEMResults], [exportToExcel()]
#' 
#' @references
#'   \insertAllCited{}
#'
#' @example inst/examples/example_predict.R
#' 
#' @export

predict <- function(
  .object                   = NULL, 
  .benchmark                = c("lm", "unit", "PLS-PM", "GSCA", "PCA", "MAXVAR"),
  .cv_folds                 = 10,
  .handle_inadmissibles     = c("stop", "ignore", "set_NA"),
  .r                        = 1,
  .test_data                = NULL,
  .approach_score_target    = c("mean", "median", "mode"),
  .sim_points               = 100,
  .disattenuate             = TRUE,
  .treat_as_continuous      = TRUE,
  .approach_score_benchmark = c("mean", "median", "mode", "round")
  ) {
  
  .benchmark            <- match.arg(.benchmark)
  .handle_inadmissibles <- match.arg(.handle_inadmissibles)
  .approach_score_target <- match.arg(.approach_score_target)
  .approach_score_benchmark <- match.arg(.approach_score_benchmark)
  
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
    
    if(sum(verify(.object))!=0) {
      stop2('The csem object is not admissible.')
    }
    
    # Stop if second order
    if(all(.object$Information$Model$structural == 0)) {
      stop2("`predict()` requires a structural model.")
    }
    
    # Stop if nonlinear. See Danks et al. (?) for how this can be addressed.
    if(.object$Information$Model$model_type != 'Linear'){
      stop2('Currently, `predict()` works only for linear models.')
    }
    
    #if(!all(.object$Information$Type_of_indicator_correlation == 'Pearson') &&
    #   is.null(.test_data) && .r > 1 && is.null(.test_data)){
    #  stop2('For categorical indicators, only one repetition can be done.')
    #}
    
    
    ## Get arguments and relevant indicators
    #  Note: It is possible that the original data set used to obtain .object
    #        contains variables that have not been used in the model. These need
    #        to be deleted. Thats why we take the column names of .object$Information$Data.
    args <- .object$Information$Arguments
    indicators <- colnames(.object$Information$Data) # the data used for the estimation 
    # (standardized and clean) with variables
    # ordered according to model$measurement.  
    
    ## Is the benchmark the same as what was used to obtain .object
    if(.benchmark == args$.approach_weights && all(.object$Information$Type_of_indicator_correlation == 'Pearson')) {
      warning2(
        "The following warning occured in the `predict()` function:\n",
        "Original estimation is based on the same approach as the benchmark approach.",
        " Target and benchmark predicitons are identical."
      )
    }
    
    if(.disattenuate && .benchmark %in% c("unit", "GSCA", "MAXVAR") && 
       any(.object$Information$Model$construct_type == "Composite")) {
      .disattenuate <- FALSE
      warning2(
        "The following warning occured in the `predict()` function:\n",
        "Disattenuation only applicable if all constructs are modeled as common factors.",
        " Results based on benchmark = `", .benchmark, "` are not disattenuated."
      )
    }
    
    if(.disattenuate & .benchmark %in% c("lm")) {
      .disattenuate <- FALSE
      warning2(
        "The following warning occured in the `predict()` function:\n",
        "Disattenuation is not applicable to benchmark `", .benchmark, "` and ignored."
      )
    }
    
    if(!all(.object$Information$Type_of_indicator_correlation == 'Pearson') && 
       .benchmark != "PLS-PM" && .treat_as_continuous == FALSE) {
      .treat_as_continuous = TRUE
      warning2(
        "The following warning occured in the `predict()` function:\n",
        "The categorical nature of the indicators can currently only be considered",
        "for .benchmark = PLS-PM, the results for benchmark = '", .benchmark, 
         "' treat the indicators as continuous."
      )
    }
    
    ## Has .test_data been supplied?
    if(!is.null(.test_data)) {
      # Is it a matrix or a data.frame?
      if(!any(class(.test_data) %in% c("data.frame", "matrix"))) {
        stop2("The following error occured in the `predict()` function:\n",
              ".test_data must be a matrix or a data frame.")
      }
      .r <- 1
      .cv_folds <- NA
      
      dat_train <- args$.data[, indicators]
      
      
      # Convert to matrix and add rownames
      # Since rownames are required further below to match the observations in the
      # k'th fold of the .r'th run with those of the r+1'th run rownames are also
      # required for the .test_data.
      if(!all(.object$Information$Type_of_indicator_correlation == 'Pearson')){
        .test_data = data.matrix(.test_data)
      }else{
      .test_data = as.matrix(.test_data)
      }
      rownames(.test_data) <- 1:nrow(.test_data)
      # Stop if .test_data does not have column names! As we need to match column
      # names between training and test data.
      if(is.null(colnames(.test_data))) {
        stop2(
          "The following error occured in the `predict()` function:\n",
          "The test data does not have column names that match the training data. "
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
        
        # For categorical indicators, data.matrix has to be used
        # For categorical indicators, both the original matrix as the data matrix
        # are saved.
        if(!all(.object$Information$Type_of_indicator_correlation == 'Pearson')){
          #X_train     <- data.matrix((do.call(rbind, dat[-j])))[, indicators]
          X_train <- do.call(rbind, dat[-j])[, indicators]
          X_test  <- data.matrix(dat[[j]][, indicators])
          #X_test <- data.matrix(dat[[j]][, indicators])
          }else{
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
        
        # If nrow = 1 sapply returns a vector, but We always need a matrix
        if(!inherits(X_test_scaled, "matrix")) {
          X_test_scaled <- matrix(X_test_scaled, nrow = 1)
        }
        
        # Keep rownames to be able to find individual observations
        colnames(X_test_scaled) <- colnames(X_test)
        rownames(X_test_scaled) <- rownames(X_test)
          }
        ### Estimate target and benchmark using training data and original arguments
        args$.data <- X_train
        
        args_target <- args_benchmark <- args
        .disattenuate_benchmark <- .disattenuate
        .disattenuate_target <- args_target$.disattenuate
        
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
          
          # For the target predictions, the parameter estimates are obtained
          # for the target data, while the estimates equal the estimates of 
          # .object if a test data is given.
          if(k == 1){
          if(is.null(.test_data)){
          Est        <- do.call(foreman, args_list[[k]])
          }else{
            Est <- .object
          }
          }else{
            # For the benchmark predictions, the parameter estimates are obtained
            # for the target data. If a test sample is given and .disattenuate
            # equals .disattenuate of .object and if all indicators are numeric
            # and should be treated in their original form, the estimates equal
            # the estimates of .object.
            if((!is.null(.test_data) && args_list[[k]]$.disattenuate == .disattenuate &&
               all(.object$Information$Type_of_indicator_correlation == 'Pearson'))|
               (!is.null(.test_data) && args_list[[k]]$.disattenuate == .disattenuate &&
               .treat_as_continuous == FALSE)){
                 Est <- .object
            }else{
              
            if(args_list[[k]]$.disattenuate != .disattenuate){
              args_list[[k]]$.disattenuate <- .disattenuate
            }
            if(.treat_as_continuous){
              args_list[[k]]$.data <- data.matrix(args_list[[k]]$.data)
            }
              Est        <- do.call(foreman, args_list[[k]])
            }
          }
          
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
            
            # For categorical indicators use prediction from OrdPLS, else 
            # normal prediction is used
            if(!all(.object$Information$Type_of_indicator_correlation == 'Pearson') && k ==1|
               !all(.object$Information$Type_of_indicator_correlation == 'Pearson') &&k == 2 && 
               .treat_as_continuous == FALSE){
              
              if(k == 1){
              # Save the categorical indicators and the continous indicators
              is_numeric_indicator <- lapply(X_train, is.numeric)
              cat_indicators <- names(is_numeric_indicator[is_numeric_indicator == FALSE])
              cont_indicators <- names(is_numeric_indicator[is_numeric_indicator == TRUE])
              }
              
              if(k == 1){
                .approach_score = .approach_score_target
              }else{
                .approach_score = .approach_score_benchmark
              }
              
              if(!all(endo_indicators%in%cat_indicators) && .approach_score == "mode"){
                stop2(
                  "The following error occured in the `predict()` function:\n",
                  "The option '.approach_score = mode' can only be applied to only\n", 
                  "categorical indicators"
                )
              }
              
              if(k == 1){
              X_train <- data.matrix(X_train)
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
              
              colnames(X_test_scaled) <- colnames(X_test)
              rownames(X_test_scaled) <- rownames(X_test)
              
              # Replace the scaled categorical indicators by their original values
              # Reason: Categorical indicators should not be scaled
              X_test_scaled[,cat_indicators] <- X_test[,cat_indicators]
              }
            # get the thresholds for the categorical indicators
            thresholds <- Est$Information$Threshold_parameter_estimates
            thresholds <- thresholds[!is.na(thresholds)]
            Tmin <- -4
            Tmax <- 4
            correction <- function(x) {
              multval <- x[duplicated(x)]
              threshold.indices.to.change <- which(x == multval)
              if (length(threshold.indices.to.change) == 0) return(x)
              x[threshold.indices.to.change] <- multval - epsilon * rev(threshold.indices.to.change - min(threshold.indices.to.change))    
              x
            }
            for (th in 1:length(thresholds)) thresholds[[th]] <- correction(thresholds[[th]])
            
            thresholds <- lapply(thresholds, function(x) c(Tmin, x, Tmax))
            
            Cov_ind <- Est$Estimates$Indicator_VCV
            X_hat <- matrix(0, nrow = nrow(X_test), ncol = length(endo_indicators),
                                     dimnames = list(rownames(X_test),endo_indicators))
            if(all(exo_indicators %in% cont_indicators)){
              # Predict scores for the exogenous constructs
              eta_hat_exo <- X_test_scaled[,exo_indicators]%*%t(W_train[cons_exo, exo_indicators, drop = FALSE])
              
              # Predict scores for the endogenous constructs (structural prediction)
              eta_hat_endo_star <- eta_hat_exo %*% t(Gamma_train) %*% t(solve(diag(nrow(B_train)) - B_train))
              
              # Predict scores for indicators of endogenous constructs (communality prediction)
              X_hat_endo_star <- eta_hat_endo_star %*% loadings_train[cons_endo, endo_indicators, drop = FALSE]
              
              X_hat <- X_hat_endo_star
              
              for(m in colnames(X_hat_endo_star)){
                if(m %in% cat_indicators){
                  for(o in 1:length(X_hat_endo_star[,1])){
                    X_hat[o,m] <- findInterval(X_hat_endo_star[o,m], thresholds[[m]])
                  }
                }
              }
              
              X_hat_rescaled <- sapply(colnames(X_hat), function(x) {
                mean_train[x] + X_hat[, x] * sd_train[x]
              })
              X_hat_rescaled[,cat_indicators] <- X_hat[,cat_indicators]
              
            }else{
              for(o in 1:nrow(X_test)){
                l <- NA
                u <- NA
                for(p in exo_indicators){
                  
                  # Lower and upper thresholds for each observation have to be defined
                  # for the categorical indicators, the estimated thresholds are used
                  # for the continous indicators, -10 and 10 are used
                  if(p %in% cat_indicators){
                    l <- c(l,thresholds[[p]][X_test[o,p]])
                    u <- c(u,thresholds[[p]][X_test[o,p]+1] )
                  }else{
                    l <- c(l, -10)
                    u <- c(u, 10)
                  }
                }
                l <- l[-1]
                u <- u[-1]
                names(l) <- exo_indicators
                names(u) <- exo_indicators
                
                # Simulation of values of the truncated normal distribution for the categorical indicators
                Xstar <- t(TruncatedNormal::mvrandn(l = l, u = u, Cov_ind[exo_indicators, exo_indicators],.sim_points))
                colnames(Xstar) <- exo_indicators
                
                # The continuous indicators are replaced through their original values of the test data
                for(z in colnames(Xstar)){
                  if(z %in% cont_indicators){
                    Xstar[,z] <- X_test_scaled[o,z]
                  }
                }
                
                # Predict scores for the exogenous constructs
                eta_hat_exo <- Xstar[,exo_indicators]%*%t(W_train[cons_exo, exo_indicators, drop = FALSE])
                
                # Predict scores for the endogenous constructs (structural prediction)
                eta_hat_endo_star <- eta_hat_exo %*% t(Gamma_train) %*% t(solve(diag(nrow(B_train)) - B_train))
                
                # Predict scores for indicators of endogenous constructs (communality prediction)
                X_hat_endo_star <- eta_hat_endo_star %*% loadings_train[cons_endo, endo_indicators, drop = FALSE]
                
                # Aggregation of the npred estimations via mean or median
                if(.approach_score == "mean"){
                  X_hat[o,] <- apply(X_hat_endo_star,2,mean)
                }else if(.approach_score == "median"){
                  X_hat[o,] <- apply(X_hat_endo_star,2,median)
                }else if(.approach_score == "mode"){
                  for(z in colnames(X_hat_endo_star)){
                    breaks <- thresholds[[z]]
                    dupl <- which(duplicated(thresholds[[z]]))
                    if (length(dupl) > 0) breaks <- thresholds[[z]][-dupl]
                    if (min(X_hat_endo_star[, z]) < breaks[1]) breaks[1] <- min(X_hat_endo_star[, z])
                    if (max(X_hat_endo_star[, z]) > breaks[length(breaks)]) breaks[length(breaks)] <- max(X_hat_endo_star[, z])
                    
                    X_hat[o,z] <- which.max(graphics::hist(X_hat_endo_star[, z], breaks = breaks, plot = FALSE)$density)
                  }
                }
                
                if(.approach_score == "mean" || .approach_score == "median"){
                  # Calculating the categorical variables for categorical indicators, 
                  # for continuous indicators, the mean/median values are saved
                  for(m in colnames(X_hat)){
                    if(m %in% cat_indicators){
                      X_hat[o,m] <- findInterval(X_hat[o,m], thresholds[[m]])
                    }
                  }
                }
                
                
              }
              # Rescale the continuous indicators
                X_hat_rescaled <- sapply(colnames(X_hat), function(x) {
                  mean_train[x] + X_hat[, x] * sd_train[x]
                })
                X_hat_rescaled[,endo_indicators[which(endo_indicators %in% cat_indicators)]] <- X_hat[,endo_indicators[which(endo_indicators %in% cat_indicators)]]
            }
            
            }else{
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

            # If nrow = 1 sapply returns a vector, but We always need a matrix
            if(!inherits(X_hat_rescaled, "matrix")) {
              X_hat_rescaled <- matrix(X_hat_rescaled, nrow = 1)
              colnames(X_hat_rescaled) <- colnames(X_hat)
              rownames(X_hat_rescaled) <- rownames(X_hat)
            }
            
            # Select only endogenous indicators
            X_hat_rescaled <- X_hat_rescaled[, endo_indicators, drop = FALSE]
            }
            # Calculate the difference between original and predicted values
            residuals_target <- X_test[, endo_indicators, drop = FALSE] - 
                                X_hat_rescaled[, endo_indicators, drop = FALSE]
            
          } else if(status_code != 0 & .handle_inadmissibles == "set_NA"){
            X_hat_rescaled  <- residuals_target <- X_test[, endo_indicators, drop = FALSE] 
            X_hat_rescaled[] <- NA
            residuals_target[] <- NA
          } else {
            stop2("Estimation based on one of the cross-validation folds yielded an inadmissible results.\n",
                  " Consider setting handle_inadmissibles = 'ignore'.")
          }
          results[[k]] <- list(X_hat_rescaled, residuals_target)
        }
        
        if(.benchmark %in% c("unit", "PLS-PM", "GSCA", "PCA", "MAXVAR")) {
          if(!all(.object$Information$Type_of_indicator_correlation == 'Pearson')){
          predictions_benchmark <- results[[2]][[1]]
          predictions_benchmark[,colnames(predictions_benchmark) %in% cat_indicators] <- round(predictions_benchmark[,colnames(predictions_benchmark) %in% cat_indicators])
          residuals_benchmark   <- results[[2]][[2]]
          residuals_benchmark[,colnames(residuals_benchmark) %in% cat_indicators] <- round(residuals_benchmark[,colnames(residuals_benchmark) %in% cat_indicators])
          
          }else{
            predictions_benchmark <- results[[2]][[1]]
            residuals_benchmark   <- results[[2]][[2]]
          }
        } else if(.benchmark == "lm") {
          ## Compute naive predictions based on a linear model that explains each
          ## endogenous indicator by all exogenous indicators
          beta_exo <- solve(t(X_train[, exo_indicators]) %*% 
                              X_train[, exo_indicators]) %*% 
            t(X_train[, exo_indicators]) %*% X_train[, endo_indicators, drop = FALSE]
          if(!all(.object$Information$Type_of_indicator_correlation == 'Pearson')){
          predictions_benchmark <- as.matrix(X_test[, exo_indicators]) %*% beta_exo
          predictions_benchmark[,colnames(predictions_benchmark) %in% cat_indicators] <- round(predictions_benchmark[,colnames(predictions_benchmark) %in% cat_indicators])
          
          residuals_benchmark   <- X_test[, endo_indicators] - predictions_benchmark
          residuals_benchmark[,colnames(residuals_benchmark) %in% cat_indicators] <- round(residuals_benchmark[,colnames(residuals_benchmark) %in% cat_indicators])
          }else{
            predictions_benchmark <- as.matrix(X_test[, exo_indicators]) %*% beta_exo
            residuals_benchmark   <- X_test[, endo_indicators, drop = FALSE] - predictions_benchmark 
          }
        }
        ## Compute naive mean-based predictions and residuals
        residuals_mb   <- t(t(X_test[, endo_indicators, drop = FALSE]) - mean_train[endo_indicators])
      
        
        ## Output
        out_cv[[j]] <- list(
          "Predictions_target"    = results[[1]][[1]],
          "Residuals_target"      = results[[1]][[2]],
          "Predictions_benchmark" = predictions_benchmark,
          "Residuals_benchmark"   = residuals_benchmark,
          "Residuals_mb"          = residuals_mb,
          "Actual"                = X_test[,endo_indicators, drop = FALSE]
        )
      } # END for j in 1:length(dat)  
      
      out_temp <- lapply(purrr::transpose(out_cv), function(x) {
        x <- do.call(rbind, x)
        x <- x[order(as.numeric(rownames(x))), , drop = FALSE]
        x
      })
      
      out_all[[i]] <- out_temp
    }
    
    #out_temp<- lapply(purrr::transpose(out_all), function(x) {
    #  x <- do.call(rbind, x)
    #  x <- x[order(as.numeric(rownames(x))), ]
    #  x
    #})
    
    out_temp <- purrr::transpose(out_all)
    
    # Compute average prediction over all .r runs that are not NA
    #out_temp <- lapply(purrr::transpose(out_all), function(x) {
      
    #  a <- apply(abind::abind(x, along = 3), 1:2, function(y) sum(y, na.rm = TRUE))
    #  b <- Reduce("+", lapply(x, function(y) !is.na(y)))
      
    #  a / b
    #})
    df_metrics <- list()
    for(q in 1: length(out_all)){
    ## Compute prediction metrics ------------------------------------------------
    Res_t <- out_all[[q]]$Residuals_target
    Res_b <- out_all[[q]]$Residuals_benchmark
    Pred_t <- out_all[[q]]$Predictions_target
    Pred_b <- out_all[[q]]$Predictions_benchmark
    act   <- out_all[[q]]$Actual

    mse2_target = calculateMSE2(pred = Pred_t, act = act, resid = Res_t)
    mse2_benchmark = calculateMSE2(pred = Pred_b, act = act, resid = Res_b)
    
    ## Create data frame
    df_metrics[[q]] <- data.frame(
      "MAE_target"     = calculateMAE(resid = Res_t),
      "MAE_benchmark"  = calculateMAE(resid = Res_b),
      "RMSE_target"    = calculateRMSE(resid = Res_t),
      "RMSE_benchmark" = calculateRMSE(resid = Res_b),
      "Q2_predict"     = calculateq2(res = Res_t, MB = out_all[[q]]$Residuals_mb),
      "misclassification_target"    = calculateMissclassification(resid = Res_t),
      "misclassification_benchmark" = calculateMissclassification(resid = Res_b),
      "MAPE_target"    = calculateMAPE(resid = Res_t, act = act),
      "MAPE_benchmark" = calculateMAPE(resid = Res_b, act = act),
      "MSE2_target"    = mse2_target,
      "MSE2_benchmark" = mse2_benchmark,
      "U1_target"      = calculateU1(act = act, mse2 = mse2_target),
      "U1_benchmark"   = calculateU1(act = act, mse2 = mse2_benchmark),
      "U2_target"      = calculateU2(act = act, resid = Res_t),
      "U2_benchmark"   = calculateU2(act = act, resid = Res_b),
      "UM_target"      = calculateUM(act = act, pred = Pred_t, mse2 = mse2_target),
      "UM_benchmark"   = calculateUM(act = act, pred = Pred_b, mse2 = mse2_benchmark),
      "UR_target"      = calculateUR(pred = Pred_t, act = act, mse2 = mse2_target),
      "UR_benchmark"   = calculateUR(pred = Pred_b, act = act, mse2 = mse2_benchmark),
      "UD_target"      = calculateUD(pred = Pred_t, act = act, mse2 = mse2_target),
      "UD_benchmark"   = calculateUD(pred = Pred_b, act = act, mse2 = mse2_benchmark),
      stringsAsFactors = FALSE
    )
    #}
    }
    df_metrics <- data.frame(
      "Name" = endo_indicators,
      Reduce("+", df_metrics)/length(df_metrics))
    rownames(df_metrics) <- NULL
    
    if(.benchmark == "PLS-PM" && !all(.object$Information$Type_of_indicator_correlation == 'Pearson')
       && .treat_as_continuous == FALSE){
      if(.approach_score_benchmark != "round"){
      .benchmark = "OrdPLS"
      }else{
        .benchmark = "OrdPLS rounded"
      }
    }
    if(.object$Information$Arguments$.approach_weights == "PLS-PM" && !all(.object$Information$Type_of_indicator_correlation == 'Pearson')){
      .target <- "OrdPLS"
    }else{
      .target = .object$Information$Arguments$.approach_weights
    }
    
    out <- list(
     # "Actual"      = if(is.null(.test_data)) {
      #  .object$Information$Arguments$.data[, endo_indicators]
      #} else {
      #  .test_data[, endo_indicators]
      #},
      "Actual"              = out_temp$Actual,
      "Predictions_target"  = out_temp$Predictions_target,
      "Residuals_target"    = out_temp$Residuals_target,
      "Residuals_benchmark" = out_temp$Residuals_benchmark,
      "Prediction_metrics"  = df_metrics,
      "Information"         = list(
        "Estimator_target"          = .target,
        "Estimator_benchmark"       = .benchmark,
        "Disattenuation_target"     = .disattenuate_target,
        "Disattenuation_benchmark"  = .disattenuate_benchmark,
        "Handle_inadmissibles"      = .handle_inadmissibles,
        "Number_of_observations_training" = nrow(X_train),
        "Number_of_observations_test" = nrow(X_test),
        "Number_of_folds"           = .cv_folds,
        "Number_of_repetitions"     = .r
      )
    )
    
    ## Return
    class(out) <- "cSEMPredict"
    out 
  }
}
