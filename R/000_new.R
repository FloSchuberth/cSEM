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


 
#' Floodlight analysis
#'
#' Calculate the the effect of an independent variable (z) on a dependent variable
#' (y) conditional on the values of a (continous) moderator variable (x) 
#' to perform a floodlight analysis \insertCite{Spiller2013}{cSEM}. Moreover, 
#' the Johnson-Neyman points are calculated, i.e. the value(s) of x for which 
#' lower or upper boundary of the confidence interval of the effect
#' estimate of z for a given x switch signs. 
#' 
#' @usage doFloodlightAnalysis(
#'  .object        = NULL,
#'  .alpha         = args_default()$.alpha,
#'  .y             = args_default()$.y, 
#'  .x             = args_default()$.x,
#'  .z             = args_default()$.z,
#'  .n_spotlights  = args_default()$.n_spotlights
#'  )
#'
#' @inheritParams csem_arguments
#'
#' @references
#'   \insertAllCited{}
#'   
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
#' 
doFloodlightAnalysis <- function(
  .object        = NULL,
  .alpha         = args_default()$.alpha,
  .y             = args_default()$.y, 
  .x             = args_default()$.x,
  .z             = args_default()$.z,
  .n_spotlights  = args_default()$.n_spotlights
){
  
  ### Errors and warnings ------------------------------------------------------
  ## Check whether .object is of class cSEMResults_resampled; if not perform
  ## standard resampling (bootstrap with .R = 499 reps)
  if(!inherits(.object, "cSEMResults_resampled")) {
    if(inherits(.object, "cSEMResults_default")) {
      args <- .object$Information$Arguments
    } else {
      args <- .object$Second_stage$Information$Arguments_original
    }
    args[".resample_method"] <- "bootstrap"
    .object <- do.call(csem, args)
  }
  
  ##  Select relevant quantities
  if(inherits(.object, "cSEMResults_default")) {
    m   <- .object$Information$Model
    est <- .object$Estimates$Estimates_resample$Estimates1$Path_estimates
    H   <- .object$Estimates$Construct_scores
  } else {
    m   <- .object$Second_stage$Information$Model
    est <- .object$Second_stage$Information$Resamples$Estimates$Estimates1$Path_estimates
    H   <- .object$Second_stage$Estimates$Construct_scores
  }
  
  # Character string containing the names of the dependent and independent variables
  dep_vars   <- rownames(m$structural[rowSums(m$structural) !=0, , drop = FALSE])
  indep_vars <- colnames(m$structural[, colSums(m$structural[.y, ,drop = FALSE]) !=0 , drop = FALSE])
  
  ## Check if model is non-linear.
  if(m$model_type != "Nonlinear"){
    stop2(
      "The following error occured in the `doFloodlightAnalysis()`` function:\n",
      "The structural model must be non-linear (i.e. contain interactions).")
  }
  
  ## Works only for one significance level, i.e., no vector of significances is allowed
  if(length(.alpha) != 1){
    stop2(
      "The following error occured in the `doFloodlightAnalysis()`` function:\n",
      "Currently only a single significance level (.alpha) is allowed.")
  }
  
  ## Check whether dependent, independent, and moderator variable are provided
  if(is.null(.y) | is.null(.x) | is.null(.z)){
    stop2(
      "The following error occured in the `doFloodlightAnalysis()`` function:\n",
      "All variables (.y, .x, and .z) must be supplied.")
  }
  
  # Check if the name of the dependent variable is valid
  if(!(.y %in% dep_vars)){
    stop2(
      "The following error occured in the `doFloodlightAnalysis()`` function:\n",
      "The dependent variable supplied to `.y` is not a dependent variable in the original model.")
  }
  
  # Check if the name of the moderator (.x) and the dependent variable (.z) supplied are used 
  # in the original model
  if(!all(c(.x, .z) %in% colnames(m$structural))){
    stop2(
      "The following error occured in the `doFloodlightAnalysis()`` function:\n",
      "Independent and/or moderator variable are not part of the original model.")
  }
  
  ### Calculation --------------------------------------------------------------
  # Possible names of the interaction terms
  possible_names <- c(paste(.x, .z, sep = '.'), paste(.z, .x, sep= '.'))
  
  # Name of the interaction term
  xz <- possible_names[possible_names %in% indep_vars]
  
  if(length(xz) != 1){
    stop2(
      "The defined interaction term does not exist in the model or is not ",
      " part of the equation of the dependent variable.")
  }
  
  
  # Effect names
  beta_z  <- paste(.y, .x, sep = ' ~ ')
  beta_xz <- paste(.y, .z, sep = ' ~ ')
  
  ## Compute spotlights (= effects of independent (z) on dependent (y) for given 
  ## values of moderator (x))
  steps <- seq(min(H[, .x]), max(H[, .x]), length.out = .n_spotlights)
  
  dataplot_temp <- lapply(steps, function(step){
    ## Note:
    #   
    #   y = a + bz + cx + dxz
    #   The marginal effect delta y = b + dx evaluated at x0 is identical to the
    #   b' of the regression that uses x' = x - x0, since
    #   y = a' + b'z + c'x' + d'x'z
    #     = (a - cx0) + (b - dx0)z + cx + dxz
    #   ---> b' = b - dx0
    #
    # Resamples of the effect of z on y at different levels of x 
    effect_resampled <- est$Resampled[ , beta_z] + est$Resampled[, beta_xz] * step 
    
    # Value of the originally estimated effect of z on y at different levels of x
    effect_original <- est$Original[beta_z] + est$Original[beta_xz] * step
    
    # Compute empirical quantile based on resamples
    bounds <- quantile(effect_resampled , c(.alpha/2, 1 - .alpha/2))
    
    # Return output
    c(effect_original, step, bounds[1],  bounds[2])
  })
  
  out <- do.call(rbind, dataplot_temp)
  colnames(out) <- c('direct_effect', 'value_z', 'lb', 'ub')
  
  # Determine Johnson-Neyman point 
  # Look for sign flips in the upper boundary
  pos_ub <- which(diff(sign(out[, 'ub'])) != 0)
  pos_lb <- which(diff(sign(out[, 'lb'])) != 0) 
  
  # Prepare and return output
  out <- list(
    "out"                   = out, 
    "Johnson_Neyman_points" = list(
      JNlb = c(x = out[,'value_z'][pos_lb], 
               y = out[,'lb'][pos_lb]),
      JNub = c(x = out[,'value_z'][pos_ub], 
               y = out[,'ub'][pos_ub]
      )
    ),
    "Information" = list(
      alpha       = .alpha, 
      dependent   = .y,
      independent = .z,
      moderator   = .x
      )
    )
  
  class(out) = "cSEMFloodlight"
  return(out)
} 

#' `cSEMFloodlight` method for `plot()`
#'
#' Plot the direct effect of an independent variable (z) on a dependent variable
#' (y) conditional on the values of a moderator variable (x), including
#' the confidence interval and the Johnson-Neyman points. 
#'
#' @usage plot.cSEMFloodlight(x, ...)
#' 
#' @param x An R object of class [cSEMResults] resulting from a call to [csem()].
#' @param ... ignored.
#' 
#' @export
plot.cSEMFloodlight <- function(x, ...) {
  
  ## Install ggplot2 if not already installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop2(
      "Package `ggplot2` required. Use `install.packages(\"ggplot2\")` and rerun.")
  }
  
  plot1 <- ggplot2::ggplot(as.data.frame(x$out), ggplot2::aes(x = x$out[, 'value_z'], 
                                                     y = x$out[, 'direct_effect'])) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = x$out[,'lb'], ymax = x$out[, 'ub']), alpha = 0.2) +
    ggplot2::labs(
      "x" = paste('Level of ', x$Information['moderator']), 
      "y" = paste('Effect of', x$Information['independent'], 'on \n', x$Information['dependent'])) +
    ggplot2::theme_bw() +
    # scale_x_continuous(breaks=seq(-3,3,0.5))+
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  
  # Add Johnson-Neyman points, if they exist in the considered range
  if(length(x$Johnson_Neyman_points$JNlb) == 2){
    JN <- x$Johnson_Neyman_points$JNlb
    plot1 <- plot1 +
      ggplot2::geom_point(x = JN['x'], y = JN['y'], size = 2)  
  }
  
  if(length(x$Johnson_Neyman_points$JNub) == 2){
    JN = x$Johnson_Neyman_points$JNub
    plot1 = plot1 +
      ggplot2::geom_point(x = JN['x'], y = JN['y'], size = 2)  
  }
  
  # Plot
  plot1
}
