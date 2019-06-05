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


 
#' Output of floodlight analysis
#'
#' Calculate the the effect of an independent variable depending on different values
#' of its moderator variable to perform a floodlight analysis \insertCite{Spiller2013}{cSEM}.
#' Moreover, the Johnson-Neyman point(s) is calculated. In doing so, it is considered 
#' whether there is a sign switch in the lower and upper confidence intervals. 
#' If so, the position of the switch is returned.
#' 
#'
#' @usage effect_moderator_two_way=function(.object=args_default()$.object,
#'                                          .steps = seq(-1.5,1.5,0.01),
#'                                          .alpha = args_default()$.alpha, 
#'                                          .dependent = NULL,
#'                                          .independent = NULL,
#'                                          .moderator = NULL )
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
                                  .steps = args_default()$.steps,
                                  .alpha = args_default()$.alpha,
                                  .dependent = args_default()$.dependent, 
                                  .independent = args_default()$.independent,
                                  .moderator = args_default()$.moderator ){
  
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
    
    bounds = quantile(effect_boot,c(.alpha/2,1-.alpha/2))

    # Prepare and return output
    out = c(bounds[1],
      bounds[2],
      de = effect_at_steps,
      value = x)
    names(out) = c('lb', 'ub', 'direct_effect', 'steps')
    
    return(out)
  })
 
  out = do.call(rbind,dataplot_temp)
  
# Determine Johnson-Neyman point 
  # Look for sign flips in the upper boundary
  pos_ub = which(diff(sign(out[,'ub']))!=0)
  pos_lb = which(diff(sign(out[,'lb']))!=0) 
# Prepare and return output
    res=list(out = out, 
             Johnson_Neyman_points = list(JNlb = c(x = out[,'steps'][pos_lb], y = out[,'lb'][pos_lb]),
                                          JNub = c(x = out[,'steps'][pos_ub], y = out[,'ub'][pos_ub])),
             Information = c(alpha=.alpha, independent = .independent, moderator= .moderator, dependent = .dependent))
    class(res) = c("Two_Way_Effect", class(res))
    return(res)
} 



#' Plot the output of floodlight analysis
#'
#' It plots the direct effect of an independent variable depending on the levels of the 
#' moderator variable. Moreover, the confidence interval are plotted including the 
#' Johnson-Neyman points. 
#'
#' @usage plot(.TWobject)
#'
#'
#' @export
plot.Two_Way_Effect = function(.TWobject){
  
  require(ggplot2)
  plot1=ggplot(as.data.frame(.TWobject$out),aes(x=.TWobject$out[,'steps'],y=.TWobject$out[,'direct_effect']))+
    geom_line()+
    geom_ribbon(aes(ymin=.TWobject$out[,'lb'], ymax=.TWobject$out[,'ub']),alpha=0.2)+
    labs(x=paste('Level of ',.TWobject$Information['moderator']) , 
    y=paste('Effect of', .TWobject$Information['independent'], 'on \n', .TWobject$Information['dependent']))+
    theme_bw()+
    # scale_x_continuous(breaks=seq(-3,3,0.5))+
    theme(panel.grid.minor = element_blank())
  
  # Plot Johnson-Neyman points, if they exist in the considered range
  if(length(.TWobject$Johnson_Neyman_points$JNlb)==2){
    JN=.TWobject$Johnson_Neyman_points$JNlb
    plot1 = plot1+
      geom_point(x=JN['x'],y=JN['y'],size=2)  
  }
  
  if(length(.TWobject$Johnson_Neyman_points$JNub)==2){
    JN=.TWobject$Johnson_Neyman_points$JNub
    plot1 = plot1+
      geom_point(x=JN['x'],y=JN['y'],size=2)  
  }
  
  # Plot
  plot1
}
  


#' Regression-based Hausman-test
#'
#' The regression based Hausman test
#'
#' @export
#' 
testHausman2 <- function(
  .object               = NULL,
  .alpha                = args_default()$.alpha,
  .eval_plan            = args_default()$.eval_plan,
  .handle_inadmissibles = args_default()$.handle_inadmissibles,
  .R                    = args_default()$.R,
  .resample_method      = args_default()$.resample_method,
  .seed                 = args_default()$.seed
  ) {
  
  # Define function to bootstrap
  fun_residuals <- function(.object) {

    ## Extract relevant quantities
    m  <- .object$Information$Model
    P  <- .object$Estimates$Construct_VCV
    
    # Dependent (LHS) variables of the structural equations
    dep_vars  <- rownames(m$structural)[rowSums(m$structural) != 0]
    
    res <- lapply(dep_vars, function(y) {
      # Which of the variables in dep_vars have instruments specified, i.e.
      # have endogenous variables on the RHS. By default: FALSE.
      endo_in_RHS <- FALSE
      
      if(!is.null(m$instruments)) {
        endo_in_RHS <- y %in% names(m$instruments)
      }
      
      # All independent variables (X) of the structural equation of construct y
      # (including the endogenous RHS variables)
      names_X <-  colnames(m$structural[y, m$structural[y, ] != 0, drop = FALSE])
      
      ## Only for equations with endogenous variables on the RHS do we need to compute
      ##  the Hausman test.
      if(endo_in_RHS & (.object$Information$Arguments$.approach_paths == "2SLS")) {
        ## First stage
        # Note: Technically, we only need to regress the P endogenous variables 
        #       (y1) on the L instruments and the K exogenous independent variables 
        #       (which must be part of Z).
        #       Therefore: y1 (N x P) and Z (N x (L + K)) and
        #       beta_1st = (Z'Z)^-1*(Z'y1)
        #       
        #       beta_1st would be ((L + K) x P), i.e. each columns represents the 
        #       first stage estimates for a regression of the instruments on
        #       on the p'th endogenous variable.
        names_endo <- rownames(m$instruments[[y]]) # names of the endogenous variables (y1's)
        names_Z <- colnames(m$instruments[[y]]) # names of the instruments 
        
        # Assuming that P (the construct correlation matrix) also contains 
        # the instruments (ensured if only internal instruments are allowed)
        
        beta_1st <- solve(P[names_Z, names_Z, drop = FALSE], 
                          P[names_Z, names_endo, drop = FALSE])
        
        ## Second stage
        # Note: we have to use (construct) correlations here since using 
        #       the plain score would yield inconsistent estimates for concepts
        #       modeled as common factors. Hence we need to "translate" the
        #       regression based second stage estimation 
        #                  y = X * beta_1 + v_hat*beta_2 + u
        #       and the estimator
        #                  beta_2nd = (W'W)^-1W'y
        #       into correlations.
        #
        #   y1_hat = Z*beta_1st
        #   v_hat  = y1 - y1_hat = y1 - Z*beta_1st
        #
        #   y = X * beta_1 + v_hat*beta_2 + u
        #
        #   Define W = [X; v_hat]
        #
        #   E(W'W ) = E[X'X ; X'v_hat
        #               v_hat'X ; v_hat'v_hat]
        #   E(W'y)  = E[X'y ; v_hat'y]'
        
        ww11 <- P[names_X, names_X, drop = FALSE]
        ww12 <- P[names_X, names_endo, drop = FALSE] - P[names_X, names_Z, drop = FALSE] %*% beta_1st
        ww21 <- t(ww12)
        ww22 <- P[names_endo, names_endo, drop = FALSE] - 
          P[names_endo, names_Z, drop = FALSE] %*% beta_1st -
          t(beta_1st) %*% P[names_Z, names_endo, drop = FALSE] +
          t(beta_1st) %*% P[names_Z, names_Z, drop = FALSE] %*% beta_1st
        
        WW <- rbind(cbind(ww11, ww12), cbind(ww21, ww22))
        
        wy1 <- P[names_X, y, drop = FALSE]
        wy2 <- P[names_endo, y, drop = FALSE] - t(beta_1st) %*% P[names_Z, y, drop = FALSE]

        Wy <- rbind(wy1, wy2)
        
        ## Estimate the second stage estimation

        beta_2nd <- solve(WW, Wy)
        rownames(beta_2nd) <- c(names_X, paste0("Resid_", names_endo))
        colnames(beta_2nd) <- y
        
        return(beta_2nd)
        
      } else {
        NA
      }
    })
    
    names(res) <- dep_vars
    res        <- Filter(Negate(anyNA), res) 
    
    # Vectorize "res" to be able to use resamplecSEMResults (which requires
    # output to be a vector or a matrix)
    # 1. Assign a unique name to each element
    res <- lapply(res, function(x) {
      rownames(x) <- paste0(colnames(x), "_", rownames(x))
      x
      })
    # 2. Combine, vectorize and assign names
    res2 <- do.call(rbind, res)
    res <- c(res2)
    names(res) <- rownames(res2)
    res
  }

  ## Resample to get standard errors
  out <- resamplecSEMResults(
    .object               = .object,
    .resample_method      = .resample_method,
    .R                    = .R,
    .handle_inadmissibles = .handle_inadmissibles,
    .user_funs            = fun_residuals,
    .eval_plan            = .eval_plan,
    .seed                 = .seed
    )
  # 
  # class(out) <- c("cSEMTestHausman", "cSEMResults_resampled")
  # out
  
  ## Get relevant quantities
  out_infer <- infer(out, .alpha = .alpha)
  beta      <- out$Estimates$Estimates_resample$Estimates1$User_fun$Original
  se        <- out_infer$User_fun$sd
  t         <- beta/se
  p_normal  <- 2*pnorm(abs(t), mean = 0, sd = 1, lower.tail = FALSE) 
  ci_percentile <- out_infer$User_fun$CI_percentile 
  
  out_data_frame <- data.frame(
    "Name"       = names(beta),
    "Estimate"   = beta,
    "Std. error" = se,
    "t-stat."    = t,
    "p-value"    = p_normal,
    "Ci_perc. L" = ci_percentile[1, ],
    "Ci_perc. H" = ci_percentile[2, ],
    stringsAsFactors = FALSE
  )
  rownames(out_data_frame) <- NULL
  return(out_data_frame)
}
