#' Do a floodlight analysis
#'
#' Calculate the effect of an independent variable (z) on a dependent variable
#' (y) conditional on the values of a (continous) moderator variable (x) 
#' to perform a floodlight analysis \insertCite{Spiller2013}{cSEM}. Moreover, 
#' the Johnson-Neyman points are calculated, i.e. the value(s) of x for which 
#' lower or upper boundary of the confidence interval of the effect
#' estimate of z for a given x switch signs. 
#' 
#' @usage doFloodlightAnalysis(
#'  .object         = NULL,
#'  .alpha          = 0.05,
#'  .dependent      = NULL, 
#'  .moderator      = NULL,
#'  .independent    = NULL,
#'  .n_steps        = 100
#'  )
#'
#' @inheritParams csem_arguments
#'
#' @references
#'   \insertAllCited{}
#' 
#' @return A list of class `cSEMFloodlight` with a corresponding method for `plot()`. 
#'   See: [plot.cSEMFloodlight()].
#' 
#' @seealso [csem()], [cSEMResults], [plot.cSEMFloodlight()]
#' @export

doFloodlightAnalysis <- function(
  .object        = NULL,
  .alpha         = 0.05,
  .dependent     = NULL, 
  .moderator     = NULL,
  .independent   = NULL,
  .n_steps  = 100
){
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, doFloodlightAnalysis, .alpha = .alpha, 
                  .dependent = .dependent, .moderator = .moderator, .independent = .independent, .n_steps = .n_steps)
    
    class(out) <- c("cSEMFloodlight", "cSEMFloodlight_multi")
    return(out)
  } else {
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
    indep_vars <- colnames(m$structural[, colSums(m$structural[.dependent, ,drop = FALSE]) !=0 , drop = FALSE])
    
    ## Check if model is nonlinear.
    if(m$model_type != "Nonlinear"){
      stop2(
        "The following error occured in the `doFloodlightAnalysis()`` function:\n",
        "The structural model must contain contain interaction terms.")
    }
    
    ## Works only for one significance level, i.e., no vector of significances is allowed
    if(length(.alpha) != 1){
      stop2(
        "The following error occured in the `doFloodlightAnalysis()`` function:\n",
        "Currently only a single significance level (.alpha) is allowed.")
    }
    
    ## Check whether dependent, independent, and moderator variable are provided
    if(is.null(.dependent) | is.null(.moderator) | is.null(.independent)){
      stop2(
        "The following error occured in the `doFloodlightAnalysis()`` function:\n",
        "All variables (.dependent, .moderator, and .independent) must be supplied.")
    }
    
    # Check if the name of the dependent variable is valid
    if(!(.dependent %in% dep_vars)){
      stop2(
        "The following error occured in the `doFloodlightAnalysis()`` function:\n",
        "The dependent variable supplied to `.dependent` is not a dependent variable in the original model.")
    }
    
    # Check if the name of the moderator (.moderator) and the dependent variable (.independent) supplied are used 
    # in the original model
    if(!all(c(.moderator, .independent) %in% colnames(m$structural))){
      stop2(
        "The following error occured in the `doFloodlightAnalysis()`` function:\n",
        "Independent and/or moderator variable are not part of the original model.")
    }
    
    ### Calculation --------------------------------------------------------------
    # Possible names of the interaction terms
    possible_names <- c(paste(.moderator, .independent, sep = '.'), paste(.independent, .moderator, sep= '.'))
    
    # Name of the interaction term
    xz <- possible_names[possible_names %in% indep_vars]
    
    if(length(xz) != 1){
      stop2(
        "The defined interaction term does not exist in the model or is not",
        " part of the equation of the dependent variable.")
    }
    
    
    # Effect names
    beta_z  <- paste(.dependent, .independent, sep = ' ~ ')
    
    beta_xz <- paste(.dependent, xz, sep = ' ~ ')
    
    ## Compute spotlights (= effects of independent (z) on dependent (y) for given 
    ## values of moderator (x))
    steps <- seq(min(H[, .moderator]), max(H[, .moderator]), length.out = .n_steps)
    
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
        dependent   = .dependent,
        independent = .independent,
        moderator   = .moderator
      )
    )
    
    class(out) <- "cSEMFloodlight"
    return(out)
  }
} 

#' Do a surface analysis
#'
#' Calculate the effect of an independent variable (z) on a dependent variable
#' (y) including their moderation and the respective quadratic terms.
#' 
#' @usage doSurfaceAnalysis(
#'  .object             = NULL,
#'  .alpha              = 0.05,
#'  .dependent          = NULL, 
#'  .independent        = NULL,
#'  .moderator          = NULL,
#'  .n_steps       = 100
#'  )
#'
#' @inheritParams csem_arguments
#'
#' @references
#'   \insertAllCited{}
#' 
#' @return A list of class `cSEMSurface` with a corresponding method for `plot()`. 
#'   See: [plot.cSEMSurface()].
#' 
#' @seealso [csem()], [cSEMResults], [plot.cSEMFloodlight()]
#' @export

doSurfaceAnalysis <- function(
    .object             = NULL,
    .alpha              = 0.05,
    .dependent          = NULL, 
    .independent        = NULL,
    .moderator          = NULL,
    .n_steps       = 100
){
  

  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, doSurfaceAnalysis, .alpha = .alpha, 
                  .dependent = .dependent, .independent = .independent,
                  .moderator = .moderator, .square = .square, .n_steps = .n_steps)
    
    class(out) <- c("cSEMFloodlight", "cSEMFloodlight_multi")
    return(out)
  } 
  # else {
    ## Check whether .object is of class cSEMResults_resampled; if not perform
    ## standard resampling (bootstrap with .R = 499 reps)
    # if(!inherits(.object, "cSEMResults_resampled")) {
    #   if(inherits(.object, "cSEMResults_default")) {
    #     args <- .object$Information$Arguments
    #   } else {
    #     args <- .object$Second_stage$Information$Arguments_original
    #   }
    #   args[".resample_method"] <- "bootstrap"
    #   .object <- do.call(csem, args)
    # }
    
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
    
    sum_object = summarize(.object = .object)
  
    # Character string containing the names of the dependent and independent variables
    dep_vars   <- rownames(m$structural[rowSums(m$structural) !=0, , drop = FALSE])
    indep_vars <- colnames(m$structural[, colSums(m$structural[.dependent, ,drop = FALSE]) !=0 , drop = FALSE])
  
    # Effects all
    effect_all <- sum_object$Estimates$Path_estimates
    ## Check if model is nonlinear.
    if(m$model_type != "Nonlinear"){
      stop2(
        "The following error occured in the `doFloodlightAnalysis()`` function:\n",
        "The structural model must contain contain nonlinear terms.")
    }
    
    ## Works only for one significance level, i.e., no vector of significances is allowed
    if(length(.alpha) != 1){
      stop2(
        "The following error occured in the `doFloodlightAnalysis()`` function:\n",
        "Currently only a single significance level (.alpha) is allowed.")
    }
    
    
    # CHECKEN OB DAS NOCH SO PASST
    ## Check whether dependent, independent, and moderator variable are provided
    if(is.null(.dependent) | is.null(.moderator) | is.null(.independent)){
      stop2(
        "The following error occured in the `doFloodlightAnalysis()`` function:\n",
        "All variables (.dependent, .moderator, and.independent) must be supplied.")
    }
    

    # Check if the name of the dependent variable is valid
    if(!(.dependent %in% dep_vars)){
      stop2(
        "The following error occured in the `doFloodlightAnalysis()`` function:\n",
        "The dependent variable supplied to `.dependent` is not a dependent variable in the original model.")
    }
    
    # Check if the name of the moderator (.x) and the dependent variable (.z) supplied are used 
    # in the original model
    if(!all(c(.moderator,.independent) %in% colnames(m$structural))){
      stop2(
        "The following error occured in the `doFloodlightAnalysis()`` function:\n",
        "Independent and/or moderator variable are not part of the original model.")
    }
    
    ### Calculation --------------------------------------------------------------
    
    ## Filter interaction, and quadratic terms:
    # Select the terms that either include x or z
    
    # if there are terms that include not only x and z but a third variable stop

    
    # Output: All interaction terms where the dependent and the moderator are 
    # included including quadratic and higher order polynomials
    
    
    # Select the corresponidng effects
    
    
    
    # Effect name of the independent variable
    beta_z  <- paste(.dependent,.independent, sep = ' ~ ')
    
    
    # Possible names of the interaction terms
    if(!is.null(.moderator)){
      # Effect name of the moderator variable
      beta_x  <- paste(.dependent,.moderator, sep = ' ~ ')
      
      possible_names_int <- c(paste(.moderator,.independent, sep = '.'), paste(.independent, .moderator, sep= '.'))
     
    # Name of the interaction term
    xz <- possible_names_int[possible_names_int %in% indep_vars]
    
    # if(length(xz) != 1){
    #   stop2(
    #     "The defined interaction term does not exist in the model or is not",
    #     " part of the equation of the dependent variable.")
    # }
    
    # effect name
    beta_xz <- paste(.dependent, xz, sep = ' ~ ')
    
    } 
    
    
      # name of the quadratic term
      zz <- paste(.quadratic,.quadratic, sep = '.')
      beta_zz <- paste(.dependent, xz, sep = ' ~ ')
      }
    
    # Hier eventuell adjustment von aussen zulassen 
    ## Compute steps for the moderator and the independent variable
    # steps <- list(
      steps_independent = seq(min(H[, .independent]), max(H[, .independent]), 
                              length.out = .n_steps)
      steps_moderator = seq(min(H[, .moderator]), max(H[, .moderator]), 
                            length.out = .n_steps)

      # Calculate the predicted values of the dependent variable
      predicted_dependent <- est$Resampled[ , beta_z] * steps_independent+
        est$Resampled[, beta_x]  * steps_moderator+
        est$Resampled[, beta_xz]  * (steps_moderator *steps_independent 
                                     - mean(steps_moderator *steps_independent))+
        est$Resampled[, beta_zz]  * (steps_independent^2-1)+
        est$Resampled[, beta_xx]  * (steps_moderator^2-1)
      
      # Value of the originally estimated effect of z on y at different levels of x
      # effect_original <- est$Original[beta_z] + est$Original[beta_xz] * step
      
      # Compute empirical quantile based on resamples
      # bounds <- quantile(effect_resampled , c(.alpha/2, 1 - .alpha/2))
      
      # Return output
      ret <- cbind(y=predicted_dependent,z=steps_independent,x=steps_moderator)


    # Determine Johnson-Neyman point 
    # Look for sign flips in the upper boundary
    # pos_ub <- which(diff(sign(out[, 'ub'])) != 0)
    # pos_lb <- which(diff(sign(out[, 'lb'])) != 0) 
    
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
        dependent   = .dependent,
        independent =.independent,
        moderator   = .moderator
      )
    )
    
    class(out) <- "cSEMSurface"
    return(out)
  }
} 