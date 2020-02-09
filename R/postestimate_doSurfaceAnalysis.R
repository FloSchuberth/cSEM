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
  effects_all <- sum_object$Estimates$Path_estimates
  
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
  
  
  
  # Select all relevant effect names 
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
  # }
  
  # Hier eventuell adjustment von aussen zulassen 
  ## Compute steps for the moderator and the independent variable
  # steps <- list(
  steps_ind = seq(min(H[, .independent]), max(H[, .independent]), 
                  length.out = .n_steps)
  steps_mod = seq(min(H[, .moderator]), max(H[, .moderator]), 
                  length.out = .n_steps)
  
  steps_independent = rep(steps_ind,each=length(steps_mod))
  steps_moderator = rep(steps_mod,times = length(steps_ind))
  # Create list with steps
  
  
  # Calculate the predicted values of the dependent variable
  predicted_dependent <- effects_all[effects_all$Name==beta_z,]$Estimate * steps_independent+
    effects_all[effects_all$Name==beta_x,]$Estimate  * steps_moderator
  
  effects_all[effects_all$Name==beta_xz,]$Estimate  * (steps_moderator *steps_independent 
                                                       - mean(steps_moderator *steps_independent))+
    effects_all[effects_all$Name==beta_zz,]$Estimate  * (steps_independent^2-1)+
    effects_all[effects_all$Name==beta_xx,]$Estimate  * (steps_moderator^2-1)
  
  # Value of the originally estimated effect of z on y at different levels of x
  # effect_original <- est$Original[beta_z] + est$Original[beta_xz] * step
  
  # Compute empirical quantile based on resamples
  # bounds <- quantile(effect_resampled , c(.alpha/2, 1 - .alpha/2))
  
  # Return output
  ret <- list(y=matrix(predicted_dependent,nrow=length(steps_ind),ncol=length(steps_mod)),
              z=steps_ind,x=steps_mod)
  
  
  # Determine Johnson-Neyman point 
  # Look for sign flips in the upper boundary
  # pos_ub <- which(diff(sign(out[, 'ub'])) != 0)
  # pos_lb <- which(diff(sign(out[, 'lb'])) != 0) 
  
  # # Prepare and return output
  # out <- list(
  #   "out"                   = out, 
  #   "Johnson_Neyman_points" = list(
  #     JNlb = c(x = out[,'value_z'][pos_lb], 
  #              y = out[,'lb'][pos_lb]),
  #     JNub = c(x = out[,'value_z'][pos_ub], 
  #              y = out[,'ub'][pos_ub]
  #     )
  #   ),
  #   # NEEDS TO BE DONE ADD ALL THE NAMES
  #   "Information" = list(
  #     alpha       = .alpha, 
  #     dependent   = .dependent,
  #     independent =.independent,
  #     moderator   = .moderator
  #   )
  # )
  
  class(ret) <- "cSEMSurface"
  return(ret)
}

# Needs to be added to dependencies
library(plotly)
# For plotting 3D figures options(viewer=NULL) must be set, perhaps there is a more elegant way

plot_ly( x = ret$x, y = ret$z, z = ret$y, type = "surface")

p1 <- plot_ly(x= ret$x, y=ret$y,z = ret$z, scene='scene1', lighting = list(ambient = 0.2)) %>%
  +   add_surface(showscale=FALSE)

p1 = plot_ly(ret, x = ~x, y = ~y, z = ~z) 
%>%
  +   add_surface(p1,showscale=FALSE)
