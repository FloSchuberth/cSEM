#' Do a nonlinear effects analysis
#'
#' Calculate the expected value of the dependent variable conditional on the values of 
#' an independent variables and a moderator variable. All other variables in the model
#' are assumed to be zero, i.e., they are fixed at their mean levels. Moreover, it produces
#' the input for the floodlight analysis.
#' 
#' @usage doNonlinearEffectsAnalysis(
#'  .object           = NULL,
#'  .dependent        = NULL, 
#'  .independent      = NULL,
#'  .moderator        = NULL,
#'  .n_steps          = 100,
#'  .values_moderator = c(-2,-1,0,1,2),
#'  .alpha            = 0.05
#'  )
#'
#' @inheritParams csem_arguments
#'
#' @references
#'   \insertAllCited{}
#' 
#' @return A list of class `cSEMNonlinearEffects` with a corresponding method 
#'   for `plot()`. 
#'   See: [plot.cSEMNonlinearEffects()].
#' 
#' @seealso [csem()], [cSEMResults], [plot.cSEMNonlinearEffects()]
#' 
#' @example inst/examples/example_doNonlinearEffectsAnalysis.R
#' 
#' @export

doNonlinearEffectsAnalysis <- function(
  .object             = NULL,
  .dependent          = NULL, 
  .independent        = NULL,
  .moderator          = NULL,
  .n_steps            = 100,
  .values_moderator   = c(-2,-1,0,1,2),
  .alpha              = 0.05
){
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, doSurfaceAnalysis, 
                  .dependent = .dependent, .independent = .independent,
                  .moderator   = .moderator, .n_steps = .n_steps)
    
    class(out) <- c("cSEMNonlinearEffects", "cSEMNonlinearEffects_multi")
    return(out)
  }else {
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
    Q  <- .object$Estimates$Reliabilities
  } 
  else {
    m   <- .object$Second_stage$Information$Model
    est <- .object$Second_stage$Information$Resamples$Estimates$Estimates1$Path_estimates
    H   <- .object$Second_stage$Estimates$Construct_scores
    Q <- .object$Second_stage$Estimates$Reliabilities
  }
  
  # Character string containing the names of the dependent and independent variables
  dep_vars   <- rownames(m$structural[rowSums(m$structural) != 0, , drop = FALSE])
  indep_vars <- colnames(m$structural[, colSums(m$structural[.dependent, ,drop = FALSE]) != 0 , drop = FALSE])
  
  ## Check if model is nonlinear.
  if(m$model_type != "Nonlinear"){
    stop2(
      "The following error occured in the `doNonlinearEffectsAnalysis()` function:\n",
      "The structural model must contain contain nonlinear terms.")
  }
  
  ## Works only for one significance level, i.e., no vector of significances is allowed
  if(length(.alpha) != 1){
    stop2(
      "The following error occured in the `doNonlinearEffectsAnalysis()` function:\n",
      "Currently only a single significance level (.alpha) is allowed.")
  }
  
  
  ## Check whether dependent, and two independent  variables are provided
  if(is.null(.dependent) | is.null(.independent) | is.null(.moderator)){
    stop2(
      "The following error occured in the `doNonlinearEffectsAnalysis()` function:\n",
      "All variables (.dependent, .independent , and .moderator) must be supplied.")
  }
  
  
  # Check if the name of the dependent variable is valid
  if(!(.dependent %in% dep_vars)){
    stop2(
      "The following error occured in the `doNonlinearEffectsAnalysis()` function:\n",
      "The dependent variable supplied to `.dependent` is not a dependent variable in the original model.")
  }
  
  # Check if the name of the moderator (.moderator) and the dependent variable (.independent) supplied are used 
  # in the original model
  if(!all(c(.moderator, .independent) %in% colnames(m$structural))){
    stop2(
      "The following error occured in the `doNonlinearEffectsAnalysis()` function:\n",
      "Independent and/or moderator variable are not part of the original model.")
  }
  
  ### Calculation Floodlight Analysis--------------------------------------------
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
  steps_flood <- seq(min(H[, .moderator]), max(H[, .moderator]), length.out = .n_steps)
  
  dataplot_temp_flood <- lapply(steps_flood, function(step){
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
  
  res_flood <- do.call(rbind, dataplot_temp_flood)

  colnames(res_flood) <- c('direct_effect', 'value_z', 'lb', 'ub')
  
  # Determine Johnson-Neyman point 
  # Look for sign flips in the upper boundary
  pos_ub <- which(diff(sign(res_flood[, 'ub'])) != 0)
  pos_lb <- which(diff(sign(res_flood[, 'lb'])) != 0) 
  
  out_flood <- list(
    "out"                   = res_flood, 
    "Johnson_Neyman_points" = list(
      JNlb = c(x = res_flood[,'value_z'][pos_lb], 
               y = res_flood[,'lb'][pos_lb]),
      JNub = c(x = res_flood[,'value_z'][pos_ub], 
               y = res_flood[,'ub'][pos_ub]
      )
    )
  )
  
  
  # Calculate information necessary for print -----------------------------
  dataplot_temp_print <- lapply(c(-2,-1,0,1,2), function(step){
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
  
  out_print <- do.call(rbind, dataplot_temp_print)
  
  colnames(out_print) <- c('direct_effect', 'value_z', sprintf("%.6g%%L", 1-.alpha),
                           sprintf("%.6g%%U", 1-.alpha))
  
  
  
  ### Calculation Simple effects and surface analysis---------------------------------------
  if(inherits(.object, "cSEMResults_default")) {
    sum_object = summarize(.object = .object)
  } else if(inherits(.object, "cSEMResults_2ndorder")){
    sum_object = summarize(.object = .object$Second_stage)
  }

  # Effects all
  effects_all <- sum_object$Estimates$Path_estimates
  
  nonlinear_vars <- indep_vars[grepl(pattern = '.',x = indep_vars,fixed = T)]
  
  # interactions of the .independent and the .moderator 
  pointer=sapply(indep_vars,function(x){
    temp=unlist(strsplit(x,'\\.'))
    if(sum(grepl(paste(c(.moderator, .independent),collapse = '|'  ),temp))==length(temp)){
      TRUE
    }else{
      FALSE
    }
  })
  
  if(!all(pointer)){
    warning2(paste0("The considered equation contains the following variables that do not\n",
                    "only involve ",.independent, " and ", .moderator, ":\n\n",
                    paste(indep_vars[!pointer],collapse =', '),"\n\n",
                    "They will be ignored in calculating the predicted values of ", .dependent,"i.e., \n",
                    "the other variables are considered at their mean values."))
  }
  
  vars_rel=indep_vars[pointer]
  
  sum_rel=sum_object$Estimates$Path_estimates[sum_object$Estimates$Path_estimates$Name%in% 
                                                paste(.dependent,"~",vars_rel),]
  
  steps_ind_simple = matrix(seq(min(H[, .independent ]), max(H[, .independent ]),
                         length.out = .n_steps),ncol=1)
  colnames(steps_ind_simple)=.independent
  
  # Calculation of simple effects analysis
  steps_mod_simple = .values_moderator

  
  y_pred_list_simple=lapply(steps_mod_simple,function(stepmode){
    temp<-lapply(sum_rel$Name,function(x){
      temp <- sum_rel[sum_rel$Name==x,]
      
      indep_names <- unlist(strsplit(temp$Name,paste(.dependent,"~ ")))[2]
      indep_rel <- unlist(strsplit(indep_names,'\\.'))
      
      # Check which is the moderator and which is the independent variable
      moderat <- indep_rel[.moderator==indep_rel]
      indep <- indep_rel[.independent==indep_rel]
      
      #Reliabilities of composites are 1 so nothing happens, except in case of 
      # SOC of the type composite of common factor
      
      # The formula for the average is fine, as the expected value are only required 
      # for terms with a mximum power of 3. In this case, the formula is appropriate
      average <- mean(matrixStats::rowProds(H[,indep_rel,drop=FALSE]))/prod(Q[indep_rel])
      
      temp$Estimate*(matrixStats::rowProds(steps_ind_simple[,indep,drop=FALSE])*
                       stepmode^length(moderat)-
                       average)  
    }
    )
    temp
    Reduce('+',temp)
  })
  
  out_simpleeffects <- data.frame(values_dep = unlist(y_pred_list_simple),
                    values_ind = rep(steps_ind_simple,length(steps_mod_simple)),
                    values_mod = paste0(rep(steps_mod_simple,each=length(steps_ind_simple))))
  
  # Calculation of surface analsyis
  steps_ind_surface = seq(min(H[, .independent ]), max(H[, .independent ]),
                   length.out = .n_steps)
  steps_mod_surface = seq(min(H[, .moderator ]), max(H[, .moderator ]),
                          length.out = .n_steps)
  
  steps_independent = rep(steps_ind_surface,each=length(steps_mod_surface))
  steps_moderator = rep(steps_mod_surface,times = length(steps_ind_surface))
  
  steps_surface=cbind(steps_independent,steps_moderator)
  colnames(steps_surface)=c(.independent ,.moderator )
  
  y_pred_list_surface=lapply(sum_rel$Name,function(x){
    temp <- sum_rel[sum_rel$Name==x,]
    
    indep_names <- unlist(strsplit(temp$Name,paste(.dependent,"~ ")))[2]
    indep_rel <- unlist(strsplit(indep_names,'\\.'))
    
    #Reliabilities of composites are 1 so nothing happens, except in case of 
    # SOC of the type composite of common factor
    
    # The formula for the average is fine, as the expected value are only required 
    # for terms with a mximum power of 3. In this case, the formula is appropriate
    average <- mean(matrixStats::rowProds(H[,indep_rel,drop=FALSE]))/prod(Q[indep_rel])
    
    temp$Estimate*(matrixStats::rowProds(steps_surface[,indep_rel,drop=FALSE])-average)  
  }
  )
  
  y_pred_surface = Reduce('+',y_pred_list_surface)
  
  out_surface <- list(values_dep=matrix(y_pred_surface,ncol=length(steps_ind_surface),nrow=length(steps_mod_surface)),
              values_ind1=steps_ind_surface,values_ind2=steps_mod_surface)
  
  # Prepare and return output
  out <- list(
    "out_floodlight"  = out_flood,
    "out_surface" = out_surface,
    "out_simpleeffects" = out_simpleeffects,  
    # "out1" = out1,
    "Information" = list(
      dependent       = .dependent,
      independent   = .independent,
      moderator   = .moderator,
      all_independent = vars_rel,
      values_moderator = .values_moderator,
      alpha = .alpha
    ),
    "Information_print" = out_print
  )
  
  class(out) <- "cSEMNonlinearEffects"
  return(out)
  }
}

