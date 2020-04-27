#' Do a simple effects analysis
#'
#' Calculate the effect of an independent variable on a dependent variable
#' conditional on the values of a (continous) moderator variable. All other variables in the model
#' are assumed to be zero, i.e., at their mean level. 
#' 
#' @usage doSimpleEffectsAnalysis(
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
#' @return A list of class `cSEMSimpleEffects` with a corresponding method for `plot()`. 
#'   See: [plot.cSEMFloodlight()].
#' 
#' @seealso [csem()], [cSEMResults], [plot.cSEMFloodlight()]
#' @export

doSimpleEffectsAnalysis <- function(
  .object             = NULL,
  .dependent          = NULL, 
  .independent      = NULL,
  .moderator      = NULL,
  .n_steps            = 100
){
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, doSurfaceAnalysis, 
                  .dependent = .dependent, .independent = .independent,
                  .moderator   = .moderator, .n_steps = .n_steps)
    
    class(out) <- c("cSEMSurface", "cSEMSurface_multi")
    return(out)
  } 
  
  # Might be relevant for upcoming CIs
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
    Q  <- .object$Estimates$Reliabilities
  } 
  else {
    m   <- .object$Second_stage$Information$Model
    est <- .object$Second_stage$Information$Resamples$Estimates$Estimates1$Path_estimates
    H   <- .object$Second_stage$Estimates$Construct_scores
    Q <- .object$Second_stage$Estimates$Reliabilities
  }
  
  ## Check if model is nonlinear.
  if(m$model_type != "Nonlinear"){
    stop2(
      "The following error occured in the `doSimpleEffectsAnalysis()`` function:\n",
      "The structural model must contain contain nonlinear terms.")
  }
  

  ## Check whether dependent, and two independent  variables are provided
  if(is.null(.dependent) | is.null(.independent) | is.null(.moderator)){
    stop2(
      "The following error occured in the `doSimpleEffectsAnalysis()` function:\n",
      "All variables (.dependent, .independent , and .moderator) must be supplied.")
  }
  
  
  # Check if the name of the dependent variable is valid
  if(!(.dependent %in% rownames(m$structural[rowSums(m$structural) !=0, , drop = FALSE]))){
    stop2(
      "The following error occured in the `doSimpleEffectsAnalysis()` function:\n",
      "The dependent variable supplied to `.dependent` is not a dependent variable in the original model.")
  }
  
  # Check if the names of the two independent variables (.independent_1 and .independent_2)
  # supplied are used in the original model
  if(!all(c(.independent ,.moderator) %in% colnames(m$structural))){
    stop2(
      "The following error occured in the `doSimpleEffectsAnalysis()` function:\n",
      "Either the independent variable or the moderator are not part of the original model.")
  }
  
  ### Calculation --------------------------------------------------------------
  if(inherits(.object, "cSEMResults_default")) {
    sum_object = summarize(.object = .object)
  } else if(inherits(.object, "cSEMResults_2ndorder")){
    sum_object = summarize(.object = .object$Second_stage)
  }
  # Character string containing the names of the dependent and independent variables
  dep_vars   <- rownames(m$structural[rowSums(m$structural) !=0, , drop = FALSE])
  indep_vars <- colnames(m$structural[, colSums(m$structural[.dependent, ,drop = FALSE]) !=0 , drop = FALSE])
  
  # Effects all
  effects_all <- sum_object$Estimates$Path_estimates
  
  nonlinear_vars=indep_vars[grepl(pattern = '.',x = indep_vars,fixed = T)]
  
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
                    "They will be ignored in calculating the predicted values of ", .dependent,".\n",
                    "Hence, the other variables are considered at their mean values."))
  }
  
  vars_rel=indep_vars[pointer]
  
  sum_rel=sum_object$Estimates$Path_estimates[sum_object$Estimates$Path_estimates$Name%in% 
                                                paste(.dependent,"~",vars_rel),]
  
  steps_ind = matrix(seq(min(H[, .independent ]), max(H[, .independent ]),
                   length.out = .n_steps),ncol=1)
  colnames(steps_ind)=.independent
  
  steps_mod = c(-2,-1,0,1,0)

  # steps_independent1 = rep(steps_ind1,each=length(steps_ind2))
  # steps_independent2 = rep(steps_ind2,times = length(steps_ind1))
  
  # steps=cbind(steps_independent1,steps_independent2)
  # colnames(steps)=c(.independent_1 ,.independent_2 )
  
  y_pred_list=lapply(steps_mod,function(stepmode){
    temp<-lapply(sum_rel$Name,function(x){
    temp <- sum_rel[sum_rel$Name==x,]
    
    indep_names <- unlist(strsplit(temp$Name,paste(.dependent,"~ ")))[2]
    indep_rel <- unlist(strsplit(indep_names,'\\.'))
    
    # Check which is the moderator and which is the independent variable
    moderat <- indep_rel[.moderator==indep_rel]
    indep <- indep_rel[.independent==indep_rel]
    
    #Reliabilities of composites are 1 so nothing happens, except in case of 
    # SOC of the type composite of common factor
    
    # !!! It might be that the calculation for the expected value needs to be adjusted.
    average <- mean(matrixStats::rowProds(H[,indep_rel,drop=FALSE]))/prod(Q[indep_rel])
    
    temp$Estimate*(matrixStats::rowProds(steps_ind[,indep,drop=FALSE])*
                     steps_mod^length(moderat)-
                     average)  
  }
  )
    temp
  })
  
  y_pred = Reduce('+',y_pred_list)
  
  # Return output
  # The returned matrix has in the columns the value of the independent variable
  # and in the rows the values of the moderator
  out <- list(values_dep=matrix(y_pred,ncol=length(steps_ind1),nrow=length(steps_ind2)),
              values_ind1=steps_ind1,values_ind2=steps_ind2)
  
  # out1 <- list(y=y_pred, z=steps_independent1,x=steps_independent2)
  
  # Prepare and return output
  out <- list(
    "out"  = out,
    # "out1" = out1,
    "Information" = list(
      alpha           = .alpha,
      dependent       = .dependent,
      independent_1   = .independent_1,
      independent_2   = .independent_2,
      all_independent = vars_rel 
    )
  )
  
  class(out) <- "cSEMSimpleEffects"
  return(out)
}

