#' Do a surface analysis
#'
#' Based on a nonlinear model, the dependent variable of a certain equation is
#' predicted by a certain independent variable and a certain second variable (moderator)
#' including their higher-order terms.
#' 
#' @usage doSurfaceAnalysis(
#'  .object             = NULL,
#'  .alpha              = 0.05,
#'  .dependent          = NULL, 
#'  .independent        = NULL,
#'  .moderator          = NULL,
#'  .n_steps            = 100
#'  )
#'
#' @inheritParams csem_arguments
#'
#' @return A list of class `cSEMSurface` with a corresponding method for `plot()`. 
#'   See: [plot.cSEMSurface()].
#' 
#' @seealso [csem()], [cSEMResults], [plot.cSEMSurface()]
#' @export

doSurfaceAnalysis <- function(
  .object             = NULL,
  .alpha              = 0.05,
  .dependent          = NULL, 
  .independent        = NULL,
  .moderator          = NULL,
  .n_steps            = 100
){
  
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, doSurfaceAnalysis, .alpha = .alpha, 
                  .dependent = .dependent, .independent = .independent,
                  .moderator = .moderator, .n_steps = .n_steps)
    
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
  } else {
    m   <- .object$Second_stage$Information$Model
    est <- .object$Second_stage$Information$Resamples$Estimates$Estimates1$Path_estimates
    H   <- .object$Second_stage$Estimates$Construct_scores
    Q <- .object$Second_stage$Estimates$Reliabilities
    }
  
  ## Check if model is nonlinear.
  if(m$model_type != "Nonlinear"){
    stop2(
      "The following error occured in the `doSurfaceAnalysis()`` function:\n",
      "The structural model must contain contain nonlinear terms.")
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
      "All variables (.dependent, .moderator, and.independent) must be supplied.")
  }
  
  
  # Check if the name of the dependent variable is valid
  if(!(.dependent %in% rownames(m$structural[rowSums(m$structural) !=0, , drop = FALSE]))){
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
  sum_object = summarize(.object = .object)
  
  # Character string containing the names of the dependent and independent variables
  dep_vars   <- rownames(m$structural[rowSums(m$structural) !=0, , drop = FALSE])
  indep_vars <- colnames(m$structural[, colSums(m$structural[.dependent, ,drop = FALSE]) !=0 , drop = FALSE])
  
  # Effects all
  effects_all <- sum_object$Estimates$Path_estimates
  
  nonlinear_vars=indep_vars[grepl(pattern = '.',x = indep_vars,fixed = T)]
  
  # Grab quadratic independent
  quad_ind=nonlinear_vars[grepl(paste0(.independent,'.',.independent),nonlinear_vars)]
  quad_mod=nonlinear_vars[grepl(paste0(.moderator,'.',.moderator),nonlinear_vars)]
  
  # interactions of the .independent and the .moderator
  pointer=sapply(indep_vars,function(x){
    temp=unlist(strsplit(x,'\\.'))
    if(sum(grepl(paste(c(.moderator, .independent),collapse = '|'  ),temp))==length(temp)){
      TRUE
    }else{
      FALSE
    }
  })
  
  vars_rel=indep_vars[pointer]
  
  sum_rel=sum_object$Estimates$Path_estimates[sum_object$Estimates$Path_estimates$Name%in% 
                                        paste(.dependent,"~",vars_rel),]
  
  steps_ind = seq(min(H[, .independent]), max(H[, .independent]), 
                  length.out = .n_steps)
  steps_mod = seq(min(H[, .moderator]), max(H[, .moderator]), 
                  length.out = .n_steps)
  
  steps_independent = rep(steps_ind,each=length(steps_mod))
  steps_moderator = rep(steps_mod,times = length(steps_ind))
  
  steps=cbind(steps_independent,steps_moderator)
  colnames(steps)=c(.independent,.moderator)
  
  y_pred_list=lapply(sum_rel$Name,function(x){
    temp <- sum_rel[sum_rel$Name==x,]
  
    indep_names <- unlist(strsplit(temp$Name,paste(.dependent,"~ ")))[2]
    indep_rel <- unlist(strsplit(indep_names,'\\.'))
    
    #Reliabilities of composites are 1 so nothing happens, except in case of 
    # SOC of the type composite of common factor
    average <- mean(matrixStats::rowProds(H[,indep_rel,drop=FALSE]))/prod(Q[indep_rel])

    temp$Estimate*(matrixStats::rowProds(steps[,indep_rel,drop=FALSE])-average)  
    }
    )
  
  y_pred = Reduce('+',y_pred_list)
  
  # Return output
  out <- list(y=matrix(y_pred,nrow=length(steps_ind),ncol=length(steps_mod)),
              z=steps_ind,x=steps_mod)
  
  out1 <- list(y=y_pred, z=steps_independent,x=steps_moderator)
  
  # Prepare and return output
  out <- list(
    "out"                   = out,
    "out1" = out1,
    "Information" = list(
      alpha       = .alpha,
      dependent   = .dependent,
      independent =.independent,
      moderator   = .moderator,
      all_independent = vars_rel 
    )
  )
  
  class(out) <- "cSEMSurface"
  return(out)
}

