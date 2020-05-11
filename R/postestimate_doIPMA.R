#' Do an importance-performance matrix analysis
#'
#' Performs an importance-performance matrix analysis (IPMA).
#' 
#' To calculate the performance and importance, the weights of the indicators 
#' are unstandardized using the standard deviation of the original 
#' indicators but normed to have a length of 1.
#' Normed construct scores are calculated based on the original indicators and the 
#' unstandardized weights. 
#' 
#' The importance is calculated as the mean of 
#' the original indicators or the unstandardized construct scores, respectively. 
#' The performance is calculated as the unstandardized total effect if
#' `.level == "construct"` and as the normed weight times the unstandardized 
#' total effect if `.level == "indicator"`. The literature recommend to use an 
#' estimation approach as input for `[doIPMA()] that is based on normed 
#' indicators, e.g., by scaling all indicators to 0 to 100, 
#' see e.g., \insertCite{Henseler2020,Ringle2016;textual}{cSEM}.
#' 
#' Indicators are not normed internally, as theoretical maximum/minimum can 
#' differ from the empirical maximum/minimum which leads to an incorrect normalization.
#' 
#' @usage doIPMA(.object)
#'
#' @param .object A `cSEMResults` object.`
#'
#' @return A list of class `cSEMIPA` with a corresponding method for `plot()`. 
#'   See: [plot.cSEMIPMA()].
#' 
#' @seealso [csem()], [cSEMResults], [plot.cSEMIPMA()]
#' @export

doIPMA <- function(.object){
  
  # Checks ----
  # Check whether .object is a cSEMResults object
  if(inherits(.object, "cSEMResults_multi")) {
    out <- lapply(.object, doIPMA)
    return(out)
    class(out) <- c("cSEMIPMA", "cSEMIPMA_multi")
  } else if(inherits(.object, "cSEMResults_default")) {
    
  } else if(inherits(.object, "cSEMResults_2ndorder")){
    stop2("Currently, `doIPMA()` is not implemented for models containing", 
          " higher-order constructs.")
  } else {
    stop2("`.object` must be of class cSEMREsults.")
  }
  
  
  # Check whether Pearson was applied 
  if(.object$Information$Type_of_indicator_correlation != 'Pearson'){
    stop2("Currently, IPMA can only be applied if Pearson correlation was used.")
  }
  
  # Check whether PLS was applied
  if(.object$Information$Arguments$.approach_weights != 'PLS-PM'){
    stop2("Currently, IPMA can only be applied if PLS-PM was used as weighting scheme.")
  }
  
  # Check whether model is linear, i.e., no higher-order terms
  if(.object$Information$Model$model_type != 'Linear'){
    stop2("Currently, IPMA can only applied to linear models, i.e., models without higher-order terms.")
  }
  
  # Check whether model was estimated by OLS.
  if(.object$Information$Arguments$.approach_paths != 'OLS'){
    stop2("Currently, IPMA can only applied to models that were estimated by OLS.")
  }
  
  # Check whether all constructs are composite
  if(!all(.object$Information$Model$construct_type %in% "Composite")){
    warning2("At least one construct is modeled as common factor,\n",
    " IPMA's results are not trustworthy in this case.")
  }
  
  # Check whether all weights are positive
  if(!all(.object$Estimates$Weight_estimates >= 0)){
    warning2("At least one weight is negative,\n",
    " which makes the performance importance matrix analysis questionable.")
  }
  
  
  out_scores <- getConstructScores(.object = .object, .standardized = FALSE)
  
  # Peformance indicators ----
  performance_indicator <- colMeans(out_scores$Indicators_used)
  
  scores <- out_scores$Construct_scores

  # Calculate the performance values for constructs ----
  # In principle, we can allow for other location measures
  performance_construct <- colMeans(scores)

  # Calculate importance values for constructs, see Eq. 13.4 in Henseler (2020) ----
  effects <- calculateEffects(.object = .object, .output_type = 'matrix')
  total_effects_std <- effects$Total_effect
  
  # calculate unstandardized total effects
  total_effects_unstd <- total_effects_std
  for(i in rownames(total_effects_unstd)){
    for(j in colnames(total_effects_unstd)){
      total_effects_unstd[i,j] <- total_effects_std[i,j]* sd(scores[,i])/sd(scores[,j])
    }
  }

  importance_construct <- total_effects_unstd 

  # Obtain importance values for indicators ----
  importance_indicator <- out_scores$W_used
  for(ind in colnames(importance_indicator)){
    for(cons in rownames(importance_indicator)){
      # if the indicator is not connected to this construct, 
      # its normed weight is multiplied with the unstandardized total effect
      # Else the importance equals the normed weight
      if(.object$Information$Model$measurement[cons,ind] == 0){
        cons_names <- rownames(.object$Information$Model$measurement)
        ind_belongs <- cons_names[.object$Information$Model$measurement[,ind]==1]
        importance_indicator[cons,ind] <- out_scores$W_used[ind_belongs,ind]*total_effects_unstd[cons,ind_belongs]
      }  
    }
  }
  
  # collect output
  out<-list()
  out[["Construct"]]<-list("Importance"=importance_construct,
                           "Performance" = performance_construct)
  out[["Indicator"]]<-list("Importance"= importance_indicator,
                           "Performance" = performance_indicator)
  out[["Construct_names"]] <- c(.object$Information$Model$cons_exo,.object$Information$Model$cons_endo)
  out[["Indicator_names"]] <- .object$Information$Model$indicators

  class(out) <- "cSEMIPMA"
  
  return(out)
}
