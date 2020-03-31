#' Do a importance performance analysis
#'
#' Description
#' 
#' 
#' @usage doIPA(
#'  .object             = NULL
#'  )
#'
#' @inheritParams csem_arguments
#'
#' @return A list of class `cSEMIPA` with a corresponding method for `plot()`. 
#'   See: [plot.cSEMIPA()].
#' 
#' @seealso [csem()], [cSEMResults], [plot.cSEMIPA()]
#' @export

doIPA = function(.object){
  # Check whether Pearson was applied 
  if(.object$Information$Type_of_indicator_correlation != 'Pearson'){
    cSEM:::stop2("Importance-performance matrix can currently only be obtained if Pearson correlation is applied.")
  }
  
  # Check whether all constructs are composite
  if(!all(.object$Information$Model$construct_type %in% "Composite")){
    cSEM:::stop2("Importance-performance matrix can only be obtained if all constructs are composites.")
  }
  
  # Für multi verbieten
  
  # Redo the estimation with scaled indicators, see Eq. 13.1 in Henseler (2020).
  # collect arguments
  arguments <- .object$Information$Arguments
  
  data_org <- arguments$.data
  
  # scale the indicators
  data_scaled <- apply(data_org,2,function(x){
    (x-min(x))*100/(max(x)-min(x))
  })
  
  # Overwrite the original dataset
  arguments$.data <- data_scaled
  
  est_new<-do.call(csem,arguments)
  
  
  # get unstandardized construct scores
  scores <- getConstructScores(.object = est_new,.standardized = FALSE)
  scale(scores)
  
  # Calculate importance values, see Eq. 13.4 in Henseler (2020)
  effects <- cSEM:::calculateEffects(.object = est_new,.output_type = 'matrix')
  total_effects_std <- effects$Total_effect
  
  # calculate unstandardized total effects
  total_effects_unstd <- total_effects_std
  for(i in rownames(total_effects_unstd)){
    for(j in colnames(total_effects_unstd)){
     total_effects_unstd[i,j] <- total_effects_std[i,j]* sd(scores[,j])/sd(scores[,i])
    }
  }
  
  # Calculate the performance values
  # In principle, we can allow for other location measures
  perf<-colMeans(scores)
  
  # collect output
  out<-list(Importance=total_effects_unstd,
            Performance=perf)
  
  class(out) <- "cSEMIPA"
  
  return(out)
}
