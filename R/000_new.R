
# Function that calculates unstandardized construct scores
calculateUnstandardizedConstructScores=function(.object){
  
  dataScaled=.object$Information$Data
  
  # Scale data back
  dataUnscaled=dataScaled * attr(dataScaled, 'scaled:scale') + attr(dataScaled, 'scaled:center') 
  
  #Unstandardized scores equal X_i%*%w_i/sum(w_i)
  
  w_transformed <- t(.object$Estimates$Weight_estimates)%*%solve(diag(rowSums(.object$Estimates$Weight_estimates)))
  
  ScoresUnstandardized <- dataUnscaled%*%w_transformed
  
  colnames(ScoresUnstandardized) =  rownames(.object$Estimates$Weight_estimates)
  
  return(ScoresUnstandardized)
}