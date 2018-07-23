status=function(object,...){
  # 0: Everything is fine
  # 1 PLS-PM has not not converged
  # 2: at least one absolute standardized loading estimate is larger than 1,
  # which implies either a negative veriance of the measurement error or
  # a correlation larger than 1
  # 3: construct VCV is not positive semi-definit
  # 4 model-implied indicators VCV is not positive semi-definit
  stat=NULL
  
  # Only if PLS is used, it is necessary to check 
  if(!object$Information$Converged){
  stat=c(stat,1)    
    
  }
  if(max(abs(bject$Estimates$Cross_loadings))>1){
    stat=c(stat,2)
    }
  if(!matrixcalc::is.positive.semi.definite(object$Estimates$Construct_VCV)){
    stat = c(stat,3)
    
  }
  
  if(!matrixcalc::is.positive.semi.definite(fitted(object))){
    stat = c(stat,4)
  }
  
# if no problem occured, it seems that the estimation is fine
  if(is.null(stat)){stat=0}
  
  return(stat)
}