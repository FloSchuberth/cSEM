#' Function to check the status of a cSEMResults object
#'
#' It returns a vector where the number indicates the following:
#' \itemize{
#' \item 0: Everything is fine
#' \item 1 Algorithm has not not converged
#' \item 2: at least one absolute standardized loading estimate is larger than 1,
#' which implies either a negative veriance of the measurement error or
#' a correlation larger than 1
#' \item 3: construct VCV is not positive semi-definit
#' \item 4 model-implied indicators VCV is not positive semi-definit
#' }
#' @usage status(object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem], [cSEMResults]
#'
#' @export
#'

status <- function(object,...){
  # 0: Everything is fine
  # 1 Algorithm has not not converged
  # 2: at least one absolute standardized loading estimate is larger than 1,
  # which implies either a negative veriance of the measurement error or
  # a correlation larger than 1
  # 3: construct VCV is not positive semi-definit
  # 4 model-implied indicators VCV is not positive semi-definit
  stat=NULL
  
  # Only if PLS is used, it is necessary to check 
  
  if(!object$Information$Convergence_status){
    
    stat=c(stat,1)    
    
  }
  if(max(abs(object$Estimates$Cross_loadings))>1){
    stat=c(stat,2)
  }
  if(!matrixcalc::is.positive.semi.definite(object$Estimates$Construct_VCV)){
    stat = c(stat,3)
    
  }
  
  # if(!matrixcalc::is.positive.semi.definite(fitted(object))){
  #   stat = c(stat,4)
  # }
  
  # if no problem occured, it seems that the estimation is fine
  if(is.null(stat)){stat=0}
  
  return(stat)
}