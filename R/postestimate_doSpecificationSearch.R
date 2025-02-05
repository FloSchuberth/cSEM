#' Specification search 
#' 
#' \lifecycle{experimental}
#'
#' The `doSpecificationSearch()` 
#' 
#' @usage doSpecificationSearch(
#' .object = NULL) 
#'
#' @inheritParams csem_arguments
#' 
#' 
#' @return List
#' 
#' @seealso [csem()], [cSEMResults], \link[DiagrammeR]{grViz}
#' 
#' 
#' @example inst/examples/example_doSpecificationSearch.R
#' 
#' @export


doSpecificationSearch <- function(.object = NULL){
  
  out=list()
  
  
  
  out[["Estimates"]]=.object$Estimates
  out[[2]]=.object$Information$Model
  out[[3]]=.object$Information$Data
  out[[4]]=.object$Information$Arguments$.data
  
  # ga(X)

  out
  }