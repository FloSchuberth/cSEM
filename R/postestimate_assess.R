#' Assess model
#'
#' Assess model using common evaluation criteria and fit measures.
#' 
#' The following evaulation criteria and fit measures are used (TODO).
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#'
#' @return (TODO)
#' @export

assess <- function(.object){
  
  ## Check if cSEMResults object
  if(class(.object) != "cSEMResults") {
    stop("`.object` must be of class `cSEMResults`.", call. = FALSE)
  }
  
  paste("not yet implemented")
}

