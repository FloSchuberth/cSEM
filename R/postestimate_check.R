#' Check fit
#'
#' Check model fit using common fit measures.
#' 
#' The following fit measures are used (TODO).
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#'
#' @return (TODO)
#' @export

check <- function(.object){
  
  ## Check if cSEMResults object
  if(class(.object) != "cSEMResults") {
    stop("`.object` must be of class `cSEMResults`.", call. = FALSE)
  }
  
  paste("not yet implemented")
}

