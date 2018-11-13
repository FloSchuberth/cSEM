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

assess <- function(.object) {
  UseMethod("assess")
}

#' @describeIn assess (TODO)
#' @export

assess.cSEMResults_default <- function(.object){
  
  paste("not yet implemented")
}

#' @describeIn assess (TODO)
#' @export

assess.cSEMResults_multi <- function(.object){
  
  paste("not yet implemented")
}

#' @describeIn assess (TODO)
#' @export

assess.cSEMResults_2ndorder <- function(.object){
  
  paste("not yet implemented")
}