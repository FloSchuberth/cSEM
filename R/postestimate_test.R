#' Diverse tests
#'
#' Test estimated model (TODO).
#' 
#' More details here (TODO).
#'
#' @inheritParams csem_arguments
#' @param ... Further arguments.
#'
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#'
#' @return TODO
#'
#' @export

test <- function(.object, .which_test = "all", ...) {
  
  ## Select tests to compute
  if(.which_test %in% c("all", "testOMF")) {
    out_omf <- testOMF(.object, ...)
  } else {
    out_omf <- "Not computed"
  }
  
  if(.which_test %in% c("all", "testMGA")) {
    out_mga <- testMGA(.object, ...)
  } else {
    out_mga <- "Not computed"
  }
  
  out <- list(
    "OMF" = out_omf,
    "MGA" = out_mga
  )
  
  return(out)
}