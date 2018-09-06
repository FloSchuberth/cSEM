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
  
  ## Check if cSEMResults object
  if(class(.object) != "cSEMResults") {
    stop("`.object` must be of class `cSEMResults`.", call. = FALSE)
  }
  
  if(attr(.object, "single") == FALSE) {
    stop("`test()` not applicable to multiple groups or data sets.\n",
         "Use `lapply(.object, test)` instead.",
         call. = FALSE)
  }
  
  ## Select tests to compute
  if(.which_test %in% c("all", "testOMF")) {
    out_omf <- testOMF(.object, ...)
  } else {
    out_omf <- "Not computed"
  }
  
  if(.which_test %in% c("all", "testMGA")) {
    out_mgd <- testMGD(.object, ...)
  } else {
    out_mgd <- "Not computed"
  }
  
  out <- list(
    "OMF" = out_omf,
    "MGA" = out_mga
  )
  
  return(out)
}