#' Internal: get significance stars
#'
#' Transforms a p-value into stars. 
#' 
#' `.pvalue` Numeric. A p-value that is transformed into a star.
#' 
#' 
#' @usage get_significance_stars(
#' .pvalue
#' )
#' 
#' @inheritParams csem_arguments
#' 
#' @return Character string. A p-value transformed into a star.
#' 
#' @keywords internal
#'
get_significance_stars <- function(.pvalue) {
  if (is.na(.pvalue)) return("")
  else if (.pvalue < 0.001) return("***")
  else if (.pvalue < 0.01) return("**")
  else if (.pvalue < 0.05) return("*")
  else return("")
}
