#' `cSEMTestHausman` method for `print()`
#'
#' The `cSEMTestHausman` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cSEMResults], [testHausman()]
#'
#' @export
#' @keywords internal
print.cSEMTestHausman <- function(x, ...) {
  
  x1 <- x$Regressions
  x2 <- x$Information
  
  cat2(
    rule2(type = 2), "\n",
    rule(center = "Regression-based Hausman test",
         width = 80)
  )
  
  ## Null hypothesis -----------------------------------------------------------
  cat2(
    "\n\nNull hypothesis:\n\n", 
    boxx(c("H0: Variable(s) suspected to be endogenous are uncorrelated with the",
           "error term (no endogeneity)."), float = "center",width=80)
  )
  
  ## Regression output ---------------------------------------------------------
  cat("\n\nRegression output: \n\n\t", sep = "")
  
  for(i in seq_along(x1)) {
    cat2("\n  Dependent construct: '", names(x1)[i], "'\n")
    l <- max(nchar(x1[[i]]$Name))
    cat2(
      "\n\t", 
      col_align("Independent construct", max(l, nchar("Independent construct")) + 2), 
      col_align("Estimate", 10, align = "right"), 
      col_align("Std. error", 12, align = "right"),
      col_align("t-stat.", 10, align = "right"), 
      col_align("p-value", 10, align = "right")
    )
    for(j in 1:nrow(x1[[i]])) {
      cat2(
        "\n\t", 
        col_align(x1[[i]][j, "Name"], max(l, nchar("Independent construct")) + 2), 
        col_align(sprintf("%.4f", x1[[i]][j, "Estimate"]), 10, align = "right"),
        col_align(sprintf("%.4f", x1[[i]][j, "Std_error"]), 12, align = "right"),
        col_align(sprintf("%.4f", x1[[i]][j, "t_stat"]), 10, align = "right"),
        col_align(sprintf("%.4f", x1[[i]][j, "p_value"]), 10, align = "right")
      )
    }
    cat2("\n")
  }
  cat2(rule2(type = 2), "\n")
}