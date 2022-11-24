#' `cSEMTestCVPAT` method for `print()`
#'
#' The `cSEMTestCVAT` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cSEMResults], [testCVPAT()]
#'
#' @export
#' @keywords internal
print.cSEMTestCVPAT <- function(x, ...) {
  
  cat2(
    rule2(type = 2), "\n",
    rule2("Cross validated predictive ability test (Liengaard et al., 2020)")
  )
  
  ## Null hypothesis -----------------------------------------------------------
    cat2(
      "\n\nNull hypothesis:\n\n", 
      boxx(c("H0: Both models perform equally in predicting indicators of endogenous constructs."), float = "center",width=80)
    )
  
  ## Test statistic and p-value-------------------------------------------------
  cat("\n\nTest statistic and p-value: \n\t", sep = "")

  cat2(
  col_align("\n\tTest statistic", 35), "= ", x$test_statistic,
  col_align("\n\tP value", 35), "= ", x$p_value,
  col_align("\n\tDegrees of freedom", 35), "= ", x$Information$Degrees_of_Freedom,
  col_align("\n\tTesttype", 35), "= ", x$Testtype
  )
  
  cat("\n\nInformation:\n\t", sep = "")
  
  cat2(
    col_align("\n\tSeed", 35), ": ", x$Information$Seed,
    col_align("\n\tSample Size", 35), ": ", x$Information$`Sample Size`
  )
  
}