#' `cSEMTestOMF` method for `print()`
#'
#' The `cSEMTestOMF` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cSEMResults], [testOMF()]
#'
#' @export
#' @keywords internal
print.cSEMTestOMF <- function(x, ...) {
  
  cat2(
    rule2(type = 2), "\n",
    rule2("Test for overall model fit based on Beran & Srivastava (1985)")
  )
  
  ## Null hypothesis -----------------------------------------------------------
  if(!isTRUE(x$Information$Saturated)){
  cat2(
    "\n\nNull hypothesis:\n\n", 
    boxx(c("H0: The model-implied indicator covariance matrix equals the", 
           "population indicator covariance matrix."), float = "center",width=80)
  )
  }else if(isTRUE(x$Information$Saturated)){
    cat2(
      "\n\nNull hypothesis:\n\n", 
      boxx(c("H0: The model-implied indicator covariance matrix with a saturated", 
             "structural model equals the population indicator covariance matrix."), float = "center",width=80)
    )
  }
  
  ## Test statistic and critical value -----------------------------------------
  cat("\n\nTest statistic and critical value: \n\n\t", sep = "")
  
  cat(
    col_align("", width = 20),
    col_align("", width = 14), 
    "\t",
    col_align("Critical value", width = 8*ncol(x$Critical_value), 
              align = "center"),
    sep = ""
  )
  cat(
    "\n\t",
    col_align("Distance measure", width = 20),
    col_align("Test statistic", width = 14), 
    "\t",
    sep = ""
  )
  
  for(i in colnames(x$Critical_value)) {
    cat(col_align(i, width = 6, align = "center"), "\t", sep = "")
  }
  
  cat("\n\t")
  
  for(j in seq_along(x$Test_statistic)) {
    cat(
      col_align(names(x$Test_statistic)[j], width = 20),
      col_align(sprintf("%.4f", x$Test_statistic[j]), width = 14, 
                align = "center"), 
      "\t", 
      sep = ""
    )
    for(k in 1:ncol(x$Critical_value)) {
      cat(sprintf("%.4f", x$Critical_value[j, k]), "\t", sep = "")
    }
    cat("\n\t")
  }
  
  ## Decision ------------------------------------------------------------------
  cat("\n\nDecision: \n\n\t", sep = "")
  
  # Width of columns
  l <- apply(x$Decision, 2, function(x) {
    ifelse(any(x == TRUE), nchar("Do not reject"), nchar("reject"))
  })
  
  l1 <- max(c(sum(l) + 3*(ncol(x$Decision) - 1)), nchar("Significance level"))
  
  cat(
    col_align("", width = 20), 
    "\t",
    col_align("Significance level", 
              width = l1, 
              align = "center"),
    sep = ""
  )
  cat(
    "\n\t",
    col_align("Distance measure", width = 20), 
    "\t",
    sep = ""
  )
  
  for(i in colnames(x$Critical_value)) {
    cat(col_align(i, width = l[i], align = "center"), "\t", sep = "")
  }
  
  cat("\n\t")
  
  for(j in seq_along(x$Test_statistic)) {
    
    cat(col_align(names(x$Test_statistic)[j], width = 20), "\t", sep = "")
    
    for(k in 1:ncol(x$Critical_value)) {
      cat(
        col_align(ifelse(x$Decision[j, k], 
                         green("Do not reject"), red("reject")), 
                  width = l[k], align = "center"), 
        "\t", 
        sep = ""
      )
    }
    cat("\n\t")
  }
  
  ## Additional information ----------------------------------------------------
  cat("\nAdditional information:")
  cat2(
    "\n\n\tOut of ", x$Information$Total_runs , " bootstrap replications ", 
    x$Information$Number_admissibles, " are admissible.\n\t",
    "See ", yellow("?"), magenta("verify"), "()",
    " for what constitutes an inadmissible result.\n\n\t",
    "The seed used was: ", x$Information$Seed, "\n"
  )
  cat2(rule2(type = 2), "\n")
}