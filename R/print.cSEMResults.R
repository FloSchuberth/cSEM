#' `cSEMResults` method for `print()`
#'
#' The [cSEMResults] method for the generic function [print()].
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cSEMResults]
#'
#' @export
#' @keywords internal
print.cSEMResults <- function(x, ...) {
  
  # Uses function rule2() in file zz_utils.R
  
  cat2(
    rule2(type = 2), "\n",
    rule2("Overview", type = 1, align = "center")
    )
  
  ## Verify
  status <- verify(x)
  
  if(inherits(x, "cSEMResults_multi")) {
    cat2("\n\nEstimation status by group/data set:\n")
    
    for(i in names(x)) {
      cat2("\n\t", col_align(cyan(i), max(18, max(nchar(names(x))))))
      
      if(inherits(x[[i]], "cSEMResults_2ndorder")) {
        for(j in names(x[[i]])) {
          cat2("\n\t  ", col_align(names(x[[i]][j]),  15), ": " )
          cat2(ifelse(sum(status[[i]][[j]]) == 0, green("successful"), red("not successful")))
        }
      } else {
        cat2(": ", ifelse(sum(status[[i]]) == 0, green("successful"), red("not successful")))
      }

    }
    if(sum(unlist(status)) != 0) {
      cat2(
        "\n\nSee ", magenta("verify"), "(", cyan("<object-name>"), ")", " for details.")
    }
    cat2("\n\nThe result for each group/data set is a list of class " %+% bold("cSEMResults") %+%"",
        "\nwith list elements:\n\n\t")
    
  } else if(inherits(x, "cSEMResults_2ndorder")) {
    
    cat2("\n\nEstimation status by stage:\n")
    for(i in names(x)) {
      cat2("\n\t", col_align(cyan(i), 18), ": ",
          ifelse(sum(status[[i]]) == 0, green("successful"), red("not successful")))
    }
    if(sum(unlist(status)) != 0) {
      cat2("\n\nSee ", magenta("verify"), "(", cyan("<object-name>"), ")", 
          " for details.")
    }
    cat2("\n\nThe result for each stage is a list of class " %+% bold("cSEMResults") %+%"",
        "\nwith list elements:\n\n\t")
  } else {
    cat2(
      "\n\nEstimation was ", ifelse(sum(status) == 0, 
                                    green("successful"), red("not successful")), ".")
    if(sum(unlist(status)) != 0) {
      cat2(" See ", magenta("verify"), "(", cyan("<object-name>"), ")", 
          " for details.")
    }
    cat("\n\nThe result is a list of class " %+% bold("cSEMResults") %+%" with list elements:\n\n\t")
  }
  if(inherits(x, "cSEMResults_multi") & inherits(x, "cSEMResults_2ndorder")) {
    cat2(
      "- ", green("First_stage\n\t"),
      "\t- ", green("Estimates\n\t"),
      "\t- ", green("Information\n\t"),
      "- ", green("Second_stage\n\t"),
      "\t- ", green("Estimates\n\t"),
      "\t- ", green("Information\n\n"))
  } else {
    cat2(
      "- ", green("Estimates\n\t"),
      "- ", green("Information\n\n"))
  }
  cat2("To get an overview or help type:\n\n\t",
      "- ", yellow("?"), cyan("cSEMResults"),"\n\t",
      "- ", magenta("str"), "(", cyan("<object-name>"), ")\n\t",
      "- ", magenta("listviewer"), yellow("::"), magenta("jsondedit"),
      "(", cyan("<object-name>"), ", ", red("mode"), " = ", cyan("'view'"), ")\n\n")
  cat2("If you wish to access the list elements directly type e.g. \n\n\t",
      "- ", cyan("<object-name>"), yellow("$"), green("Estimates"), "\n\n")
  cat2("Available postestimation commands:\n\n\t",
      "- ", magenta("assess"), "(", cyan("<object-name>"), ")\n\t",
      "- ", magenta("infer"), "(", cyan("<object-name"), ")\n\t",
      "- ", magenta("predict"), "(", cyan("<object-name>"), ")\n\t",
      "- ", magenta("summarize"), "(", cyan("<object-name>"), ")\n\t",
      "- ", magenta("verify"), "(", cyan("<object-name>"), ")\n")
  cat(rule2(type = 2), "\n")
}