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
      cat2("\n\t", ansi_align(col_cyan(i), max(18, max(nchar(names(x))))))
      
      if(inherits(x[[i]], "cSEMResults_2ndorder")) {
        for(j in names(x[[i]])) {
          cat2("\n\t  ", ansi_align(names(x[[i]][j]),  15), ": " )
          cat2(ifelse(sum(status[[i]][[j]]) == 0, col_green("successful"), col_red("not successful")))
        }
      } else {
        cat2(": ", ifelse(sum(status[[i]]) == 0, col_green("successful"), col_red("not successful")))
      }

    }
    if(sum(unlist(status)) != 0) {
      cat2(
        "\n\nSee ", col_magenta("verify"), "(", col_cyan("<object-name>"), ")", " for details.")
    }
    cat2("\n\nThe result for each group/data set is a list of class ", style_bold("cSEMResults"), "",
        "\nwith list elements:\n\n\t")
    
  } else if(inherits(x, "cSEMResults_2ndorder")) {
    
    cat2("\n\nEstimation status by stage:\n")
    for(i in names(x)) {
      cat2("\n\t", ansi_align(col_cyan(i), 18), ": ",
          ifelse(sum(status[[i]]) == 0, col_green("successful"), col_red("not successful")))
    }
    if(sum(unlist(status)) != 0) {
      cat2("\n\nSee ", col_magenta("verify"), "(", col_cyan("<object-name>"), ")", 
          " for details.")
    }
    cat2("\n\nThe result for each stage is a list of class ", style_bold("cSEMResults"), "",
        "\nwith list elements:\n\n\t")
  } else {
    cat2(
      "\n\nEstimation was ", ifelse(sum(status) == 0, 
                                    col_green("successful"), col_red("not successful")), ".")
    if(sum(unlist(status)) != 0) {
      cat2(" See ", col_magenta("verify"), "(", col_cyan("<object-name>"), ")", 
          " for details.")
    }
    cat("\n\nThe result is a list of class ", style_bold("cSEMResults"), " with list elements:\n\n\t", sep = "")
  }
  if(inherits(x, "cSEMResults_multi") & inherits(x, "cSEMResults_2ndorder")) {
    cat2(
      "- ", col_green("First_stage\n\t"),
      "\t- ", col_green("Estimates\n\t"),
      "\t- ", col_green("Information\n\t"),
      "- ", col_green("Second_stage\n\t"),
      "\t- ", col_green("Estimates\n\t"),
      "\t- ", col_green("Information\n\n"))
  } else {
    cat2(
      "- ", col_green("Estimates\n\t"),
      "- ", col_green("Information\n\n"))
  }
  cat2("To get an overview or help type:\n\n\t",
      "- ", col_yellow("?"), col_cyan("cSEMResults"),"\n\t",
      "- ", col_magenta("str"), "(", col_cyan("<object-name>"), ")\n\t",
      "- ", col_magenta("listviewer"), col_yellow("::"), col_magenta("jsondedit"),
      "(", col_cyan("<object-name>"), ", ", col_red("mode"), " = ", col_cyan("'view'"), ")\n\n")
  cat2("If you wish to access the list elements directly type e.g. \n\n\t",
      "- ", col_cyan("<object-name>"), col_yellow("$"), col_green("Estimates"), "\n\n")
  cat2("Available postestimation commands:\n\n\t",
      "- ", col_magenta("assess"), "(", col_cyan("<object-name>"), ")\n\t",
      "- ", col_magenta("infer"), "(", col_cyan("<object-name"), ")\n\t",
      "- ", col_magenta("predict"), "(", col_cyan("<object-name>"), ")\n\t",
      "- ", col_magenta("summarize"), "(", col_cyan("<object-name>"), ")\n\t",
      "- ", col_magenta("verify"), "(", col_cyan("<object-name>"), ")\n")
  cat(rule2(type = 2), "\n")
}