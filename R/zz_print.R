#' `cSEMResults` method for `print()`
#'
#' The [cSEMResults] method for the generic function [print()].
#'
#' @usage print(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#'
#' @export
#' 
print.cSEMResults <- function(.object) {
    
  cat(
    rule(line = "bar2"), "\n",
    rule(center = "Overview"), 
    sep = "")
  
  if(attr(.object, "single") == FALSE) {
    cat("\n\nEstimation status by group/data set:\n", sep = "")
    for(i in names(.object)) {
      cat("\n\t", col_align(cyan(i), 15), ": ",
          ifelse(.object[[i]]$Information$Weight_info$Convergence_status == TRUE, 
                 green("successful"), red("not successful")), sep = "")
    }
    cat("\n\nThe result for each group/data set is a list of class " %+% bold("cSEMResults") %+%"", 
        "\nwith list elements:\n\n\t", sep = "")
    
  } else {
    cat(
      "\n\nEstimation was ", ifelse(.object$Information$Weight_info$Convergence_status == TRUE, 
                                    green("successful."), red("not successful.")), sep = "") 
    cat("\n\nThe result is a list of class " %+% bold("cSEMResults") %+%" with list elements:\n\n\t")
  }
  cat(
      "- ", green("Estimates\n\t"),
      "- ", green("Information\n\n"), sep = "")
  cat("To get an overview or help type:\n\n\t",
      "- ", yellow("?"), cyan("cSEMResults"),"\n\t",
      "- ", magenta("str"), "(", cyan("<object-name>"), ")\n\t",
      "- ", magenta("listviewer"), yellow("::"), magenta("jsondedit"),
      "(", cyan("<object-name>"), ", ", red("mode"), " = ", cyan("'view'"), ")\n\n", sep = "")
  cat("If you wish to access the list elements directly type e.g. \n\n\t",
      "- ", cyan("<object-name>"), yellow("$"), green("Estimates"), "\n\n", sep = "")
  cat("Available postestimation commands:\n\n\t",
      "- ", magenta("check"), "(", cyan("<object-name>"), ")\n\t",
      "- ", magenta("fit"), "(", cyan("<object-name>"), ")\n\t",
      "- ", magenta("summarize"), "(", cyan("<object-name>"), ")\n\t"  ,
      "- ", magenta("test"), "(", cyan("<object-name>"), ")\n\t"  ,
      "- ", magenta("verify"), "(", cyan("<object-name>"), ")\n", sep = "")
  cat(rule(line = "bar2"), "\n")
}

#' `cSEMSummary` method for `print()`
#'
#' The [cSEMSummary] method for the generic function [print()]. 
#'
#' @usage print(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults], [summarize()]
#'
#' @export
#'

print.cSEMSummarize <- function(.object) {
  
  x1 <- .object$Estimates
  x2 <- .object$Information
  
  cat(
    rule(line = "bar2"), "\n",
    rule(center = "Overview"), 
    "\n", sep = "")
  
  cat(
    col_align("\n\tNumber of Observations", 25), "= ", nrow(x2$Arguments$.data),
    col_align("\n\tWeight estimator", 25), "= ", x2$Arguments$.approach_weights, 
    sep = "")
  
  if(x2$Arguments$.approach_weights == "PLS") {
    cat(
      col_align("\n\tInner Weighting Scheme", 25), "= ", 
      x2$Arguments$.PLS_weight_scheme_inner, 
      sep = "")
  }

  cat(
    col_align("\n\tPath estimator", 25), "= ", x2$Arguments$.approach_paths,
    col_align("\n\tType of path model", 25), "= ", x2$Model$model_type,
    col_align("\n\tDisattenuated", 25), "= ", ifelse(x2$Arguments$.disattenuate, "Yes", "no"),
    sep = "")
  
  cat("\n\n\tConstruct modeled as:\n\t","---------------------", sep = "")
  l <- max(nchar(names(x2$Model$construct_type)))

  for(i in names(x2$Model$construct_type)) {
    cat("\n\t", col_align(i, l + 2), ": ", x2$Model$construct_type[i], sep = "")
  }

  if(x2$Arguments$.approach_weights == "PLS") {
    cat("\n\n\tPLS Modes:\n\t","----------", sep = "")

    for(i in names(x2$Weight_info$Modes)) {
      cat("\n\t", col_align(i, l + 2), ": ", x2$Weight_info$Modes[i], sep = "")
    }
  }

  cat("\n\n", rule(center = "Estimates"), "\n\n", sep = "")

  cat("Estimated Path Coefficients:\n============================", sep = "")
  l <- max(nchar(x1$Path_estimates[, 1]))
  
  cat("\n\t", col_align("", l + 2), col_align("Estimate"))
  for(i in 1:nrow(x1$Path_estimates)) {
    cat("\n\t", col_align(x1$Path_estimates[i, 1], l + 2), 
        col_align(sprintf("%.4f", x1$Path_estimates[i, 2]), 10, align = "right"), 
        sep = "")
  }
  
  cat("\n\nEstimated Loadings:\n===================", sep = "")
  l <- max(nchar(x1$Loading_estimates[, 1]))
  
  cat("\n\t", col_align("", l + 2), col_align("Estimate"))
  for(i in 1:nrow(x1$Loading_estimates)) {
    cat("\n\t", col_align(x1$Loading_estimates[i, 1], l + 2), 
        col_align(sprintf("%.4f", x1$Loading_estimates[i, 2]), 10, align = "right"), 
        sep = "")
  }


  cat("\n\nEstimated Weights:\n==================\n", sep = "")
  l <- max(nchar(x1$Weight_estimates[, 1]))
  
  cat("\n\t", col_align("", l + 2), col_align("Estimate"))
  for(i in 1:nrow(x1$Weight_estimates)) {
    cat("\n\t", col_align(x1$Weight_estimates[i, 1], l + 2), 
        col_align(sprintf("%.4f", x1$Weight_estimates[i, 2]), 10, align = "right"), 
        sep = "")
  }

  # if(x$Weight_estimator == "PLS") {
  #   cat("\nEstimated Correction Factors:\n=============================\n", sep = "")
  #   print(x$Correction_factors)
  # }

  # cat("\n\n", cli::rule(center = "Other output"), "\n\n\t", sep = "")
  # 
  # cat("<not yet implemented>")
  # 
  # cat("\n\n", cli::rule(center = "Fit Indices"), "\n\n\t", sep = "")
  # 
  # cat(crayon::col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\t",
  #     crayon::col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\t",
  #     crayon::col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\n",
  #     sep = "")
  
  cat("\n", rule(line = "bar2"), sep = "")
}

#' `cSEMVerify` method for `print()`
#'
#' The `cSEMVerify` method for the generic function [print()]. 
#'
#' @usage print(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#'
#' @export
#'
print.cSEMVerify <- function(.object) {
  
  cat(rule(line = "bar2"), sep = "")
  
  cat("\n\nVerify admissibility of the model estimates:\n", sep = "")
  
  if(sum(.object) == 0) {
    cat(green("\n\t admissible\n"), sep = "")
  } else {
    cat(red("\n\t inadmissible\n"), sep = "")
  }

  text <- c("1" = "Convergence", 
            "2" = "At least one standardized loadings > 1", 
            "3" = "Construct VCV not positive semi-definite", 
            "4" = "Model-implied VCV not positive semi-definite",
            "5" = "At least one construct reliability > 1")
  cat("\nStatus code, status, and description:\n\n", sep = "")
  for(i in names(.object)) {
    cat("\t", i, ": ", col_align(ifelse(.object[i] == FALSE, green("ok"), red("not ok")), 8), text[i], "\n", sep = "")
  }
  cat(rule(line = "bar2"), sep = "")
}

#' `cSEMTestOMF` method for `print()`
#'
#' The `cSEMTestOMF` method for the generic function [print()]. 
#'
#' @usage print(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#'
#' @export
#'
print.cSEMTestOMF <- function(.object) {
  
  cat(
    rule(line = "bar2"), "\n",
    rule(center = "Test for Overall Model Fit"), 
    sep = "")
  cat("\n\nNull Hypothesis:\n\n", 
      boxx(c("H0: No significant difference between empirical and", 
       "model-implied indicator covariance matrix."), padding = 1, float = "center"), sep = "")
  
  cat("\n\nTest statistic and critical value: \n\n\t", sep = "")
  
  cat(col_align("", width = 20, align = "left"),
      col_align("", 18, align = "left"), "\t",
      col_align("Critical value", 8*nrow(.object$Critical_value), align = "center"), "\n\t",
      col_align("Distance measure", width = 20, align = "left"),
      col_align("Test statistic", 18, align = "left"), "\t",
      sep = "")
  
  for(i in rownames(.object$Critical_value)) {
    cat(col_align(i, 6, align = "center"), "\t", sep = "")
  }
  cat("\n\t")
  for(j in seq_along(.object$Test_statistic)) {
    cat(col_align(names(.object$Test_statistic)[j], 20),
        col_align(sprintf("%.4f", .object$Test_statistic[j]), 18), "\t", sep = "")
    for(k in 1:nrow(.object$Critical_value)) {
      cat(sprintf("%.4f", .object$Critical_value[k, j]), "\t",
          sep = "")
    }
    cat("\n\t")
  }
  
  cat("\n\nDecision: \n\n\t", sep = "")
  
  cat(col_align("", width = 20, align = "left"), "\t",
      col_align("Significance level", 8*nrow(.object$Critical_value), align = "center"), "\n\t",
      col_align("Distance measure", width = 20, align = "left"), "\t",
      sep = "")
  
  for(i in rownames(.object$Critical_value)) {
    cat(col_align(i, 8, align = "center"), "\t", sep = "")
  }
  cat("\n\t")
  for(j in seq_along(.object$Test_statistic)) {
    cat(col_align(names(.object$Test_statistic)[j], 20), "\t", sep = "")
    for(k in 1:nrow(.object$Critical_value)) {
      cat(col_align(ifelse(.object$Decision[k, j], green("Do not reject"), red("reject")), 8), "\t", sep = "")
    }
    cat("\n\t")
  }

  cat("\nAdditonal information:")
  cat("\n\n\tThere are ", .object$Number_admissibles, " admissibles bootstrap results.\n\t",
      "See ", yellow("?"), magenta("verify"), "()",
      " for what constitutes an inadmissible result.\n", sep = "")
  cat(rule(line = "bar2"), sep = "")
}

#' `cSEMTestMGD` method for `print()`
#'
#' The `cSEMTestMGD` method for the generic function [print()]. 
#'
#' @usage print(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#'
#' @export
#'
print.cSEMTestMGD <- function(.object) {
  
  cat(
    rule(line = "bar2"), "\n",
    rule(center = "Test for Multigroup Differences"), 
    sep = "")
  cat("\n\nNull Hypothesis:\n\n", 
      boxx("H0: No significant difference between groups.", 
           padding = 1, float = "center"), sep = "")
  
  cat("\n\nTest statistic and critical value: \n\n\t", sep = "")
  
  cat(col_align("", width = 20, align = "left"),
      col_align("", 18, align = "left"), "\t",
      col_align("Critical value", 8*nrow(.object$Critical_value), align = "center"), "\n\t",
      col_align("Distance measure", width = 20, align = "left"),
      col_align("Test statistic", 18, align = "left"), "\t",
      sep = "")
  
  for(i in rownames(.object$Critical_value)) {
    cat(col_align(i, 6, align = "center"), "\t", sep = "")
  }
  cat("\n\t")
  for(j in seq_along(.object$Test_statistic)) {
    cat(col_align(names(.object$Test_statistic)[j], 20),
        col_align(sprintf("%.4f", .object$Test_statistic[j]), 18), "\t", sep = "")
    for(k in 1:nrow(.object$Critical_value)) {
      cat(sprintf("%.4f", .object$Critical_value[k, j]), "\t",
          sep = "")
    }
    cat("\n\t")
  }
  
  cat("\n\nDecision: \n\n\t", sep = "")
  
  cat(col_align("", width = 20, align = "left"), "\t",
      col_align("Significance level", 8*nrow(.object$Critical_value), align = "center"), "\n\t",
      col_align("Distance measure", width = 20, align = "left"), "\t",
      sep = "")
  
  for(i in rownames(.object$Critical_value)) {
    cat(col_align(i, 8, align = "center"), "\t", sep = "")
  }
  cat("\n\t")
  for(j in seq_along(.object$Test_statistic)) {
    cat(col_align(names(.object$Test_statistic)[j], 20), "\t", sep = "")
    for(k in 1:nrow(.object$Critical_value)) {
      cat(col_align(ifelse(.object$Decision[k, j], green("Do not reject"), red("reject")), 8), "\t", sep = "")
    }
    cat("\n\t")
  }
  
  cat("\nAdditonal information:")
  cat("\n\n\tThere are ", .object$Number_admissibles, " admissibles bootstrap results.\n\t",
      "See ", yellow("?"), magenta("verify"), "()",
      " for what constitutes an inadmissible result.\n", sep = "")
  cat(rule(line = "bar2"), sep = "")
}