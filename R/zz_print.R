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
  
  cat(rule(line = "bar2"), "\n")
  cat(rule(center = "Overview"), "\n\n")
  
  if(!all(names(.object) %in% c("Estimates", "Information"))) {
    cat("Estimation status by group/data set:\n\n\t")
    for(i in names(.object)) {
      cat(col_align(cyan(i), 15), ": ",
          ifelse(.object[[i]]$Information$Weight_info$Convergence_status == TRUE, 
                 green("successful\n\t"), red("not successful\n\t")), sep = "")
    }
  } else {
    cat("Estimation was ", ifelse(.object$Information$Weight_info$Convergence_status == TRUE, 
                                  green("successful.\n"), red("not successful.\n")), sep = "") 
  }
  cat("\nThe result is a list of class " %+% bold("cSEMResults") %+%" with list elements:\n\n\t",
      "- ", green("Estimates\n\t"),
      "- ", green("Information\n\n"), sep = "")
  cat("To get an overview or help type:\n\n\t",
      "- ", magenta("?"), cyan("cSEMResults"),"\n\t",
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
print.cSEMSummary <- function(.object) {
  cat(cli::rule(line = "bar2"), "\n",
      cli::rule(center = "General Information"), "\n\n\t", sep = "")
  
  cat(crayon::col_align("Number of Observations", 25), "= ", x$Number_of_observations, "\n\t",
      crayon::col_align("Weight estimator", 25), "= ", x$Weight_estimator, "\n\t", sep = "")
  if(x$Weight_estimator == "PLS") {
    cat(crayon::col_align("Inner Weighting Scheme ", 25), "= ", x$PLS_weight_scheme_inner, "\n\t", sep = "")
  }
  cat(
    crayon::col_align("Path estimator", 25), "= ", x$Path_estimator, "\n\t",
    crayon::col_align("Convergence Status", 25), "= ", c("not yet implemented"), "\n\t",
    crayon::col_align("Overall model fit", 25), "= ", c("not yet implemented"), "\n\t",
    crayon::col_align("Degrees of Freedom", 25), "= ", c("not yet implemented"), "\n\t",
    crayon::col_align("Computation Time", 25), "= ", c("not yet implemented"), "\n\n\t",
    sep = "")
  
  cat("Construct Types:\n\t","----------------","\n\t", sep = "")
  
  for(i in seq_along(x$Construct_types$Name)) {
    cat(crayon::col_align(x$Construct_types$Name[i], 10), ": ", x$Construct_types$Type[i],"\n\t", sep = "")
  }
  cat("\n")
  
  if(x$Weight_estimator == "PLS") {
    cat("\tPLS Modes:\n\t","----------------","\n\t", sep = "")
    
    for(i in seq_along(x$Construct_type$Name)) {
      cat(crayon::col_align(x$Construct_types$Name[i], 10), ": ", x$PLS_modes[i],"\n\t", sep = "")
    }
    cat("\n")
  }
  
  cat(cli::rule(center = "Estimates"), "\n\n", sep = "")
  
  cat("Estimated Path Coefficients:\n============================\n", sep = "")
  print(x$Path_estimates, row.names = FALSE)
  
  
  cat("\nEstimated Loadings:\n===================\n", sep = "")
  print(x$Loading_estimates, row.names = FALSE)
  
  cat("\nEstimated Weights:\n==================\n", sep = "")
  print(x$Weight_estimates, row.names = FALSE)
  
  if(x$Weight_estimator == "PLS") {
    cat("\nEstimated Correction Factors:\n=============================\n", sep = "")
    print(x$Correction_factors)
  }
  
  cat("\n\n", cli::rule(center = "Other output"), "\n\n\t", sep = "")
  
  cat("<not yet implemented>")
  
  cat("\n\n", cli::rule(center = "Fit Indices"), "\n\n\t", sep = "")
  
  cat(crayon::col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\t",
      crayon::col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\t",
      crayon::col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\n",
      sep = "")
  
  cat(cli::rule(line = "bar2"))
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
  
  cat("Verify admissibility of the model estimates:\n", sep = "")
  
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
  cat("Null Hypothesis:\n\n", 
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
      cat(col_align(ifelse(.object$Decision[k, j], "Do not reject", "reject"), 8), "\t", sep = "")
    }
    cat("\n\t")
  }

  cat("\nAdditonal information:")
  cat("\n\n\tThere are ", .object$Number_admissibles, " admisibles bootstrap results.\n\tSee `?verfiy`",
      " for what constitutes an inadmissible result.", sep = "")
}
