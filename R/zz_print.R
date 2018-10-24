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
  
  if(class(.object)[2] == "cSEMResults_multi") {
    cat("\n\nEstimation status by group/data set:\n", sep = "")
    for(i in names(.object)) {
      cat("\n\t", col_align(cyan(i), 15), ": ",
          ifelse(sum(verify.cSEMResults_default(.object[[i]])) == 0, 
                 green("successful"), red("not successful")), ".", sep = "")
    }
    if(sum(verify.cSEMResults_default(.object[[i]])) != 0) {
      cat("\n\nSee ", magenta("verify"), "(", cyan("<object-name>"), ")", 
          " for details.", sep = "")
    }
    cat("\n\nThe result for each group/data set is a list of class " %+% bold("cSEMResults") %+%"",
        "\nwith list elements:\n\n\t", sep = "")
    
  } else if(class(.object)[2] == "cSEMResults_2ndorder") {
    x <- sapply(verify.cSEMResults_2ndorder(.object), sum)
    cat("\n\nEstimation status by stage:\n", sep = "")
    for(i in names(.object)) {
      cat("\n\t", col_align(cyan(i), 15), ": ",
          ifelse(sum(x) == 0, green("successful"), red("not successful")), ".", sep = "")
    }
    if(sum(x) != 0) {
      cat("\n\nSee ", magenta("verify"), "(", cyan("<object-name>"), ")", 
          " for details.", sep = "")
    }
    cat("\n\nThe result for each stage set is a list of class " %+% bold("cSEMResults") %+%"",
        "\nwith list elements:\n\n\t", sep = "")
  } else {
    cat(
      "\n\nEstimation was ", ifelse(sum(verify.cSEMResults_default(.object)) == 0, 
                                    green("successful"), red("not successful")), ".", sep = "")
    if(sum(verify.cSEMResults_default(.object)) != 0) {
      cat(" See ", magenta("verify"), "(", cyan("<object-name>"), ")", 
          " for details.", sep = "")
    }
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
      "- ", magenta("assess"), "(", cyan("<object-name>"), ")\n\t",
      "- ", magenta("fit"), "(", cyan("<object-name>"), ")\n\t",
      "- ", magenta("summarize"), "(", cyan("<object-name>"), ")\n\t"  ,
      "- ", magenta("test"), "(", cyan("<object-name>"), ")\n\t"  ,
      "- ", magenta("verify"), "(", cyan("<object-name>"), ")\n", sep = "")
  cat(rule(line = "bar2"), "\n")
}

#' `cSEMSummarize_default` method for `print()`
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
print.cSEMSummarize_default <- function(.object) {
  
  x1 <- .object$Estimates
  x2 <- .object$Information
  
  cat(
    rule(line = "bar2"), "\n",
    rule(center = "Overview"), 
    "\n", sep = "")
  
  ### Overview -----------------------------------------------------------------
  cat(
    col_align("\n\tNumber of observations", 25), "= ", nrow(x2$Arguments$.data),
    col_align("\n\tWeight estimator", 25), "= ", x2$Arguments$.approach_weights, 
    sep = "")
  
  if(x2$Arguments$.approach_weights == "PLS-PM") {
    cat(
      col_align("\n\tInner weighting scheme", 25), "= ", 
      x2$Arguments$.PLS_weight_scheme_inner, 
      sep = "")
  }
  cat(
    col_align("\n\tType of correlation", 25), "= ", 
    ifelse(x2$Arguments$.approach_cor_robust == "none", 
           "Bravais-Pearson", x2$Arguments$.approach_cor_robust),
    sep = "")
  cat(
    col_align("\n\tPath model estimator", 25), "= ", x2$Arguments$.approach_paths,
    col_align("\n\tType of path model", 25), "= ", x2$Model$model_type,
    col_align("\n\tDisattenuated", 25), "= ", 
    ifelse(x2$Arguments$.disattenuate, 
           ifelse(x2$Arguments$.approach_weights == "PLS-PM", "Yes (PLSc)", "Yes"), "no"),
    sep = "")
  
  cat("\n\n\tConstruct details:\n\t","------------------", sep = "")
  l <- max(nchar(names(x2$Model$construct_type)))
  
  cat("\n\t", 
      col_align("Name", max(l, nchar("Name")) + 2), 
      col_align("Modeled as", 13 + 2),
      col_align("Order", 12 + 2), sep = "")
  if(x2$Arguments$.approach_weights == "PLS-PM") {
    cat(col_align("Mode", 5), sep = "")
  }
  cat("\n")
  
  for(i in names(x2$Model$construct_type)) {
    cat("\n\t", 
        col_align(i, max(l, nchar("Name")) + 2), 
        col_align(x2$Model$construct_type[i], 13 + 2), 
        col_align(x2$Model$construct_order[i], 12 + 2), sep = "")
    if(x2$Arguments$.approach_weights == "PLS-PM") {
       cat(col_align(x2$Weight_info$Modes[i], 5), sep = "")
    }
  }

  ### Estimates ----------------------------------------------------------------
  cat("\n\n", rule(center = "Estimates"), "\n\n", sep = "")

  ## Path estimates
  cat("Estimated Path Coefficients:\n============================", sep = "")
  l <- max(nchar(x1$Path_estimates[, "Name"]))
  
  cat("\n\t", 
      col_align("Path", max(l, nchar("Path")) + 2), 
      col_align("Estimate", 10, align = "right"), 
      col_align("Std. Error", 12, align = "right"),
      col_align("t-value", 10, align = "right"), sep = "")
  
  for(i in 1:nrow(x1$Path_estimates)) {
    cat("\n\t", 
        col_align(x1$Path_estimates[i, "Name"], max(l, nchar("Path")) + 2), 
        col_align(sprintf("%.4f", x1$Path_estimates[i, "Estimate"]), 10, align = "right"),
        col_align(NA, 12, align = "right"),
        col_align(NA, 10, align = "right"),
        sep = "")
  }
  
  ## Loadings
  cat("\n\nEstimated Loadings:\n===================", sep = "")
  l <- max(nchar(x1$Loading_estimates[, "Name"]))
  
  cat("\n\t", 
      col_align("Loading", max(l, nchar("Loading")) + 2), 
      col_align("Estimate", 10, align = "right"), 
      col_align("Std. Error", 12, align = "right"),
      col_align("t-value", 10, align = "right"), sep = "")
  
  for(i in 1:nrow(x1$Loading_estimates)) {
    cat("\n\t", 
        col_align(x1$Loading_estimates[i, "Name"], max(l, nchar("Loading")) + 2), 
        col_align(sprintf("%.4f", x1$Loading_estimates[i, "Estimate"]), 10, align = "right"), 
        col_align(NA, 12, align = "right"),
        col_align(NA, 10, align = "right"),
        sep = "")
  }

  ## Only print weights, for constructs modeled as common factors
  temp_w <- x1$Weight_estimates[x1$Weight_estimates$Construct_type == "Common factor", , drop = FALSE] 
  
  cat("\n\nEstimated Weights:\n==================\n", sep = "")
  l <- max(nchar(temp_w[, "Name"]))
  
  cat("\n\t", 
      col_align("Weight", max(l, nchar("Weight")) + 2), 
      col_align("Estimate", 10, align = "right"), sep = "")
  for(i in 1:nrow(temp_w)) {
    cat("\n\t", 
        col_align(temp_w[i, "Name"], max(l, nchar("Weight")) + 2), 
        col_align(sprintf("%.4f", temp_w[i, "Estimate"]), 10, align = "right"), 
        sep = "")
  }

  # if(x$Weight_estimator == "PLS-PM") {
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

#' `cSEMSummarize_2ndorder` method for `print()`
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
print.cSEMSummarize_2ndorder <- function(.object) {
  
  ### Exctract name and objects 
  x11 <- .object$First_stage$Estimates
  x12 <- .object$First_stage$Information
  
  x21 <- .object$Second_stage$Estimates
  x22 <- .object$Second_stage$Information
  
  ## Collect all necessary sets
  # All constructs used in the first step (= all first order constructs)
  c_linear_1step        <- rownames(x12$Model$structural)
  # All second order constructs
  c_2nd_order           <- grep("_temp", rownames(x22$Model$structural), 
                                value = TRUE, invert = TRUE)
  # # All linear constructs of the original model
  c_linear_original     <- c(c_linear_1step, c_2nd_order)
  
  cat(
    rule(line = "bar2"), "\n",
    rule(center = "Overview"), 
    "\n", sep = "")
  
  ### Overview -----------------------------------------------------------------
  cat(
    col_align("\n\tNumber of observations", 25), "= ", nrow(x12$Arguments$.data),
    col_align("\n\tWeight estimator", 25), "= ", x12$Arguments$.approach_weights, 
    sep = "")
  
  if(x12$Arguments$.approach_weights == "PLS-PM") {
    cat(
      col_align("\n\tInner weighting scheme", 25), "= ", 
      x12$Arguments$.PLS_weight_scheme_inner, 
      sep = "")
  }
  cat(
    col_align("\n\tType of correlation", 25), "= ", 
    ifelse(x12$Arguments$.approach_cor_robust == "none", 
           "Bravais-Pearson", x12$Arguments$.approach_cor_robust),
    sep = "")
  cat(
    col_align("\n\tPath model estimator", 25), "= ", x12$Arguments$.approach_paths,
    # Its important to read the model type of the second stage (since the first
    # is always linear)
    col_align("\n\tType of path model", 25), "= ", x22$Model$model_type,
    col_align("\n\tDisattenuated", 25), "= ", 
    ifelse(x12$Arguments$.disattenuate, 
           ifelse(x12$Arguments$.approach_weights == "PLS-PM", "Yes (PLSc)", "Yes"), "no"),
    sep = "")
  
  cat("\n\n\tConstruct details:\n\t","------------------", sep = "")
  # First order constructs
  l <- max(nchar(c_linear_original))
  # Second order constructs
  
  cat("\n\t", 
      col_align("Name", max(l, nchar("Name")) + 2), 
      col_align("Modeled as", 13 + 2),
      col_align("Order", 12 + 2), sep = "")
  if(x12$Arguments$.approach_weights == "PLS-PM") {
    cat(col_align("Mode", 5), sep = "")
  }
  cat("\n")
  # First stage
  for(i in names(x12$Model$construct_type)) {
    cat("\n\t", 
        col_align(i, max(l, nchar("Name")) + 2), 
        col_align(x12$Model$construct_type[i], 13 + 2), 
        col_align("First order", 12 + 2), sep = "")
    if(x12$Arguments$.approach_weights == "PLS-PM") {
      cat(col_align(x12$Weight_info$Modes[i], 5), sep = "")
    }
  }
  # Second stage
  for(i in names(x22$Model$construct_type[c_2nd_order])) {
    cat("\n\t", 
        col_align(i, max(l, nchar("Name")) + 2), 
        col_align(x22$Model$construct_type[i], 13 + 2), 
        col_align("Second order", 12 + 2), sep = "")
    if(x12$Arguments$.approach_weights == "PLS-PM") {
      cat(col_align(x22$Weight_info$Modes[i], 5), sep = "")
    }
  }
  
  ### Estimates ----------------------------------------------------------------
  cat("\n\n", rule(center = "Estimates"), "\n\n", sep = "")
  
  ## Path estimates
  cat("Estimated Path Coefficients:\n============================", sep = "")
  l <- max(nchar(x21$Path_estimates[, "Name"]))
  
  cat("\n\t", 
      col_align("Path", max(l, nchar("Path")) + 2), 
      col_align("Estimate", 10, align = "right"), 
      col_align("Std. Error", 12, align = "right"),
      col_align("t-value", 10, align = "right"), sep = "")
  
  for(i in 1:nrow(x21$Path_estimates)) {
    cat("\n\t", 
        col_align(x21$Path_estimates[i, "Name"], max(l, nchar("Path")) + 2), 
        col_align(sprintf("%.4f", x21$Path_estimates[i, "Estimate"]), 10, align = "right"),
        col_align(NA, 12, align = "right"),
        col_align(NA, 10, align = "right"),
        sep = "")
  }
  
  ## Loadings
  cat("\n\nEstimated Loadings:\n===================", sep = "")
  l <- max(nchar(c(x11$Loading_estimates[, "Name"], x21$Loading_estimates[, "Name"])))
  
  cat("\n\t", 
      col_align("Loading", max(l, nchar("Loading")) + 2), 
      col_align("Estimate", 10, align = "right"), 
      col_align("Std. Error", 12, align = "right"),
      col_align("t-value", 10, align = "right"), sep = "")
  # First Stage
  for(i in 1:nrow(x11$Loading_estimates)) {
    cat("\n\t", 
        col_align(x11$Loading_estimates[i, "Name"], max(l, nchar("Loading")) + 2), 
        col_align(sprintf("%.4f", x11$Loading_estimates[i, "Estimate"]), 10, align = "right"), 
        col_align(NA, 12, align = "right"),
        col_align(NA, 10, align = "right"),
        sep = "")
  }
  # Second stage
  for(i in 1:nrow(x21$Loading_estimates)) {
    cat("\n\t", 
        col_align(x21$Loading_estimates[i, "Name"], max(l, nchar("Loading")) + 2), 
        col_align(sprintf("%.4f", x21$Loading_estimates[i, "Estimate"]), 10, align = "right"), 
        col_align(NA, 12, align = "right"),
        col_align(NA, 10, align = "right"),
        sep = "")
  }
  
  ## Only print weights, for constructs modeled as common factors
  temp_w1 <- x11$Weight_estimates[x11$Weight_estimates$Construct_type == "Common factor", , drop = FALSE] 
  temp_w2 <- x21$Weight_estimates[x21$Weight_estimates$Construct_type == "Common factor", , drop = FALSE] 
  
  cat("\n\nEstimated Weights:\n==================\n", sep = "")
  l <- max(nchar(c(temp_w1[, "Name"], temp_w2[, "Name"])))
  
  cat("\n\t", 
      col_align("Weight", max(l, nchar("Weight")) + 2), 
      col_align("Estimate", 10, align = "right"), sep = "")
  # First stage
  for(i in 1:nrow(temp_w1)) {
    cat("\n\t", 
        col_align(temp_w1[i, "Name"], max(l, nchar("Weight")) + 2), 
        col_align(sprintf("%.4f", temp_w1[i, "Estimate"]), 10, align = "right"), 
        sep = "")
  }
  # Second stage
  for(i in 1:nrow(temp_w2)) {
    cat("\n\t", 
        col_align(temp_w2[i, "Name"], max(l, nchar("Weight")) + 2), 
        col_align(sprintf("%.4f", temp_w2[i, "Estimate"]), 10, align = "right"), 
        sep = "")
  }
  # if(x$Weight_estimator == "PLS-PM") {
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

#' `cSEMVerify_default` method for `print()`
#'
#' The `cSEMVerify_default` method for the generic function [print()]. 
#'
#' @usage print(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#'
#' @export
#'
print.cSEMVerify_default <- function(.object) {
  
  cat(rule(line = "bar2"), sep = "")
  
  cat("\n\nVerify admissibility:\n", sep = "")
  
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
  cat("\nDetails:\n\n", sep = "")
  
  cat("  ", col_align("Code", 7), col_align("Status", 10), "Description\n", sep = "")
  for(i in names(.object)) {
    cat("  ", col_align(i, 7), 
        col_align(ifelse(.object[i] == FALSE, green("ok"), red("not ok")), 10), 
        col_align(text[i], max(nchar(text)) + 2, align = "left"), "\n", sep = "")
  }
  cat(rule(line = "bar2"), sep = "")
}

#' `cSEMVerify_2ndorder` method for `print()`
#'
#' The `cSEMVerify_default` method for the generic function [print()]. 
#'
#' @usage print(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#'
#' @export
#'
print.cSEMVerify_2ndorder <- function(.object) {
  
  cat(rule(line = "bar2"), sep = "")
  cat("\n\nVerify admissibility:\n", sep = "")
  
  if(sum(.object$First_stage) == 0 & sum(.object$Second_stage) == 0) {
    cat(green("\n\t admissible"), sep = "")
  } else {
    cat(red("\n\t inadmissible"), sep = "")
  }
  
  text <- c("1" = "Convergence", 
            "2" = "At least one standardized loadings > 1", 
            "3" = "Construct VCV not positive semi-definite", 
            "4" = "Model-implied VCV not positive semi-definite",
            "5" = "At least one construct reliability > 1")
  
  cat("\n\nDetails:\n\n", sep = "")
  
  cat("  ", col_align("Code", 7), col_align("Stage 1", 10), 
      col_align("Stage 2", 10), "Description\n", sep = "")
  for(i in names(.object$First_stage)) {
    cat("  ", col_align(i, 7), 
        col_align(ifelse(.object$First_stage[i] == FALSE, green("ok"), red("not ok")), 10), 
        col_align(ifelse(.object$Second_stage[i] == FALSE, green("ok"), red("not ok")), 10), 
        col_align(text[i], max(nchar(text) + 2), align = "left"), "\n", sep = "")
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
      col_align("Critical value", 8*ncol(.object$Critical_value), align = "center"), "\n\t",
      col_align("Distance measure", width = 20, align = "left"),
      col_align("Test statistic", 18, align = "left"), "\t",
      sep = "")
  
  for(i in colnames(.object$Critical_value)) {
    cat(col_align(i, 6, align = "center"), "\t", sep = "")
  }
  cat("\n\t")
  for(j in seq_along(.object$Test_statistic)) {
    cat(col_align(names(.object$Test_statistic)[j], 20),
        col_align(sprintf("%.4f", .object$Test_statistic[j]), 18), "\t", sep = "")
    for(k in 1:ncol(.object$Critical_value)) {
      cat(sprintf("%.4f", .object$Critical_value[j, k]), "\t",
          sep = "")
    }
    cat("\n\t")
  }
  
  cat("\n\nDecision: \n\n\t", sep = "")
  
  cat(col_align("", width = 20, align = "left"), "\t",
      col_align("Significance level", 8*ncol(.object$Critical_value), align = "center"), "\n\t",
      col_align("Distance measure", width = 20, align = "left"), "\t",
      sep = "")
  
  for(i in colnames(.object$Critical_value)) {
    cat(col_align(i, 8, align = "center"), "\t", sep = "")
  }
  cat("\n\t")
  for(j in seq_along(.object$Test_statistic)) {
    cat(col_align(names(.object$Test_statistic)[j], 20), "\t", sep = "")
    for(k in 1:ncol(.object$Critical_value)) {
      cat(col_align(ifelse(.object$Decision[j, k], green("Do not reject"), red("reject")), 8), "\t", sep = "")
    }
    cat("\n\t")
  }

  cat("\nAdditonal information:")
  cat("\n\n\tOut of ", .object$Total_runs , " bootstrap replications ", .object$Number_admissibles, " are admissible.\n\t",
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
      col_align("Critical value", 8*ncol(.object$Critical_value), align = "center"), "\n\t",
      col_align("Distance measure", width = 20, align = "left"),
      col_align("Test statistic", 18, align = "left"), "\t",
      sep = "")
  
  for(i in colnames(.object$Critical_value)) {
    cat(col_align(i, 6, align = "center"), "\t", sep = "")
  }
  cat("\n\t")
  for(j in seq_along(.object$Test_statistic)) {
    cat(col_align(names(.object$Test_statistic)[j], 20),
        col_align(sprintf("%.4f", .object$Test_statistic[j]), 18), "\t", sep = "")
    for(k in 1:ncol(.object$Critical_value)) {
      cat(sprintf("%.4f", .object$Critical_value[j, k]), "\t",
          sep = "")
    }
    cat("\n\t")
  }
  
  cat("\n\nDecision: \n\n\t", sep = "")
  
  cat(col_align("", width = 20, align = "left"), "\t",
      col_align("Significance level", 8*ncol(.object$Critical_value), align = "center"), "\n\t",
      col_align("Distance measure", width = 20, align = "left"), "\t",
      sep = "")
  
  for(i in colnames(.object$Critical_value)) {
    cat(col_align(i, 8, align = "center"), "\t", sep = "")
  }
  cat("\n\t")
  for(j in seq_along(.object$Test_statistic)) {
    cat(col_align(names(.object$Test_statistic)[j], 20), "\t", sep = "")
    for(k in 1:ncol(.object$Critical_value)) {
      cat(col_align(ifelse(.object$Decision[k, j], green("Do not reject"), red("reject")), 8), "\t", sep = "")
    }
    cat("\n\t")
  }
  
  cat("\nAdditonal information:")
  cat("\n\n\tOut of ", .object$Total_runs , " permutation runs, ", 
      .object$Number_admissibles, " where admissible.\n\t",
      "See ", yellow("?"), magenta("verify"), "()",
      " for what constitutes an inadmissible result.\n", sep = "")
  cat(rule(line = "bar2"), sep = "")
}

#' `cSEMTestMICOM` method for `print()`
#'
#' The `cSEMTestMICOM` method for the generic function [print()]. 
#'
#' @usage print(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cca()], [foreman()], [cSEMResults]
#'
#' @export
#'
# print.testMICOM <- function(x, ...) {
# 
#   cat(cli::rule(), "\n")
#   cat(cli::rule(center = "Overview", line = "bar3"), "\n\n",
#       crayon::col_align("\tNumber of Observations", 25), "= ", x$Meta_information$Number_of_observations[1, 2], "\n", sep = "")
#   for(i in 2:nrow(x$Meta_information$Number_of_observations)) {
#     cat("\t\t", x$Meta_information$Number_of_observations[i, "x"], " : ",
#     x$Meta_information$Number_of_observations[i, "n"], "\n")
#   }
#   cat(
#       crayon::col_align("\tNumber of Groups", 25), "= ", x$Meta_information$Number_of_Groups, "\n",
#       crayon::col_align("\tGrouping Variable", 25), "= ", x$Meta_information$Grouping_variable, "\n\n",
#       sep = "")
#   cat(cli::rule(center = "Details", line = "bar3"), "\n")
#   cat(cli::rule(center = "Step 1 - Configural invariance", line = 2), "\n\n",
#       "\tConfigural invariance is a precondition for step 2 and 3.\n",
#       "\tDo not proceed to interpret results unless\n",
#       "\tconfigural invariance has been established.\n\n",
#       sep = "")
#   cat(cli::rule(center = "Step 2 - Compositional invariance", line = 2), "\n\n",
#       cli::boxx("H0: Compositional measurement invariance holds", float = "center"), "\n\n",
#       sep = "")
# 
#   l <- max(nchar(c("Construct", rownames(x$Step2[[1]]))))
# 
#     for(i in seq_along(x$Step2)) {
# 
#     cat("Groups: ", names(x$Step2)[i], "\n\t",
#         crayon::col_align("", width = l + 2, align = "center"),
#         crayon::col_align("", 10, align = "center"), "\t",
#         crayon::col_align("Critical Value(s)", 8*(ncol(x$Step2[[1]]) - 1), align = "center"), "\n\t",
#         crayon::col_align("Construct", width = l + 2, align = "center"),
#         crayon::col_align("c", 10, align = "center"), "\t",
#         sep = "")
#     for(j in colnames(x$Step2[[1]])[-1]) {
#       cat(crayon::col_align(j, 6, align = "center"), "\t", sep = "")
#     }
#     cat("\n\t")
# 
#     for(j in 1:nrow(x$Step2[[i]])) {
#       cat(crayon::col_align(row.names(x$Step2[[i]])[j], l + 2), ": ",
#           sprintf("%7.4f", x$Step2[[i]][j, "c"]) , "\t", sep = "")
#       for(k in 2:ncol(x$Step2[[i]])) {
#         cat(sprintf("%7.4f", x$Step2[[i]][j, k]), "\t",
#             sep = "")
#       }
#       cat("\n\t")
#     }
#     cat("\n")
#   }
# 
#   cat(cli::rule(center = "Step 3 - Equality of the mean values and variances", line = 2), "\n\n",
#       cli::boxx(c("1. H0: Difference between group means is zero",
#                   "2. H0: Log of the ratio of the group variances is zero"),
#                 float = "center"),
#       sep = "")
# 
#   cat("\n\nEquality of the means:\n", "______________________", sep = "")
#   for(i in seq_along(x$Step3$Mean_diff)) {
# 
#     cat("\n\nGroups: ", names(x$Step3$Mean_diff)[i], "\n\t",
#         crayon::col_align("", width = l + 2, align = "center"),
#         crayon::col_align("", 10, align = "center"), "\t",
#         crayon::col_align("Critical Value(s)", 8*(ncol(x$Step3$Mean_diff[[1]]) - 1), align = "center"), "\n\t",
#         crayon::col_align("Construct", width = l + 2, align = "center"),
#         crayon::col_align("Mean diff.", 11, align = "center"), "\t",
#         sep = "")
# 
#     for(j in colnames(x$Step3$Mean_diff[[1]][-1])) {
#       cat(crayon::col_align(j, 6, align = "center"), "\t", sep = "")
#     }
#     cat("\n\t")
# 
#     for(j in 1:nrow(x$Step3$Mean_diff[[i]])) {
#       cat(crayon::col_align(row.names(x$Step3$Mean_diff[[i]])[j], l + 2), ": ",
#           crayon::col_align(sprintf("%7.4f", x$Step3$Mean_diff[[i]][j, "Diff_mean"]), 11), sep = "")
#       for(k in 2:ncol(x$Step3$Mean_diff[[i]])) {
#         cat(sprintf("%7.4f", x$Step3$Mean_diff[[i]][j, k]), "\t",
#             sep = "")
#       }
#       cat("\n\t")
#     }
#     cat("\n")
#   }
# 
#   cat("\n\nEquality of the variances:\n", "__________________________", sep = "")
#   for(i in seq_along(x$Step3$Var_diff)) {
# 
#     cat("\n\nGroups: ", names(x$Step3$Var_diff)[i], "\n\t",
#         crayon::col_align("", width = l + 2, align = "center"),
#         crayon::col_align("", 10, align = "center"), "\t",
#         crayon::col_align("Critical Value(s)", 8*(ncol(x$Step3$Var_diff[[1]]) - 1), align = "center"), "\n\t",
#         crayon::col_align("Construct", width = l + 2, align = "center"),
#         crayon::col_align("Var diff.", 11, align = "center"), "\t",
#         sep = "")
# 
#     for(j in colnames(x$Step3$Var_diff[[1]][-1])) {
#       cat(crayon::col_align(j, 6, align = "center"), "\t", sep = "")
#     }
#     cat("\n\t")
# 
#     for(j in 1:nrow(x$Step3$Var_diff[[i]])) {
#       cat(crayon::col_align(row.names(x$Step3$Var_diff[[i]])[j], l + 2), ": ",
#           crayon::col_align(sprintf("%7.4f", x$Step3$Var_diff[[i]][j, "Diff_log_var"]), 11), sep = "")
#       for(k in 2:ncol(x$Step3$Var_diff[[i]])) {
#         cat(sprintf("%7.4f", x$Step3$Var_diff[[i]][j, k]), "\t",
#             sep = "")
#       }
#       cat("\n\t")
#     }
#     cat("\n")
#   }
# }

