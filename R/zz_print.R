#' `cSEMResults` method for `print()`
#'
#' The [cSEMResults] method for the generic function [print()].
#'
#' @usage print(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
#' @keywords internal
print.cSEMResults <- function(.object, ...) {
    
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
          ifelse(sum(x) == 0, green("successful"), red("not successful")), sep = "")
    }
    if(sum(x) != 0) {
      cat("\n\nSee ", magenta("verify"), "(", cyan("<object-name>"), ")", 
          " for details.", sep = "")
    }
    cat("\n\nThe result for each stage is a list of class " %+% bold("cSEMResults") %+%"",
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
      "- ", magenta("summarize"), "(", cyan("<object-name>"), ")\n\t",
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
#' @seealso [csem()], [foreman()], [cSEMResults], [summarize()]
#'
#' @export
#' @keywords internal
print.cSEMSummarize_default <- function(.object, .full_output = FALSE, ...) {
  
  x1 <- .object$Estimates
  x2 <- .object$Information
  
  cat2(
    rule(line = "bar2"), "\n",
    rule(center = "Overview"), 
    "\n"
  )
  
  ### Overview -----------------------------------------------------------------
  cat2(
    col_align("\n\tNumber of observations", 35), "= ", nrow(x2$Arguments$.data),
    col_align("\n\tWeight estimator", 35), "= ", 
    ifelse(x2$Arguments$.approach_weights == "PLS-PM" && 
             x2$Type_of_indicator_correlation %in% c("Polychoric", "Polyserial"), 
           "PLS-PM (OrdPLS)", x2$Arguments$.approach_weights)
  )
  
  if(x2$Arguments$.approach_weights == "PLS-PM") {
    cat2(
      col_align("\n\tInner weighting scheme", 35), "= ", 
       x2$Arguments$.PLS_weight_scheme_inner
    )
  }
  cat2(
    col_align("\n\tType of indicator correlation", 35), "= ", 
    paste0(x2$Type_of_indicator_correlation, collapse = ", ")
    )
  cat(
    col_align("\n\tPath model estimator", 35), "= ", x2$Arguments$.approach_paths,
    col_align("\n\tType of path model", 35), "= ", x2$Model$model_type,
    col_align("\n\tDisattenuated", 35), "= ", 
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
  
  # Get the column names of the columns containing confidence intervals
  ci_colnames <- colnames(x1$Path_estimates)[-c(1:6)]
  
  # Are there confidence intervals
  if(length(ci_colnames) != 0) {
    
    cat(
      paste0("Inference based on ", x2$Information_resample$Method, " resampling."),
      sep = ""
    )
    
    if(length(ci_colnames) > 2 && !.full_output) {
      cat("\n",
        "Only the first confidence interval shown by default.\n",
        "Use `print(.object, .full_output = TRUE)` to print all confidence intervals.",
        sep = "")
      ci_colnames <- ci_colnames[1:2]
    }
    cat("\n\n")
  }

  ## Path estimates
  cat("Estimated path coefficients:\n============================", sep = "")
  l <- max(nchar(x1$Path_estimates[, "Name"]))
 
  if(length(ci_colnames) != 0) {
    xx <- regmatches(ci_colnames, regexpr("\\.", ci_colnames), invert = TRUE)
    interval_names    <- unique(sapply(xx, `[`, 1))
    sig_level_names   <- unique(gsub("[LU]", "", sapply(xx, `[`, 2)))
    
    cat("\n  ", 
        col_align("", width = max(l, nchar("Path")) + 44), 
        sep = "")
    for(i in interval_names) {
      cat(col_align(i, width = 20*length(sig_level_names), align = "center"),
          sep = "")
    }
  }
  cat("\n  ", 
      col_align("Path", max(l, nchar("Path")) + 2), 
      col_align("Estimate", 10, align = "right"), 
      col_align("Std. error", 12, align = "right"),
      col_align("t-stat.", 10, align = "right"), 
      col_align("p-value", 10, align = "right"),
      sep = "")
  if(length(ci_colnames) != 0) {
    for(i in rep(sig_level_names, length(interval_names))) {
      cat(
        col_align(i, 20, align = "center"),
        sep = "" 
      )
    } 
  }
  
  for(i in 1:nrow(x1$Path_estimates)) {
    cat("\n  ", 
        col_align(x1$Path_estimates[i, "Name"], max(l, nchar("Path")) + 2), 
        col_align(sprintf("%.4f", x1$Path_estimates[i, "Estimate"]), 10, align = "right"),
        col_align(sprintf("%.4f", x1$Path_estimates[i, "Std_err"]), 12, align = "right"),
        col_align(sprintf("%.4f", x1$Path_estimates[i, "t_stat"]), 10, align = "right"),
        col_align(ifelse(x1$Path_estimates[i, "p_value"] < 0.05, 
                     green(sprintf("%.4f", x1$Path_estimates[i, "p_value"])),
                     sprintf("%.4f", x1$Path_estimates[i, "p_value"])), 10, align = "right"),
        sep = "")
    if(length(ci_colnames) != 0) {
      for(j in seq(1, length(ci_colnames), by = 2) + 6) {
        cat(
          col_align(
            paste0("[", sprintf("%7.4f", x1$Path_estimates[i, j]), ";", 
                   sprintf("%7.4f", x1$Path_estimates[i, j+1]), "]"), 20, align = "center"),
          sep = "" 
        )
      } 
    }
  }
  
  ## Loadings
  cat("\n\nEstimated Loadings:\n===================", sep = "")
  l <- max(nchar(x1$Loading_estimates[, "Name"]))
  
  if(length(ci_colnames) != 0) {
    cat("\n  ", 
        col_align("", width = max(l, nchar("Loading")) + 44), 
        sep = "")
    for(i in interval_names) {
      cat(col_align(i, width = 20*length(sig_level_names), align = "center"),
          sep = "")
    }
  }
  
  cat("\n  ", 
      col_align("Loading", max(l, nchar("Loading")) + 2), 
      col_align("Estimate", 10, align = "right"), 
      col_align("Std. error", 12, align = "right"),
      col_align("t-stat.", 10, align = "right"), 
      col_align("p-value", 10, align = "right"),
      sep = "")
  if(length(ci_colnames) != 0) {
    for(i in rep(sig_level_names, length(interval_names))) {
      cat(
        col_align(i, 20, align = "center"),
        sep = "" 
      )
    } 
  }
  
  for(i in 1:nrow(x1$Loading_estimates)) {
    cat("\n  ", 
        col_align(x1$Loading_estimates[i, "Name"], max(l, nchar("Loading")) + 2), 
        col_align(sprintf("%.4f", x1$Loading_estimates[i, "Estimate"]), 10, align = "right"), 
        col_align(sprintf("%.4f", x1$Loading_estimates[i, "Std_err"]), 12, align = "right"),
        col_align(sprintf("%.4f", x1$Loading_estimates[i, "t_stat"]), 10, align = "right"),
        col_align(ifelse(x1$Loading_estimates[i, "p_value"] < 0.05, 
                         green(sprintf("%.4f", x1$Loading_estimates[i, "p_value"])),
                         sprintf("%.4f", x1$Loading_estimates[i, "p_value"])), 10, align = "right"),
        sep = "")
    if(length(ci_colnames) != 0) {
      for(j in seq(1, length(ci_colnames), by = 2) + 6) {
        cat(
          col_align(
            paste0("[", sprintf("%7.4f", x1$Loading_estimates[i, j]), ";", 
                   sprintf("%7.4f", x1$Loading_estimates[i, j+1]), "]"), 20, align = "center"),
          sep = "" 
        )
      } 
    }
  }

  ## Only print weights, for constructs modeled as composites
  temp_w <- x1$Weight_estimates[x1$Weight_estimates$Construct_type == "Composite", , drop = FALSE] 
  
  if(nrow(temp_w) != 0) {
    cat("\n\nEstimated Weights:\n==================\n", sep = "")
    l <- max(nchar(temp_w[, "Name"]))
    
    if(length(ci_colnames) != 0) {
      cat("\n  ", 
          col_align("", width = max(l, nchar("Weights")) + 44), 
          sep = "")
      for(i in interval_names) {
        cat(col_align(i, width = 20*length(sig_level_names), align = "center"),
            sep = "")
      }
    }
    
    cat("\n  ", 
        col_align("Weights", max(l, nchar("Loading")) + 2), 
        col_align("Estimate", 10, align = "right"), 
        col_align("Std. error", 12, align = "right"),
        col_align("t-stat.", 10, align = "right"), 
        col_align("p-value", 10, align = "right"),
        sep = "")
    if(length(ci_colnames) != 0) {
      for(i in rep(sig_level_names, length(interval_names))) {
        cat(
          col_align(i, 20, align = "center"),
          sep = "" 
        )
      } 
    }
    
    for(i in 1:nrow(temp_w)) {
      cat("\n  ", 
          col_align(temp_w[i, "Name"], max(l, nchar("Loading")) + 2), 
          col_align(sprintf("%.4f", temp_w[i, "Estimate"]), 10, align = "right"), 
          col_align(sprintf("%.4f", temp_w[i, "Std_err"]), 12, align = "right"),
          col_align(sprintf("%.4f", temp_w[i, "t_stat"]), 10, align = "right"),
          col_align(ifelse(temp_w[i, "p_value"] < 0.05, 
                           green(sprintf("%.4f", temp_w[i, "p_value"])),
                           sprintf("%.4f", temp_w[i, "p_value"])), 10, align = "right"),
          sep = "")
      if(length(ci_colnames) != 0) {
        for(j in seq(1, length(ci_colnames), by = 2) + 6) {
          cat(
            col_align(
              paste0("[", sprintf("%7.4f", temp_w[i, j]), ";", 
                     sprintf("%7.4f", temp_w[i, j+1]), "]"), 20, align = "center"),
            sep = "" 
          )
        } 
      }
    }
  }

  
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
#' @seealso [csem()], [foreman()], [cSEMResults], [summarize()]
#'
#' @export
#' @keywords internal
print.cSEMSummarize_2ndorder <- function(.object, ...) {
  
  ### Extract name and objects 
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
  cat2(
    col_align("\n\tNumber of observations", 35), "= ", nrow(x12$Arguments$.data),
    col_align("\n\tWeight estimator", 35), "= ", 
    ifelse(x12$Arguments$.approach_weights == "PLS-PM" && 
             x12$Type_of_indicator_correlation %in% c("Polychoric", "Polyserial"), 
           "PLS-PM (OrdPLS)", x12$Arguments$.approach_weights)
  )
  
  if(x12$Arguments$.approach_weights == "PLS-PM") {
    cat(
      col_align("\n\tInner weighting scheme", 35), "= ", 
      x12$Arguments$.PLS_weight_scheme_inner, 
      sep = "")
  }
  cat2(
    col_align("\n\tType of indicator correlation", 35), "= ", 
    paste0(x12$Type_of_indicator_correlation, collapse = ", ")
  )
  cat(
    col_align("\n\tPath model estimator", 35), "= ", x12$Arguments$.approach_paths,
    # Its important to read the model type of the second stage (since the first
    # is always linear)
    col_align("\n\tType of path model", 35), "= ", x22$Model$model_type,
    col_align("\n\tDisattenuated", 35), "= ", 
    ifelse(x12$Arguments$.disattenuate, 
           ifelse(x12$Arguments$.approach_weights == "PLS-PM", "Yes (PLSc)", "Yes"), "No"),
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
  
  # cat("\n\n", rule(center = "Other output"), "\n\n\t", sep = "")
  # 
  # cat("<not yet implemented>")
  # 
  # cat("\n\n", rule(center = "Fit Indices"), "\n\n\t", sep = "")
  # 
  # cat(col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\t",
  #     col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\t",
  #     col_align("<some_index>", 30), "= ", c("not yet implemented"), "\n\n",
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
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
#' @keywords internal
print.cSEMVerify_default <- function(.object, ...) {
  
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
            "5" = "At least one proxy reliability > 1")
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
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
#' @keywords internal
print.cSEMVerify_2ndorder <- function(.object, ...) {
  
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
            "5" = "At least one proxy reliability > 1")
  
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
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
#' @keywords internal
print.cSEMTestOMF <- function(.object, ...) {
  
  cat(
    rule(line = "bar2"), "\n",
    rule(center = "Test for overall model fit based on Dijkstra & Henseler (2015)"), 
    sep = ""
  )
  
  ## Null hypothesis -----------------------------------------------------------
  cat(
    "\n\nNull hypothesis:\n\n", 
    boxx(c("H0: No significant difference between empirical and", 
           "model-implied indicator covariance matrix."), float = "center"), 
    sep = ""
  )
  
  ## Test statistic and critical value -----------------------------------------
  cat("\n\nTest statistic and critical value: \n\n\t", sep = "")
  
  cat(
    col_align("", width = 20),
    col_align("", width = 14), 
    "\t",
    col_align("Critical value", width = 8*ncol(.object$Critical_value), 
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
  
  for(i in colnames(.object$Critical_value)) {
    cat(col_align(i, width = 6, align = "center"), "\t", sep = "")
  }
  
  cat("\n\t")
  
  for(j in seq_along(.object$Test_statistic)) {
    cat(
      col_align(names(.object$Test_statistic)[j], width = 20),
      col_align(sprintf("%.4f", .object$Test_statistic[j]), width = 14, 
                align = "center"), 
      "\t", 
      sep = ""
    )
    for(k in 1:ncol(.object$Critical_value)) {
      cat(sprintf("%.4f", .object$Critical_value[j, k]), "\t", sep = "")
    }
    cat("\n\t")
  }
  
  ## Decision ------------------------------------------------------------------
  cat("\n\nDecision: \n\n\t", sep = "")
  
  # Width of columns
  l <- apply(.object$Decision, 2, function(x) {
    ifelse(any(x == TRUE), nchar("Do not reject"), nchar("reject"))
  })
  
  l1 <- max(c(sum(l) + 3*(ncol(.object$Decision) - 1)), nchar("Significance level"))
  
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
  
  for(i in colnames(.object$Critical_value)) {
    cat(col_align(i, width = l[i], align = "center"), "\t", sep = "")
  }
  
  cat("\n\t")
  
  for(j in seq_along(.object$Test_statistic)) {
    
    cat(col_align(names(.object$Test_statistic)[j], width = 20), "\t", sep = "")
    
    for(k in 1:ncol(.object$Critical_value)) {
      cat(
        col_align(ifelse(.object$Decision[j, k], 
                         green("Do not reject"), red("reject")), 
                  width = l[k], align = "center"), 
        "\t", 
        sep = ""
      )
    }
    cat("\n\t")
  }

  ## Additional information ----------------------------------------------------
  cat("\nAdditonal information:")
  cat(
    "\n\n\tOut of ", .object$Information$Total_runs , " bootstrap replications ", 
    .object$Information$Number_admissibles, " are admissible.\n\t",
    "See ", yellow("?"), magenta("verify"), "()",
    " for what constitutes an inadmissible result.\n", 
    sep = ""
  )
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
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
#' @keywords internal
print.cSEMTestMGD <- function(.object, ...) {

  cat(
    rule(line = "bar2"), "\n",
    rule(center = "Test for multigroup differences based on Klesel (forthcoming)"), 
    sep = ""
  )
  
  ## Null hypothesis -----------------------------------------------------------
  cat(
    "\n\nNull hypothesis:\n\n", 
    boxx("H0: No significant difference between groups.", float = "center"), 
    sep = ""
  )
  
  ## Test statistic and critical value -----------------------------------------
  cat("\n\nTest statistic and critical value: \n\n\t", sep = "")
  
  cat(
    col_align("", width = 20),
    col_align("", width = 14), 
    "\t",
    col_align("Critical value", width = 8*ncol(.object$Critical_value), 
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
  
  for(i in colnames(.object$Critical_value)) {
    cat(col_align(i, width = 6, align = "center"), "\t", sep = "")
  }
  
  cat("\n\t")
  
  for(j in seq_along(.object$Test_statistic)) {
    cat(
      col_align(names(.object$Test_statistic)[j], width = 20),
      col_align(sprintf("%.4f", .object$Test_statistic[j]), width = 14, 
                align = "center"), 
      "\t", 
      sep = ""
    )
    for(k in 1:ncol(.object$Critical_value)) {
      cat(sprintf("%.4f", .object$Critical_value[j, k]), "\t", sep = "")
    }
    cat("\n\t")
  }
  
  ## Decision ------------------------------------------------------------------
  cat("\n\nDecision: \n\n\t", sep = "")
  
  # Width of columns
  l <- apply(.object$Decision, 2, function(x) {
    ifelse(any(x == TRUE), nchar("Do not reject"), nchar("reject"))
  })
  
  l1 <- max(c(sum(l) + 3*(ncol(.object$Decision) - 1)), nchar("Significance level"))
  
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
  
  for(i in colnames(.object$Critical_value)) {
    cat(col_align(i, width = l[i], align = "center"), "\t", sep = "")
  }
  
  cat("\n\t")
  
  for(j in seq_along(.object$Test_statistic)) {
    
    cat(col_align(names(.object$Test_statistic)[j], width = 20), "\t", sep = "")
    
    for(k in 1:ncol(.object$Critical_value)) {
      cat(
        col_align(ifelse(.object$Decision[j, k], 
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
  cat(
    "\n\n\tOut of ", .object$Information$Total_runs , " permutation runs, ", 
    .object$Information$Number_admissibles, " where admissible.\n\t",
    "See ", yellow("?"), magenta("verify"), "()",
    " for what constitutes an inadmissible result.", 
    sep = ""
  )
  
  cat("\n\n\tNumber of observations per group:")
  
  l <- max(nchar(c(.object$Information$Group_names, "Group")))
  
  cat("\n\n\t",
      col_align("Group", width = l + 6),
      col_align("No. observations", width = 15),
      sep = ""
  )
  for(i in seq_along(.object$Information$Group_names)) {
    cat(
      "\n\t",
      col_align(.object$Information$Group_names[i], width = l + 6), 
      col_align(.object$Information$Number_of_observations[i], width = 15),
      sep = ""
    )
  }
  
  cat("\n", rule(line = "bar2"), sep = "")
}

#' `cSEMTestMICOM` method for `print()`
#'
#' The `cSEMTestMICOM` method for the generic function [print()]. 
#'
#' @usage print(.object)
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
#' @keywords internal
print.cSEMTestMICOM <- function(.object, ...) {

  cat(
    rule(line = "bar2"), "\n",
    rule(center = "Test for measurement invariance based on Henseler et al (2016)"),
    "\n",
    sep = ""
  )
  
  x2 <- .object$Step2
  x3 <- .object$Step3
  
  construct_names <- names(x2$Test_statistic[[1]])
  no_groups       <- length(x2$Test_statistic)
  l_names <- max(nchar(c("Construct", construct_names)))
  
  ### Step 1 ===================================================================
  cat(
    rule(center = "Step 1 - Configural invariance", line = 2), "\n\n",
    "\tConfigural invariance is a precondition for step 2 and 3.\n",
    "\tDo not proceed to interpret results unless\n",
    "\tconfigural invariance has been established.\n\n",
    sep = ""
  )
  ### Step 2 ===================================================================
  cat(rule(center = "Step 2 - Compositional invariance", line = 2))
  
  ## Null hypothesis -----------------------------------------------------------
  cat(
    "\n\nNull hypothesis:\n\n",
    boxx("H0: Compositional measurement invariance of the constructs.", 
         float = "center"), 
    sep = ""
  )
  
  ## Test statistic and critical value -----------------------------------------
  cat("\n\nTest statistic and critical value: \n\n", sep = "")

  for(i in 1:no_groups) {
    
    cat("  Compared groups: ", names(x2$Test_statistic)[i], "\n\t", sep = "")
    cat(
      col_align("", width = l_names),
      "\t",
      col_align("", width = 14), 
      "\t",
      col_align("Critical value", width = 8*ncol(x2$Critical_value[[i]]), 
                align = "center"),
      sep = ""
    )
    cat(
      "\n\t",
      col_align("Construct", width = l_names),
      "\t",
      col_align("Test statistic", width = 14), 
      "\t",
      sep = ""
    )

    for(j in colnames(x2$Critical_value[[i]])) {
      cat(col_align(j, width = 6, align = "center"), "\t", sep = "")
    }
    
    cat("\n\t")
    
    for(j in construct_names) {
      cat(
        col_align(j, width = l_names), 
        "\t",
        col_align(sprintf("%.4f", x2$Test_statistic[[i]][j]), width = 14,
                  align = "center"),
        "\t",
        sep = ""
      )
      
      for(k in colnames(x2$Critical_value[[i]])) {
        cat(sprintf("%.4f", x2$Critical_value[[i]][j, k]), "\t", sep = "")
      } # END for j (each construct)
      cat("\n\t")
    } # END for k (each significance level)
    cat("\n\n")
  } # END for i (each group)
  
  ## Decision ------------------------------------------------------------------
  cat("Decision: \n\n", sep = "")
  
  for(i in 1:no_groups) {
    # Width of columns
    l <- apply(x2$Decision[[i]], 2, function(x) {
      ifelse(any(x == TRUE), nchar("Do not reject"), nchar("reject"))
    })
    
    l1 <- max(c(sum(l) + 3*(ncol(x2$Decision[[i]]) - 1)), nchar("Significance level"))
    cat("  Compared groups: ", names(x2$Test_statistic)[i], "\n\t")
    cat(
      col_align("", width = l_names), 
      "\t",
      col_align("Significance level", 
                width = l1, 
                align = "center"),
      sep = ""
    )
    cat(
      "\n\t",
      col_align("Construct", width = l_names), 
      "\t",
      sep = ""
    )
    
    for(j in colnames(x2$Critical_value[[i]])) {
      cat(col_align(j, width = l[j], align = "center"), "\t", sep = "")
    }
    
    cat("\n\t")
    
    for(j in seq_along(x2$Test_statistic[[i]])) {
      
      cat(col_align(names(x2$Test_statistic[[i]])[j], width = l_names), "\t", sep = "")
      
      for(k in 1:ncol(x2$Critical_value[[i]])) {
        cat(
          col_align(ifelse(x2$Decision[[i]][j, k], 
                           green("Do not reject"), red("reject")), 
                    width = l[k], align = "center"), 
          "\t", 
          sep = ""
        )
      }
      cat("\n\t")
    }  
    cat("\n\n")
  } # END for i (each group)
  
  ### Step 3 ===================================================================
  cat(rule(center = "Step 3 - Equality of the means and variances", line = 2))
  
  ## Null hypothesis -----------------------------------------------------------
  cat(
    "\n\nNull hypothesis:\n\n",
    boxx(c("1. H0: Difference between group means is zero",
           "2. H0: Log of the ratio of the group variances is zero"),
         float = "center"),
    sep = ""
  )

  ## Test statistic and critical value -----------------------------------------
  cat("\n\nTest statistic and critical values: \n\n", sep = "")
  
  for(i in 1:no_groups) {
    
    cat("  Compared groups: ", names(x2$Test_statistic)[i], "\n\n\t", sep = "")
    for(type in names(x3)) {
      x <- x3[[type]]
      cat(type, "\n\t", sep = "")
      cat(
        col_align("", width = l_names),
        "\t",
        col_align("", width = 14), 
        "\t",
        col_align("Critical values", width = 8*ncol(x$Critical_value[[i]]), 
                  align = "center"),
        sep = ""
      )
      cat(
        "\n\t",
        col_align("Construct", width = l_names),
        "\t",
        col_align("Test statistic", width = 14), 
        "\t",
        sep = ""
      )
      
      for(j in colnames(x$Critical_value[[i]])) {
        cat(col_align(j, width = 6, align = "center"), "\t", sep = "")
      }
      
      cat("\n\t")
      
      for(j in construct_names) {
        cat(
          col_align(j, width = l_names), 
          "\t",
          col_align(sprintf("%.4f", x$Test_statistic[[i]][j]), width = 14,
                    align = "center"),
          "\t",
          sep = ""
        )
        
        for(k in colnames(x$Critical_value[[i]])) {
          cat(sprintf("%.4f", x$Critical_value[[i]][j, k]), "\t", sep = "")
        } # END for j (each construct)
        cat("\n\t")
      } # END for k (each significance level)
      cat("\n\t")
    } # END for type (one of "Mean" and "Var")
    cat("\n")
  } # END for i (each group)
  
  ## Decision ------------------------------------------------------------------
  cat("Decision: \n\n", sep = "")
  
  for(i in 1:no_groups) {

    cat("  Compared groups: ", names(x2$Test_statistic)[i], "\n\n\t")
    
    for(type in names(x3)) {
      x <- x3[[type]]
      
      # Width of columns
      l <- apply(x$Decision[[i]], 2, function(x) {
        ifelse(any(x == TRUE), nchar("Do not reject"), nchar("reject"))
      })
      
      n1 <- colnames(x$Critical_value[[i]])
      n2 <- paste0("[", n1[seq(1, (length(n1) - 1), by = 2)], ";", 
                   n1[seq(2, length(n1), by = 2)], "]")
      
      l1 <- max(c(sum(l) + 3*(ncol(x$Decision[[i]]) - 1)), 
                nchar("Significance levels"),
                nchar(n2))
      
      cat(type, "\n\t", sep = "")
      cat(
        col_align("", width = l_names), 
        "\t",
        col_align("Significance level", 
                  width = l1, 
                  align = "center"),
        sep = ""
      )
      cat(
        "\n\t",
        col_align("Construct", width = l_names), 
        "\t",
        sep = ""
      )
      
      for(j in seq_along(n2)) {
        cat(col_align(n2[j], width = max(l[j], nchar(n2[j])), align = "center"),
            "\t", sep = "")
      }
      
      cat("\n\t")
      
      for(j in seq_along(x$Test_statistic[[i]])) {
        
        cat(col_align(names(x$Test_statistic[[i]])[j], width = l_names), "\t", sep = "")
        
        for(k in 1:ncol(x$Decision[[i]])) {
          cat(
            col_align(ifelse(x$Decision[[i]][j, k], 
                             green("Do not reject"), red("reject")),
                      width = max(l[k], nchar(n2[k])), align = "center"), 
            "\t", 
            sep = ""
          )
        } # END for k (each construct)
        cat("\n\t")
      } # END for j (each significance level)
      cat("\n\t")
    } # END for type (one of "Mean" and "Var")
    cat("\n")
  } # END for i (each group)
  
  ## Additional information ----------------------------------------------------
  cat("Additional information:")
  cat(
    "\n\n\tOut of ", .object$Information$Total_runs , " permutation runs, ", 
    .object$Information$Number_admissibles, " where admissible.\n\t",
    "See ", yellow("?"), magenta("verify"), "()",
    " for what constitutes an inadmissible result.", 
    sep = ""
  )
  
  cat("\n\n\tNumber of observations per group:")
  
  l <- max(nchar(c(.object$Information$Group_names, "Group")))
  
  cat("\n\n\t",
      col_align("Group", width = l + 6),
      col_align("No. observations", width = 15),
      sep = ""
  )
  for(i in seq_along(.object$Information$Group_names)) {
    cat(
      "\n\t",
      col_align(.object$Information$Group_names[i], width = l + 6), 
      col_align(.object$Information$Number_of_observations[i], width = 15),
      sep = ""
    )
  }
  
  cat("\n", rule(line = "bar2"), sep = "")
}

