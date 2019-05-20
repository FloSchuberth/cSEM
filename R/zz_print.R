#' `cSEMResults` method for `print()`
#'
#' The [cSEMResults] method for the generic function [print()].
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
#' @keywords internal
print.cSEMResults <- function(x, ...) {
  
  cat(
    rule(line = "bar2", width = 80), "\n",
    rule(center = "Overview", width = 80), 
    sep = "")
  
  if(class(x)[2] == "cSEMResults_multi") {
    cat("\n\nEstimation status by group/data set:\n", sep = "")
    for(i in names(x)) {
      cat("\n\t", col_align(cyan(i), 15), ": ",
          ifelse(sum(verify.cSEMResults_default(x[[i]])) == 0, 
                 green("successful"), red("not successful")), ".", sep = "")
    }
    if(sum(verify.cSEMResults_default(x[[i]])) != 0) {
      cat("\n\nSee ", magenta("verify"), "(", cyan("<object-name>"), ")", 
          " for details.", sep = "")
    }
    cat("\n\nThe result for each group/data set is a list of class " %+% bold("cSEMResults") %+%"",
        "\nwith list elements:\n\n\t", sep = "")
    
  } else if(class(x)[2] == "cSEMResults_2ndorder") {
    x <- sapply(verify.cSEMResults_2ndorder(x), sum)
    cat("\n\nEstimation status by stage:\n", sep = "")
    for(i in names(x)) {
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
      "\n\nEstimation was ", ifelse(sum(verify.cSEMResults_default(x)) == 0, 
                                    green("successful"), red("not successful")), ".", sep = "")
    if(sum(verify.cSEMResults_default(x)) != 0) {
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
      "- ", magenta("predict"), "(", cyan("<object-name>"), ")\n\t",
      "- ", magenta("summarize"), "(", cyan("<object-name>"), ")\n\t",
      "- ", magenta("verify"), "(", cyan("<object-name>"), ")\n", sep = "")
  cat(rule(line = "bar2", width = 80), "\n")
}

#' `cSEMSummarize` method for `print()`
#'
#' The [cSEMSummary] method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [foreman()], [cSEMResults], [summarize()]
#'
#' @export
#' @keywords internal
print.cSEMSummarize <- function(x, .full_output = TRUE, ...) {
  
  ## Check the class
  if(inherits(x, "cSEMSummarize_2ndorder")) {
    x11 <- x$First_stage$Estimates
    x12 <- x$First_stage$Information
    
    x21 <- x$Second_stage$Estimates
    x22 <- x$Second_stage$Information
  } else {
    
    x21 <- x$Estimates
    x22 <- x$Information
  }
  
  cat2(
    rule(line = "bar2", width = 80), "\n",
    rule(center = "Overview", width = 80), 
    "\n"
  )
  
  ### Overview -----------------------------------------------------------------
  ## General information + resample information
  printSummarizeOverview(x)
  
  ## Construct details
  cat2("\n\n\tConstruct details:\n\t","------------------")
  
  printSummarizeConstructDetails(x)
  
  ### Estimates ----------------------------------------------------------------
  cat2("\n\n", rule(center = "Estimates", width = 80), "\n\n")
  
  ## Confidence intervals
  # Get the column names of the columns containing confidence intervals
  ## Check the class
  ci_colnames <- colnames(x21$Path_estimates)[-c(1:6)]
  
  # Are there more confidence intervals than the default (the 95% percentile CI)
  # Inform the user to use xxx instead.
  if(length(ci_colnames) > 2) {
    cat2(
      "By default, only one confidence interval supplied to `.ci` is printed.\n",
      "Use `xxx` to print all confidence intervals (not yet implemented)."
    )
    ci_colnames <- ci_colnames[1:2]
    cat("\n\n")
  }
  
  ## Path estimates
  cat2("Estimated path coefficients:\n============================")
  printSummarizePath(x, .ci_colnames = ci_colnames)
  
  ## Loadings and Weights
  printSummarizeLoadingsWeights(x, .ci_colnames = ci_colnames)
  

  if(.full_output && x22$Model$model_type == "Linear") {
    ### Effects ----------------------------------------------------------------
    cat2("\n\n", rule(center = "Effects", width = 80), "\n\n")
    ## Path estimates
    cat2("Estimated total effects:\n========================")
    
    printSummarizePath(x, .ci_colnames = ci_colnames, .what = "Total effect")
    
    cat2("\n\nEstimated indirect effects:\n===========================")
    
    printSummarizePath(x, .ci_colnames = ci_colnames, .what = "Indirect effect")
  }

  cat2("\n", rule(line = "bar2", width = 80))
}

#' `cSEMVerify_default` method for `print()`
#'
#' The `cSEMVerify_default` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
#' @keywords internal
print.cSEMVerify_default <- function(x, ...) {
  
  cat(rule(line = "bar2", width = 80), sep = "")
  
  cat("\n\nVerify admissibility:\n", sep = "")
  
  if(sum(x) == 0) {
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
  for(i in names(x)) {
    cat("  ", col_align(i, 7), 
        col_align(ifelse(x[i] == FALSE, green("ok"), red("not ok")), 10), 
        col_align(text[i], max(nchar(text)) + 2, align = "left"), "\n", sep = "")
  }
  cat(rule(line = "bar2", width = 80), sep = "")
}

#' `cSEMVerify_2ndorder` method for `print()`
#'
#' The `cSEMVerify_default` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
#' @keywords internal
print.cSEMVerify_2ndorder <- function(x, ...) {
  
  cat(rule(line = "bar2", width = 80), sep = "")
  cat("\n\nVerify admissibility:\n", sep = "")
  
  if(sum(x$First_stage) == 0 & sum(x$Second_stage) == 0) {
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
  for(i in names(x$First_stage)) {
    cat("  ", col_align(i, 7), 
        col_align(ifelse(x$First_stage[i] == FALSE, green("ok"), red("not ok")), 10), 
        col_align(ifelse(x$Second_stage[i] == FALSE, green("ok"), red("not ok")), 10), 
        col_align(text[i], max(nchar(text) + 2), align = "left"), "\n", sep = "")
  }
  cat(rule(line = "bar2", width = 80), sep = "")
}

#' `cSEMTestOMF` method for `print()`
#'
#' The `cSEMTestOMF` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
#' @keywords internal
print.cSEMTestOMF <- function(x, ...) {
  
  cat(
    rule(line = "bar2", width = 80), "\n",
    rule(center = "Test for overall model fit based on Beran & Srivastava (1985)",
         width = 80), 
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
  cat("\nAdditonal information:")
  cat(
    "\n\n\tOut of ", x$Information$Total_runs , " bootstrap replications ", 
    x$Information$Number_admissibles, " are admissible.\n\t",
    "See ", yellow("?"), magenta("verify"), "()",
    " for what constitutes an inadmissible result.\n", 
    sep = ""
  )
  cat(rule(line = "bar2", width = 80), sep = "")
}

#' `cSEMTestMGD` method for `print()`
#'
#' The `cSEMTestMGD` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
#' @keywords internal
print.cSEMTestMGD <- function(x, ...) {
  
  cat(
    rule(line = "bar2", width = 80), "\n",
    rule(center = "Test for multigroup differences based on Klesel (forthcoming)",
         width = 80), 
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
  cat(
    "\n\n\tOut of ", x$Information$Total_runs , " permutation runs, ", 
    x$Information$Number_admissibles, " where admissible.\n\t",
    "See ", yellow("?"), magenta("verify"), "()",
    " for what constitutes an inadmissible result.", 
    sep = ""
  )
  
  cat("\n\n\tNumber of observations per group:")
  
  l <- max(nchar(c(x$Information$Group_names, "Group")))
  
  cat("\n\n\t",
      col_align("Group", width = l + 6),
      col_align("No. observations", width = 15),
      sep = ""
  )
  for(i in seq_along(x$Information$Group_names)) {
    cat(
      "\n\t",
      col_align(x$Information$Group_names[i], width = l + 6), 
      col_align(x$Information$Number_of_observations[i], width = 15),
      sep = ""
    )
  }
  
  cat("\n", rule(line = "bar2", width = 80), sep = "")
}

#' `cSEMTestMICOM` method for `print()`
#'
#' The `cSEMTestMICOM` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [foreman()], [cSEMResults]
#'
#' @export
#' @keywords internal
print.cSEMTestMICOM <- function(x, ...) {
  
  cat(
    rule(line = "bar2", width = 80), "\n",
    rule(center = "Test for measurement invariance based on Henseler et al (2016)",
         width = 80),
    "\n",
    sep = ""
  )
  
  x2 <- x$Step2
  x3 <- x$Step3
  
  construct_names <- names(x2$Test_statistic[[1]])
  no_groups       <- length(x2$Test_statistic)
  l_names <- max(nchar(c("Construct", construct_names)))
  
  ### Step 1 ===================================================================
  cat(
    rule(center = "Step 1 - Configural invariance", line = 2, width = 80), "\n\n",
    "\tConfigural invariance is a precondition for step 2 and 3.\n",
    "\tDo not proceed to interpret results unless\n",
    "\tconfigural invariance has been established.\n\n",
    sep = ""
  )
  ### Step 2 ===================================================================
  cat(rule(center = "Step 2 - Compositional invariance", line = 2, width = 80))
  
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
  cat(rule(center = "Step 3 - Equality of the means and variances", line = 2, width = 80))
  
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
    "\n\n\tOut of ", x$Information$Total_runs , " permutation runs, ", 
    x$Information$Number_admissibles, " where admissible.\n\t",
    "See ", yellow("?"), magenta("verify"), "()",
    " for what constitutes an inadmissible result.", 
    sep = ""
  )
  
  cat("\n\n\tNumber of observations per group:")
  
  l <- max(nchar(c(x$Information$Group_names, "Group")))
  
  cat("\n\n\t",
      col_align("Group", width = l + 6),
      col_align("No. observations", width = 15),
      sep = ""
  )
  for(i in seq_along(x$Information$Group_names)) {
    cat(
      "\n\t",
      col_align(x$Information$Group_names[i], width = l + 6), 
      col_align(x$Information$Number_of_observations[i], width = 15),
      sep = ""
    )
  }
  
  cat("\n", rule(line = "bar2", width = 80), sep = "")
}

