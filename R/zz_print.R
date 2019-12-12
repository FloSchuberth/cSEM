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


#' `cSEMAssess` method for `print()`
#'
#' The `cSEMAssess` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [assess()]
#'
#' @export
#' @keywords internal
print.cSEMAssess <- function(x, ...) {
  
  cat2(rule2(type = 2))
  nn <- intersect(names(x), c("AVE", "RhoC", "RhoC_weighted", "RhoT", "RhoT_weighted", "R2", "R2_adj"))
  
  if(length(nn) > 0) {
    # Max name length 
    c_names <- names(x[[nn[1]]])
    if(length(c_names) > 0) {
      l <- max(nchar(c_names)) 
      
      ## If more than 4 quality criteria are to be printed: open a second block
      ## otherwise print one block
      if(length(nn) > 4) {
        nn1 <- nn[1:4]
        nn2 <- setdiff(nn, nn1)
        nn <- list(nn1, nn2)
      } else {
        nn <- list(nn)
      }
      
      for(j in seq_along(nn)) {
        cat2(
          "\n\n\t", 
          col_align("Construct", max(l, nchar("Construct")) + 2)
        )
        for(k in nn[[j]]) {
          cat2(
            col_align(k, 14, align = "center")
          )
        }
        
        for(i in c_names) {
          cat2(
            "\n\t", 
            col_align(i, max(l, nchar("Construct")) + 2)
          )
          
          for(k in nn[[j]]) {
            cat2(col_align(sprintf("%.4f",x[[k]][i]), 14, align = "center"))
          }
        } 
      }
    }
  }
  
  if(any(names(x) %in% c("CFI", "GFI", "IFI", "NFI", "NNFI", "RMSEA", 
                         "RMS_theta", "SRMR", "DG", "DL", "DML"))) {
    cat2("\n\n", rule2("Distance and fit measures"), "\n")
    
    if(any(names(x) == "DG")) {
      cat2(col_align("\n\tGeodesic distance", 30), "= ", x$DG)
    }
    if(any(names(x) == "DL")) {
      cat2(col_align("\n\tSquared Euclidian distance", 30), "= ", x$DL)
    }
    if(any(names(x) == "DML")) {
      cat2(col_align("\n\tML distance", 30), "= ", x$DML)
    }
    
    cat2("\n")
    
    if(any(names(x) == "CFI")) {
      cat2(col_align("\n\tCFI", 15), "= ", x$CFI)
    }
    if(any(names(x) == "GFI")) {
      cat2(col_align("\n\tGFI", 15), "= ", x$GFI)
    }
    if(any(names(x) == "IFI")) {
      cat2(col_align("\n\tIFI", 15), "= ", x$IFI)
    }
    if(any(names(x) == "NFI")) {
      cat2(col_align("\n\tNFI", 15), "= ", x$NFI)
    }
    if(any(names(x) == "NNFI")) {
      cat2(col_align("\n\tNNFI", 15), "= ", x$NNFI)
    }
    if(any(names(x) == "RMSEA")) {
      cat2(col_align("\n\tRMSEA", 15), "= ", x$RMSEA)
    }
    if(any(names(x) == "RMS_theta")) {
      cat2(col_align("\n\tRMS_theta", 15), "= ", x$RMS_theta)
    }
    if(any(names(x) == "SRMR")) {
      cat2(col_align("\n\tSRMR", 15), "= ", x$SRMR)
    }
    
    if(any(names(x) == "Df")) {
      cat2(col_align("\n\n\tDegrees of freedom", 25), "= ", x$Df)
    }
  }
  
  if(any(names(x) == "VIF")) {
    cat2("\n\n", rule2("Variance inflation factors (VIFs)"))
    for(i in rownames(x$VIF)) {
      cat2("\n\n  Dependent construct: '", i, "'\n")
      cat2(
        "\n\t", 
        col_align("Independent construct", max(l, nchar("Independent construct")) + 2), 
        col_align("VIF value", 12, align = "center")
      )
      for(j in names(x$VIF[i, ])) {
        cat2(
          "\n\t", 
          col_align(j, max(l, nchar("Independent construct")) + 2), 
          col_align(sprintf("%.4f", x$VIF[i, j]), 12, align = "center")
        )  
      }
    }
  }
  
  if(any(names(x) == "Effect_size")) {
    cat2("\n\n", rule2("Effect sizes (f_squared)"))
    for(i in rownames(x$Effect_size)) {
      cat2("\n\n  Dependent construct: '", i, "'\n")
      cat2(
        "\n\t", 
        col_align("Independent construct", max(l, nchar("Independent construct")) + 2), 
        col_align("Effect size", 12, align = "center")
      )
      for(j in colnames(x$Effect_size[i, x$Effect_size[i, ] != 0, drop = FALSE])) {
        cat2(
          "\n\t", 
          col_align(j, max(l, nchar("Independent construct")) + 2), 
          col_align(sprintf("%.4f", x$Effect_size[i, j]), 12, align = "center")
        )  
      }
    }
  }
  
  if(any(names(x) %in% c("Fornell-Larcker", "HTMT", "RA"))) {
    cat2("\n\n", rule2("Validity assessment"))
    if(any(names(x) == "HTMT") && !anyNA(x$HTMT)) {
      cat2("\n\n\tHeterotrait-montrait ratio of correlation matrix (HTMT matrix)\n\n")
      print(x$HTMT)
    }
    
    if(any(names(x) == "Fornell-Larcker")) {
      cat2("\n\n\tFornell-Larcker matrix\n\n")
      print(x$`Fornell-Larcker`)
    }
    
    if(any(names(x) == "RA") && !anyNA(x$RA)) {
      cat2("\n\n\tRedundancy analysis")
      cat2(
        "\n\n\t", 
        col_align("Construct", max(l, nchar("Construct")) + 2),
        col_align("Value", 14, align = "center")
      )
      for(i in names(x$RA)) {
        cat2(
          "\n\t", 
          col_align(i, max(l, nchar("Construct")) + 2),
          col_align(sprintf("%.4f",x$RA[i]), 14, align = "center")
        ) 
      }
    }
  }
  
  cat2("\n", rule2(type = 2))
}


#' `cSEMPredict` method for `print()`
#'
#' The `cSEMPredict` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [predict()]
#'
#' @export
#' @keywords internal
print.cSEMPredict <- function(x, ...) {
  
  x1 <- x$Prediction_metrics
  x2 <- x$Information
  
  cat2(
    rule2(type = 2), "\n",
    rule2("Overview"), 
    "\n"
  )
  
  # cat2("\n\tGeneral information:\n\t","------------------------")
  cat2(
    col_align("\n\tNumber of obs. training", 35), "= ", x2$Number_of_observations_training,
    col_align("\n\tNumber of obs. test", 35), "= ", x2$Number_of_observations_test,
    col_align("\n\tNumber of cv folds", 35), "= ", x2$Number_of_folds,
    col_align("\n\tNumber of repetitions", 35), "= ", x2$Number_of_repetitions,
    col_align("\n\tHandle inadmissibles", 35), "= ", x2$Handle_inadmissibles,
    col_align("\n\tBenchmark", 35), "= ", paste0("'", x2$Benchmark, "'")
    )
  
  ### Prediction metricts-------------------------------------------------------
  cat2("\n\n", rule2("Prediction metrics"), "\n\n")

  l <- max(nchar(x1[, "Name"]))
  
  cat2(
    "\n  ", 
    col_align("Name", l + 2), 
    col_align("MAE target", 13, align = "right"), 
    col_align("MAE benchmark", 15, align = "right"), 
    col_align("RMSE target", 13, align = "right"),
    col_align("RMSE benchmark", 15, align = "right"),
    col_align("Q2_predict", 13, align = "right")
  )
  
  for(i in 1:nrow(x1)) {
    cat2(
      "\n  ", 
      col_align(x1[i, "Name"], l + 2), 
      col_align(sprintf("%.4f", x1[i, "MAE_target"]), 13, align = "right"),
      col_align(sprintf("%.4f", x1[i, "MAE_benchmark"]), 15, align = "right"),
      col_align(sprintf("%.4f", x1[i, "RMSE_target"]), 13, align = "right"),
      col_align(sprintf("%.4f", x1[i, "RMSE_benchmark"]), 15, align = "right"),
      col_align(sprintf("%.4f", x1[i, "Q2_predict"]), 13, align = "right")
    )
  }

  
  cat2("\n", rule2(type = 2))
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
    
    # Correlations
    construct_cor <- x21$Exo_construct_correlation
    res_cor  <- x11$Residual_correlation
    indi_cor <- x11$Indicator_correlation 
  } else {
    
    x21 <- x$Estimates
    x22 <- x$Information
    
    # Correlation
    construct_cor <- x21$Exo_construct_correlation
    res_cor  <- x21$Residual_correlation
    indi_cor <- x21$Indicator_correlation 
  }
  
  cat2(
    rule2(type = 2), "\n",
    rule2("Overview"), 
    "\n"
  )
  
  ### Overview -----------------------------------------------------------------
  ## General information + resample information
  printSummarizeOverview(x)
  
  ## Construct details
  cat2("\n\n\tConstruct details:\n\t","------------------")
  
  printSummarizeConstructDetails(x)
  
  ### Estimates ----------------------------------------------------------------
  cat2("\n\n", rule2("Estimates"), "\n\n")
  
  ## Confidence intervals
  # Get the column names of the columns containing confidence intervals
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
  printSummarizePathCorrelation(x, .ci_colnames = ci_colnames)
  
  ## Loadings and Weights
  printSummarizeLoadingsWeights(x, .ci_colnames = ci_colnames)
  
  ## Exogenous construct correlation
  if(.full_output && nrow(construct_cor) != 0) {
    cat2("\n\nEstimated construct correlations:\n=================================")
    printSummarizePathCorrelation(x, .ci_colnames = ci_colnames, 
                                  .what = "Construct correlation")
  }
  
  ## Residual correlation
  if(.full_output && nrow(res_cor) != 0) {
    cat2("\n\nEstimated measurement error correlations:\n=========================================")
    printSummarizePathCorrelation(x, .ci_colnames = ci_colnames, 
                                  .what = "Residual correlation")
  }
  
  ## Indicator correlation
  if(.full_output && nrow(indi_cor) != 0) {
    cat2("\n\nEstimated indicator correlations:\n=================================")
    printSummarizePathCorrelation(x, .ci_colnames = ci_colnames, 
                                  .what = "Indicator correlation")
  }

  if(.full_output && x22$Model$model_type == "Linear") {
    ### Effects ----------------------------------------------------------------
    cat2("\n\n", rule2("Effects"), "\n\n")
    ## Path estimates
    cat2("Estimated total effects:\n========================")
    
    printSummarizePathCorrelation(x, .ci_colnames = ci_colnames, 
                                  .what = "Total effect")
    
    cat2("\n\nEstimated indirect effects:\n===========================")
    
    printSummarizePathCorrelation(x, .ci_colnames = ci_colnames, 
                                  .what = "Indirect effect")
  }

  cat2("\n", rule2(type = 2))
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
print.cSEMVerify <- function(x, ...) {
  
  cat2(rule2(type = 2))
  
  cat2("\n\nVerify admissibility:\n")
  
  if(inherits(x, "cSEMVerify_multi")) {
    l <- max(nchar(names(x)), nchar("Dataset")) + 2
    
    cat2(
      "\n\t", 
      col_align("Dataset", l, align = "left"),
      "Status"
      )

    for(j in seq_along(x)) {
      n_defects <- sum(sapply(x[[j]], sum))
      
      cat2(
        "\n\t", 
        col_align(names(x)[j], l, align = "left"),
        ifelse(n_defects == 0, green("admissible"), red("inadmissible"))
        )
    }
  } else {
    n_defects <- sum(sapply(x, sum))
    
    if(n_defects == 0) {
      cat2(green("\n\t admissible"))
    } else {
      cat2(red("\n\t inadmissible"))
    }
  }
  
  cat2("\n\nDetails:\n")
  
  if(inherits(x, "cSEMVerify_multi")) {
    for(j in seq_along(x)) {
      cat2("\n", names(x)[j], "\n\n")
      printVerifyDetails(x[[j]])
    } 
  } else {
    cat2("\n")
    printVerifyDetails(x)
  }
  
  cat2(rule2(type = 2))
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
  
  cat2(
    rule(line = "bar2", width = 80), "\n",
    rule(center = "Test for overall model fit based on Beran & Srivastava (1985)",
         width = 80)
  )
  
  ## Null hypothesis -----------------------------------------------------------
  cat2(
    "\n\nNull hypothesis:\n\n", 
    boxx(c("H0: Population indicator covariance matrix is equal to", 
           "model-implied indicator covariance matrix."), float = "center")
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
  cat2(
    "\n\n\tOut of ", x$Information$Total_runs , " bootstrap replications ", 
    x$Information$Number_admissibles, " are admissible.\n\t",
    "See ", yellow("?"), magenta("verify"), "()",
    " for what constitutes an inadmissible result.\n\n\t",
    "The seed used was: ", x$Information$Seed, "\n"
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

  info <- x$Information
  approaches <- c("Klesel", "Sarstedt", "Chin", "Keil", "Nitzl", "Henseler")
  
  if(any(info$Approach == "all")) {
    info$Approach <- approaches
  } else {
    info$Approach <- info$Approach[match(info$Approach, intersect(approaches, info$Approach))]
  }
  ## Additional information ----------------------------------------------------
  cat2(
    rule2(type = 2), "\n",
    rule2("Overview")
  )
  cat2(
    col_align("\n\n\tTotal runs (permutation)", width = 37), "= ", info$Information_permutation$Total_runs,
    col_align("\n\tAdmissible results (permutation)", width = 36), "= ", info$Information_permutation$Number_admissibles,
    col_align("\n\tPermutation seed", width = 36), "= ", info$Information_permutation$Permutation_seed,
    "\n\n\tNumber of observations per group:"
  )
  
  l <- max(nchar(c(info$Group_names, "Group")))
  
  cat2("\n\n\t",
       col_align("Group", width = l + 6),
       col_align("No. Obs.", width = 10, align = "center")
  )
  for(i in seq_along(info$Group_names)) {
    cat2(
      "\n\t",
      col_align(info$Group_names[i], width = l + 6), 
      col_align(info$Number_of_observations[i], width = 10, align = "center")
    )
  }
  
  ## Overall descision only for Sarstedt, Chin and Keil
  approach <- intersect(info$Approach, c("Sarstedt", "Chin", "Keil", "Nitzl", "Henseler"))
  if(length(approach) > 0) {
    cat2("\n\n\tOverall decision (based on alpha = ", paste0(info$Alpha[1] * 100, "%):"))
    cat2("\n\n\t",col_align("", width = 10))
    for(j in names(x[[approach[1]]]$Decision_overall)) {
      cat2(
        col_align(paste0("p_adjust = '", j, "'"), width = 20, align = "right")
      )
    }
    for(i in approach) {
      cat2(
        "\n\t",
        col_align(i, width = 10)
      )
      for(j in seq_along(x[[i]]$Decision_overall)) {
        cat2(
          col_align(ifelse(x[[i]]$Decision_overall[[j]][[1]], green("Do not reject"), red("reject")),
                    width = 20, align = "right")
        )
      }
    }
  }
  
  cat2("\n")
  
  # If alpha contains more than one element, inform the user that only one alpha 
  # is printed
  if(length(info$Alpha)  > 1) {
    cat2(
      "\n\tNote: Due to space constraits of the console only results for\n ",
      "\t      alpha = ", paste0(info$Alpha[1] * 100, "%", " are shown.")
      )
    cat2("\n")
  }
  
  ## Klesel et al. (2019) ======================================================
  if(any(info$Approach == "Klesel")) {
    xk <- x$Klesel
    
    cat2(
      rule2(type = 2), "\n",
      rule2("Test for multigroup differences based on Klesel et al. (2019)")
    )
    
    ## Null hypothesis ---------------------------------------------------------
    cat2(
      "\n\nNull hypothesis:\n\n", 
      boxx(paste0("H0: Model-implied ", xk$VCV_type, " covariance matrix is equal across groups."),
           float = "center")
    )
    
    ## Test statistic and p-value ----------------------------------------------
    cat2("\n\nTest statistic and p-value: \n\n\t")
    # Are several .alphas given? Inform the user that only the first .alpha is
    # is used for decision
    if(length(info$Alpha) > 1) {
      cat2(
        "Decision is based on alpha = ", names(xk$Decision)[1],
        "\n\n\t")
    }
    
    cat2(
      col_align("Distance measure", width = 20),
      col_align("Test statistic", width = 14, align = "right"), 
      col_align("p-value", width = 16, align = "right"),
      col_align("Decision", width = 16, align = "right")
    )
    
    for(j in seq_along(xk$Test_statistic)) {
      
      cat2(
        "\n\t",
        col_align(names(xk$Test_statistic)[j], width = 20),
        col_align(sprintf("%.4f", xk$Test_statistic[j]), width = 14, 
                  align = "right"), 
        col_align(sprintf("%.4f", xk$P_value[j]), width = 16, align = "right"),
        col_align(ifelse(xk$Decision[[1]][j], green("Do not reject"), red("reject")),
                  width = 16, align = "right")
      )
    }
    cat2("\n")
  }
  
  ## Sarstedt et al. (2011) ====================================================
  if(any(info$Approach == "Sarstedt")) {
    xs <- x$Sarstedt
    
    cat2(
      rule2(type = 2), "\n",
      rule2("Test for multigroup differences based on Sarstedt et al. (2011)")
    )
    
    cat2(
      red("\n\n\t!WARNING: the test is unreliable. Usage is discouraged!")
    )
    ## Null hypothesis ---------------------------------------------------------
    cat2(
      "\n\nNull hypothesis:\n\n",
      boxx("H0: Parameter k is equal across all groups.", float = "center")
    )
    
    ## Test statistic and p-value ----------------------------------------------
    cat2("\n\nTest statistic and p-value: \n\n")
    
    # Are several .alphas given? Inform the user that only the first .alpha is
    # is used for decision
    if(length(info$Alpha) > 1) {
      cat2("\tDecision is based on the alpha = ", names(xs$Decision[[1]])[1])
    }
    
    l <- max(10, nchar(names(xs$Test_statistic)))
    
    # Create table for every p-value adjustment method
    for(i in seq_along(xs$P_value)) {
      
      cat2("\n\tMultiple testing adjustment: '", names(xs$P_value)[i], "'\n\n\t")
      
      cat2(
        col_align("Parameter", width = l + 4),
        col_align("Test statistic", width = 14, align = "right"), 
        col_align("p-value", width = 16, align = "right"),
        col_align("Decision", width = 16, align = "right")
      )
      
      for(j in seq_along(xs$Test_statistic)) {
        
        cat2(
          "\n\t",
          col_align(names(xs$Test_statistic)[j], width = l + 4),
          col_align(sprintf("%.4f", xs$Test_statistic[j]), width = 14, 
                    align = "right"), 
          col_align(sprintf("%.4f", xs$P_value[[i]][j]), width = 16, align = "right"),
          col_align(ifelse(xs$Decision[[i]][[1]][j], green("Do not reject"), red("reject")),
                    width = 16, align = "right")
        )
      }
      cat2("\n")
    }
  }
  
  ## Chin & Dibbern (2010) =====================================================
  if(any(info$Approach == "Chin")) {
    printTestMGDResults(.x = x, .approach = "Chin", .info = info)
  }
  
  ## Keil et. al (2000)=========================================================
  if(any(info$Approach == "Keil")) {
    printTestMGDResults(.x = x, .approach = "Keil", .info = info)
  }
  
  ## Nitzl (2010) ==============================================================
  if(any(info$Approach == "Nitzl")) {
    printTestMGDResults(.x = x, .approach = "Nitzl", .info = info)
  }
  
  ## Henseler (2009) ===========================================================
  if(any(info$Approach == "Henseler")) {
    printTestMGDResults(.x = x, .approach = "Henseler", .info = info)
  }
  cat2(rule2(type = 2))
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
      x3_type <- x3[[type]]
      cat(type, "\n\t", sep = "")
      cat(
        col_align("", width = l_names),
        "\t",
        col_align("", width = 14), 
        "\t",
        col_align("Critical values", width = 8*ncol(x3_type$Critical_value[[i]]), 
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
      
      for(j in colnames(x3_type$Critical_value[[i]])) {
        cat(col_align(j, width = 6, align = "center"), "\t", sep = "")
      }
      
      cat("\n\t")
      
      for(j in construct_names) {
        cat(
          col_align(j, width = l_names), 
          "\t",
          col_align(sprintf("%.4f", x3_type$Test_statistic[[i]][j]), width = 14,
                    align = "center"),
          "\t",
          sep = ""
        )
        
        for(k in colnames(x3_type$Critical_value[[i]])) {
          cat(sprintf("%.4f", x3_type$Critical_value[[i]][j, k]), "\t", sep = "")
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
      x3_type <- x3[[type]]
      
      # Width of columns
      l <- apply(x3_type$Decision[[i]], 2, function(x) {
        ifelse(any(x == TRUE), nchar("Do not reject"), nchar("reject"))
      })
      
      n1 <- colnames(x3_type$Critical_value[[i]])
      n2 <- paste0("[", n1[seq(1, (length(n1) - 1), by = 2)], ";", 
                   n1[seq(2, length(n1), by = 2)], "]")
      
      l1 <- max(c(sum(l) + 3*(ncol(x3_type$Decision[[i]]) - 1)), 
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
      
      for(j in seq_along(x3_type$Test_statistic[[i]])) {
        
        cat(col_align(names(x3_type$Test_statistic[[i]])[j], width = l_names), "\t", sep = "")
        
        for(k in 1:ncol(x3_type$Decision[[i]])) {
          cat(
            col_align(ifelse(x3_type$Decision[[i]][j, k], 
                             green("Do not reject"), red("reject")),
                      width = max(l[k], nchar(n2[k])), align = "center"), 
            "\t", 
            sep = ""
          )
        } # END for k (each construct)
        cat("\n\t")
      } # END for j (each significance level)
    } # END for type (one of "Mean" and "Var")
    cat("\n")
  } # END for i (each group)
  
  ## Additional information ----------------------------------------------------
  cat("Additional information:")
  cat(
    "\n\n\tOut of ", x$Information$Total_runs , " permutation runs, ", 
    x$Information$Number_admissibles, " where admissible.\n\t",
    "See ", yellow("?"), magenta("verify"), "()",
    " for what constitutes an inadmissible result.\n\n\t",
    "The seed used was: ", x$Information$Seed, "\n", 
    sep = ""
  )
  
  cat("\n\tNumber of observations per group:")
  
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

