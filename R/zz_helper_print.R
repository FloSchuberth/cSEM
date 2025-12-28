#' Helper for print.cSEMSummarize
#' @noRd
#' 
printSummarizeOverview <- function(.summarize_object) {
  
  ## Check the class
  x <- if(inherits(.summarize_object, "cSEMSummarize_2ndorder")) {
    .summarize_object$Second_stage$Information
  } else {
    .summarize_object$Information
  }
  
  cat2("\n\tGeneral information:\n\t","------------------------")
  cat2(
    col_align("\n\tEstimation status", 35), "= ", ifelse(sum(x$Estimation_status) == 0, green("Ok"), 
                                                         c(red("Not ok!"), "Use: verify()."))
  )
  cat2(
    col_align("\n\tNumber of observations", 35), "= ", nrow(x$Arguments$.data),
    col_align("\n\tWeight estimator", 35), "= ", 
    ifelse(x$Arguments$.approach_weights == "PLS-PM" && 
             !all(x$Type_of_indicator_correlation == 'Pearson'), 
           "PLS-PM (OrdPLS)", x$Arguments$.approach_weights)
  )
  
  if(x$Arguments$.approach_weights == "PLS-PM") {
    cat2(
      col_align("\n\tInner weighting scheme", 35), "= ", '"', 
      x$Arguments$.PLS_weight_scheme_inner, '"'
    )
  }
  cat2(
    col_align("\n\tType of indicator correlation", 35), "= ", 
    paste0(x$Type_of_indicator_correlation, collapse = ", ")
  )
  
  if(!is.null(x$Model$instruments) && x$Arguments$.approach_paths == "2SLS") {
    tmp <- setdiff(names(which(rowSums(x$Model$structural) != 0)), names(x$Model$instruments))
    cat2(
      "\n\tPath model estimator",
      "\n\t\tOLS  : ", paste0(tmp, collapse = ", "),
      "\n\t\t2SLS : ", paste0(names(x$Model$instruments), collapse = ", ")
    )
  } else {
    cat2(
      col_align("\n\tPath model estimator", 35), "= ", x$Arguments$.approach_paths
    )
  }

  cat2(
    col_align("\n\tSecond-order approach", 35), "= ", x$Approach_2ndorder
  )
  
  
  # 
  
  
  cat2(
    col_align("\n\tType of path model", 35), "= ", x$Model$model_type,
    col_align("\n\tDisattenuated", 35), "= ", 
    ifelse(x$Arguments$.disattenuate & any(x$Model$construct_type == "Common factor"), 
           ifelse(x$Arguments$.approach_weights == "PLS-PM", 
                  # If the model contains second-order constructs, the mixed or two-stage approach 
                  # is applied. Hence, the model is estimated in two stages and 
                  # consequently .disattenuate is reported for both stages.
                  # If further approaches are implemented that require only one stage 
                  # such as the repeated indicators approach we need to do a more detailed check here
                  if(!is.na(x$Approach_2ndorder)){
                    paste("First stage:", ifelse(.summarize_object$First_stage$Information$Arguments$.disattenuate==TRUE,
                                                 "Yes",
                                                 "No"), "\n",col_align("= Second stage:", 55, "right"), 
                          ifelse(x$Arguments$.disattenuate == TRUE, "Yes", "No"))
                  }else{
                         "Yes (PLSc)"
           },
                  ifelse(x$Arguments$.approach_weights == "GSCA", "Yes (GSCAm)", "Yes")
           ), "No")
  )
  
  ## Check the class
  xx <- if(inherits(.summarize_object, "cSEMSummarize_2ndorder")) {
    x$Resample
  } else {
    x
  }
  
  ## Resample information
  if(inherits(.summarize_object, "cSEMSummarize_resampled")) {
    cat2("\n\n\tResample information:\n\t","---------------------")
    cat2(
      col_align("\n\tResample method", 35), "= ", '"', xx$Information_resample$Method, '"',
      col_align("\n\tNumber of resamples", 35), "= ", xx$Information_resample$Number_of_runs
    )
    if(xx$Information_resample$Method2 %in% c("bootstrap", "jackknife")) {
      cat2(
        col_align("\n\tResample of resample method", 35), "= ", '"',xx$Information_resample$Method2, '"',
        col_align("\n\tNumber of resamples per resample", 35), "= ", xx$Information_resample$Number_of_runs2
      ) 
    }
    cat2(
      col_align("\n\tNumber of admissible results ", 35), "= ", xx$Information_resample$Number_of_admissibles,
      col_align("\n\tApproach to handle inadmissibles ", 35), "= ", '"', xx$Information_resample$Handle_inadmissibles, '"', 
      col_align("\n\tSign change option", 35), "= ", '"', xx$Information_resample$Sign_change_option, '"'
    )
    if(!isFALSE(xx$Information_resample$Seed)) {
      cat2(
        col_align("\n\tRandom seed", 35), "= ", xx$Information_resample$Seed
      )
    }
  }
}

#' Helper for print.cSEMSummarize
#' @noRd
#' 
printSummarizeConstructDetails <- function(.summarize_object) {
  
  ## Check the class
  if(inherits(.summarize_object, "cSEMSummarize_2ndorder")) {
    
    x1 <- .summarize_object$First_stage$Information
    x2 <- .summarize_object$Second_stage$Information
    
    ## Collect all necessary sets
    # All constructs used in the first step (= all first order constructs)
    c_linear_1step        <- rownames(x1$Model$structural)
    # All second order constructs
    # c_2nd_order           <- grep("_temp", rownames(x2$Model$structural), 
    #                               value = TRUE, invert = TRUE)
    c_2nd_order <- names(x1$Model$construct_order[x1$Model$construct_order == "Second order"])
    # # All linear constructs of the original model
    c_linear_original     <- c(c_linear_1step, c_2nd_order)
    
    # Max name length for all constructs
    l <- max(nchar(c_linear_original))
    
  } else {
    x2 <- .summarize_object$Information
    
    # Max name length for all constructs
    l <- max(nchar(names(x2$Model$construct_type)))
  }
  
  cat2(
    "\n\t", 
    col_align("Name", max(l, nchar("Name")) + 2), 
    col_align("Modeled as", 13 + 2),
    col_align("Order", 12 + 2)
  )
  if(x2$Arguments$.approach_weights == "PLS-PM") {
    cat2(col_align("Mode", 10))
  }
  cat("\n")
  
  if(inherits(.summarize_object, "cSEMSummarize_2ndorder")) {
    
    # First stage
    for(i in names(x1$Model$construct_type)) {
      cat2("\n\t", 
          col_align(i, max(l, nchar("Name")) + 2), 
          col_align(x1$Model$construct_type[i], 13 + 2), 
          col_align("First order", 12 + 2))
      if(x1$Arguments$.approach_weights == "PLS-PM") {
        cat2(col_align(paste0('"', x1$Weight_info$Modes[i],'"'), 10))
      }
    }
    # Second stage
    for(i in names(x2$Model$construct_type[c_2nd_order])) {
      cat2("\n\t", 
          col_align(i, max(l, nchar("Name")) + 2), 
          col_align(x2$Model$construct_type[i], 13 + 2), 
          col_align("Second order", 12 + 2))
      if(x2$Arguments$.approach_weights == "PLS-PM") {
        cat2(col_align(paste0('"', x2$Weight_info$Modes[i],'"'), 10))
      }
    }
  } else {
    for(i in names(x2$Model$construct_type)) {
      cat2(
        "\n\t", 
        col_align(i, max(l, nchar("Name")) + 2), 
        col_align(x2$Model$construct_type[i], 13 + 2), 
        col_align(x2$Model$construct_order[i], 12 + 2)
      )
      if(x2$Arguments$.approach_weights == "PLS-PM") {
        cat2(col_align(paste0('"', x2$Weight_info$Modes[i],'"'), 10))
      }
    }
  }
}

#' Helper for print.cSEMSummarize
#' @noRd
#' 
printSummarizePathCorrelation <- function(.summarize_object, .ci_colnames, .what = "Path") {
  
  ## Check the class
  x <- if(inherits(.summarize_object, "cSEMSummarize_2ndorder")) {
    switch (.what,
      "Path" = {x <- .summarize_object$Second_stage$Estimates$Path_estimates},
      "Total effect" = {.summarize_object$Second_stage$Estimates$Effect_estimates$Total_effect},
      "Indirect effect" = {.summarize_object$Second_stage$Estimates$Effect_estimates$Indirect_effect},
      "Residual correlation" = {.summarize_object$First_stage$Estimates$Residual_correlation},
      "Indicator correlation" = {.summarize_object$First_stage$Estimates$Indicator_correlation},
      "Construct correlation" = {.summarize_object$Second_stage$Estimates$Exo_construct_correlation}
    )

  } else {
    switch (.what,
      "Path" = {.summarize_object$Estimates$Path_estimates},
      "Total effect" = {.summarize_object$Estimates$Effect_estimates$Total_effect},
      "Indirect effect" = {.summarize_object$Estimates$Effect_estimates$Indirect_effect},
      "Residual correlation" = {.summarize_object$Estimates$Residual_correlation},
      "Indicator correlation" = {.summarize_object$Estimates$Indicator_correlation},
      "Construct correlation" = {.summarize_object$Estimates$Exo_construct_correlation}
    )
  }
  
  # Rename .what
  .what <- ifelse(.what %in% c("Residual correlation", "Indicator correlation", "Construct correlation"), 
         "Correlation", .what)
  
  l <- max(nchar(x[, "Name"]), nchar(.what))
  
  if(length(.ci_colnames) != 0) {
    xx <- regmatches(.ci_colnames, regexpr("\\.", .ci_colnames), invert = TRUE)
    interval_names    <- unique(sapply(xx, `[`, 1))
    sig_level_names   <- unique(gsub("[LU]", "", sapply(xx, `[`, 2)))
    
    cat2("\n  ",  col_align("", width = max(l, nchar(.what)) + 44))
    for(i in interval_names) {
      cat2(col_align(i, width = 20*length(sig_level_names), align = "center"))
    }
  }
  cat2(
    "\n  ", 
    col_align(.what, l + 2), 
    col_align("Estimate", 10, align = "right"), 
    col_align("Std. error", 12, align = "right"),
    col_align("t-stat.", 10, align = "right"), 
    col_align("p-value", 10, align = "right")
  )
  if(length(.ci_colnames) != 0) {
    for(i in rep(sig_level_names, length(interval_names))) {
      cat2(col_align(i, 20, align = "center"))
    } 
  }
  
  for(i in 1:nrow(x)) {
    cat2(
      "\n  ", 
      col_align(x[i, "Name"], l + 2), 
      col_align(sprintf("%.4f", x[i, "Estimate"]), 10, align = "right"),
      col_align(sprintf("%.4f", x[i, "Std_err"]), 12, align = "right"),
      col_align(sprintf("%.4f", x[i, "t_stat"]), 10, align = "right"),
      col_align(sprintf("%.4f", x[i, "p_value"]), 10, align = "right")
    )
    if(length(.ci_colnames) != 0) {
      for(j in seq(1, length(.ci_colnames), by = 2) + 
          ifelse(.what == "Correlation", 5, 6)) {
        cat2(
          col_align(
            paste0("[", sprintf("%7.4f", x[i, j]), ";", 
                   sprintf("%7.4f", x[i, j+1]), " ]"), 20, align = "center")
        )
      } 
    }
  }
}

#' Helper for print.cSEMSummarize
#' @noRd
#' 
printSummarizeLoadingsWeights <- function(.summarize_object, .ci_colnames) {
  
  printLoadingsWeights <- function(x, .what) {
  
    x <- if(.what == "Loadings") {
      x$Loading_estimates
    } else {
      x$Weight_estimates
    }
    
    for(i in 1:nrow(x)) {
      # Note (19.05.2021): infer() also computes standard deviations and confidence
      # intervals for constant values. Because of floating point impressions
      # its possible that the sd is not exactly zero in some instances. This
      # messes up the print method. In this case, set to NA
      if(isTRUE(all.equal(x[i, "Std_err"], 0))) {
        x[i, 4:ncol(x)] <- NA
      }
      cat2(
        "\n  ", 
        col_align(x[i, "Name"], max(l, nchar(.what)) + 2), 
        col_align(sprintf("%.4f", x[i, "Estimate"]), 10, align = "right"), 
        col_align(sprintf("%.4f", x[i, "Std_err"]), 12, align = "right"),
        col_align(sprintf("%.4f", x[i, "t_stat"]), 10, align = "right"),
        col_align(sprintf("%.4f", x[i, "p_value"]), 10, align = "right")
      )
      if(length(.ci_colnames) != 0) {
        for(j in seq(1, length(.ci_colnames), by = 2) + 6) {
          cat2(
            col_align(
              paste0("[", sprintf("%7.4f", x[i, j]), ";", 
                     sprintf("%7.4f", x[i, j+1]), " ]"), 20, align = "center")
          )
        } 
      }
    }
  }
  
  ## Check the class
  if(inherits(.summarize_object, "cSEMSummarize_2ndorder")) {
    
    x1 <- .summarize_object$First_stage$Estimates
    x2 <- .summarize_object$Second_stage$Estimates
    
    # Max name length for all constructs
    l <- max(nchar(c(x1$Loading_estimates[, "Name"], x2$Loading_estimates[, "Name"])))
  } else {
    x2 <- .summarize_object$Estimates
    
    # Max name length for all constructs
    l <- max(nchar(x2$Loading_estimates[, "Name"]))
  }
  
  if(length(.ci_colnames) != 0) {
    xx <- regmatches(.ci_colnames, regexpr("\\.", .ci_colnames), invert = TRUE)
    interval_names    <- unique(sapply(xx, `[`, 1))
    sig_level_names   <- unique(gsub("[LU]", "", sapply(xx, `[`, 2)))
    }

  cat2("\n\nEstimated loadings:\n===================")
  
  if(length(.ci_colnames) != 0) {
    cat2("\n  ",  col_align("", width = max(l, nchar("Loadings")) + 44))
    for(i in interval_names) {
      cat2(col_align(i, width = 20*length(sig_level_names), align = "center"))
    }
  }
    
  cat2(
    "\n  ", 
    col_align("Loading", max(l, nchar("Loading")) + 2), 
    col_align("Estimate", 10, align = "right"), 
    col_align("Std. error", 12, align = "right"),
    col_align("t-stat.", 10, align = "right"), 
    col_align("p-value", 10, align = "right")
  )
  
  if(length(.ci_colnames) != 0) {
    for(i in rep(sig_level_names, length(interval_names))) {
      cat2(col_align(i, 20, align = "center"))
    } 
  }
  
  if(inherits(.summarize_object, "cSEMSummarize_2ndorder")) {
    printLoadingsWeights(x1, "Loadings")
  } 
  printLoadingsWeights(x2, "Loadings")

  cat2("\n\nEstimated weights:\n==================")
  
  if(length(.ci_colnames) != 0) {
    cat2("\n  ",  col_align("", width = max(l, nchar("Weights")) + 44))
    for(i in interval_names) {
      cat2(col_align(i, width = 20*length(sig_level_names), align = "center"))
    }
  }
  
  cat2(
    "\n  ", 
    col_align("Weight", max(l, nchar("Weights")) + 2), 
    col_align("Estimate", 10, align = "right"), 
    col_align("Std. error", 12, align = "right"),
    col_align("t-stat.", 10, align = "right"), 
    col_align("p-value", 10, align = "right")
  ) 
  
  if(length(.ci_colnames) != 0) {
    for(i in rep(sig_level_names, length(interval_names))) {
      cat2(col_align(i, 20, align = "center"))
    } 
  }
  
  if(inherits(.summarize_object, "cSEMSummarize_2ndorder")) {
    # ## Only print weights, for constructs modeled as composites
    # x1$Weight_estimates <- x1$Weight_estimates[x1$Weight_estimates$Construct_type == "Composite", , drop = FALSE] 
    # 
    printLoadingsWeights(x1, "Weights")
  } 
  # x2$Weight_estimates <- x2$Weight_estimates[x2$Weight_estimates$Construct_type == "Composite", , drop = FALSE] 
  printLoadingsWeights(x2, "Weights")
}

#' Helper for print.cSEMVerify
#' @noRd
#'
printVerifyDetails <- function(.x) {
  
  text <- c("1" = "Convergence achieved", 
            "2" = "All absolute standardized loading estimates <= 1", 
            "3" = "Construct VCV is positive semi-definite", 
            "4" = "All reliability estimates <= 1",
            "5" = "Model-implied indicator VCV is positive semi-definite")
  
  if(inherits(.x, "cSEMVerify_2ndorder")) {
    cat2(
      "  ", 
      col_align("Code", 7), 
      col_align("Stage 1", 10), 
      col_align("Stage 2", 10), 
      "Description\n"
      )
    
    for(i in names(.x$First_stage)) {
      cat2("  ", col_align(i, 7), 
           col_align(ifelse(.x$First_stage[i] == FALSE, green("ok"), red("not ok")), 10), 
           col_align(ifelse(.x$Second_stage[i] == FALSE, green("ok"), red("not ok")), 10), 
           col_align(text[i], max(nchar(text) + 2), align = "left"), "\n")
    }
  } else {
    cat2("  ", 
         col_align("Code", 7), 
         col_align("Status", 10), 
         "Description\n")
    
    for(i in names(.x)) {
      cat2("  ", col_align(i, 7), 
           col_align(ifelse(.x[i] == FALSE, green("ok"), red("not ok")), 10), 
           col_align(text[i], max(nchar(text)) + 2, align = "left"), "\n")
    }
  }
}

#' Helper for print.cSEMTestMGD
#' @noRd
#'
printTestMGDResults <- function(.x, .approach, .info) {
  
  switch(.approach,
    "Chin"  = {
      x <- .x$Chin
      
      cat2(
        rule2(type = 2), "\n",
        rule2("Test for multigroup differences based on Chin & Dibbern (2010)")
      )
      },
    "Keil"  = {
      x <- .x$Keil
      
      cat2(
        rule2(type = 2), "\n",
        rule2("Test for multigroup differences based on Keil et al. (2000)")
      )
      },
    "Nitzl" = {
      x <- .x$Nitzl
      cat2(
        rule2(type = 2), "\n",
        rule2("Test for multigroup differences based on Nitzl (2010)")
      )
    },
    "Henseler" = {
      x <- .x$Henseler
      cat2(
        rule2(type = 2), "\n",
        rule2("Test for multigroup differences based on Henseler (2007)")
      )
    },
    "CI_para" = {
      x <- .x$CI_para
      cat2(
        rule2(type = 2), "\n",
        rule2("Test for multigroup differences: CI_para")
      )
    },
    "CI_overlap" = {
      x <- .x$CI_overlap
      cat2(
        rule2(type = 2), "\n",
        rule2("Test for multigroup differences: CI_overlap")
      )
    }
  )
  
  ## Null hypothesis -----------------------------------------------------------
  # if(.approach == "Henseler") {
  #   cat2(
  #     "\n\nNull hypothesis:\n\n",
  #     boxx(c("(1) H0: Parameter k of group 1 is smaller than that of group 2.", 
  #          "(2) H0: Parameter k of group 1 is larger than that of group 2."),
  #          float = "center")
  #   )
  # } else {
  #   cat2(
  #     "\n\nNull hypothesis:\n\n",
  #     boxx("H0: Parameter k is equal across two groups.", float = "center")
  #   )
  # }
  
  cat2(
    "\n\nNull hypothesis:\n\n",
    boxx("H0: Parameter k is equal across two groups.", float = "center",width=80)
  )

  
  ## Test statistic and p-value ------------------------------------------------
  cat2("\n\nTest statistic and p-value: \n\n")
  # Are several .alphas given? Inform the user that only the first .alpha is
  # is used for decision
  
  if(.approach == "CI_para") {
    # Are several .alphas given? Inform the user that only the first .alpha is
    # is used for decision
    if(length(.info$Alpha) > 1) {
      cat2(
        "\tDecision is based on the ", names(x$Decision)[1], " confidence interval",
        "\n\tType of confidence interval: ", names(x$Decision[[1]][[1]])[1]
      )
    }
    
    l <- max(10, nchar(x$Decision[[1]][[1]][[1]]$Name))
    
    # Print
    for(i in seq_along(x$Decision[[1]])) {
      
      cat2("\n\n  Compared groups: ", names(x$Decision[[1]])[i], "\n\n")
      
      cat2(
        "  ",
        col_align("Parameter", width = l + 2),
        col_align("Est. group 1", width = 14, align = "right"), 
        col_align("CI group 2", width = 20, align = "center"),
        col_align("Est. group 2", width = 14, align = "right"), 
        col_align("CI group 1", width = 20, align = "center")
        # col_align("Decision", width = 10, align = "right")
      )
      xx1 <- x$Decision[[1]][[i]][[1]]
      
      for(j in 1:nrow(xx1)) {
        cat2(
          "\n  ", 
          col_align(xx1[j, "Name"], l + 2), 
          col_align(sprintf("%.4f", xx1[j, 2]), 14, align = "right"),
          col_align(
            paste0("[", sprintf("%7.4f", xx1[j, 3]), ";", 
                   sprintf("%7.4f", xx1[j, 4]), " ]"), 20, align = "center"),
          col_align(sprintf("%.4f", xx1[j, 5]), 14, align = "right"),
          col_align(
            paste0("[", sprintf("%7.4f", xx1[j, 6]), ";", 
                   sprintf("%7.4f", xx1[j, 7]), " ]"), 20, align = "center")
          # col_align(ifelse(xx1[j, "Decision"], green("Do not reject"), red("reject")), 10, 
          #                  align = "right")
        )
      }
    }
    cat2("\n")
  } else if(.approach == "CI_overlap") {
    # Are several .alphas given? Inform the user that only the first .alpha is
    # is used for decision
    if(length(.info$Alpha) > 1) {
      cat2(
        "\tDecision is based on the ", names(x$Decision)[1], " confidence interval",
        "\n\tType of confidence interval: ", names(x$Decision[[1]][[1]])[1]
      )
    }
    
    l <- max(10, nchar(x$Decision[[1]][[1]][[1]]$Name))
    
    # Print
    for(i in seq_along(x$Decision[[1]])) {
      
      cat2("\n\n  Compared groups: ", names(x$Decision[[1]])[i], "\n\n")
      
      cat2(
        "  ",
        col_align("Parameter", width = l + 2),
        col_align("CI group 1", width = 20, align = "center"),
        col_align("CI group 2", width = 20, align = "center"),
        col_align("Decision", width = 14, align = "right")
      )
      xx1 <- x$Decision[[1]][[i]][[1]]
      
      for(j in 1:nrow(xx1)) {
        cat2(
          "\n  ", 
          col_align(xx1[j, "Name"], l + 2), 
          col_align(
            paste0("[", sprintf("%7.4f", xx1[j, 2]), ";", 
                   sprintf("%7.4f", xx1[j, 3]), " ]"), 20, align = "center"),
          col_align(
            paste0("[", sprintf("%7.4f", xx1[j, 4]), ";", 
                   sprintf("%7.4f", xx1[j, 5]), " ]"), 20, align = "center"),
          col_align(ifelse(xx1[j, "Decision"], green("Do not reject"), red("reject")), 14,
                           align = "right")
        )
      }
    }
    cat2("\n")
  } else {
    # Are several .alphas given? Inform the user that only the first .alpha is
    # is used for decision
    if(length(.info$Alpha) > 1) {
      cat2(
        "\tDecision is based on alpha = ", names(x$Decision[[1]])[1]
      )
    }
    
    l <- max(10, nchar(names(x$Test_statistic[[1]])))
    
    # Create table for every p-value adjustment method
    for(p in seq_along(x$P_value)) {
      cat2("\n\tMultiple testing adjustment: '", names(x$P_value)[p], "'")
      for(i in seq_along(x$Test_statistic)) {
        
        cat2("\n\n  Compared groups: ", names(x$Test_statistic)[i], "\n\n\t")
        
        cat2(
          col_align("Parameter", width = l),
          col_align("Test statistic", width = 14, align = "right"), 
          col_align("p-value", width = 16, align = "right"),
          col_align("Decision", width = 16, align = "right")
        )
        
        for(j in seq_along(x$Test_statistic[[i]])) {
          
          cat2(
            "\n\t",
            col_align(names(x$Test_statistic[[i]])[j], width = l),
            col_align(sprintf("%.4f", x$Test_statistic[[i]][j]), width = 14, 
                      align = "right"), 
            col_align(sprintf("%.4f", x$P_value[[p]][[i]][j]), width = 16, align = "right"),
            col_align(ifelse(x$Decision[[p]][[1]][[i]][j], green("Do not reject"), red("reject")),
                      width = 16, align = "right")
          )
        }
      } 
      cat2("\n")
    }
  }
}