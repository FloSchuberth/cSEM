#' Helper for print.summarize_default and print.summerize_2ndorder
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
    col_align("\n\tNumber of observations", 35), "= ", nrow(x$Arguments$.data),
    col_align("\n\tWeight estimator", 35), "= ", 
    ifelse(x$Arguments$.approach_weights == "PLS-PM" && 
             x$Type_of_indicator_correlation %in% c("Polychoric", "Polyserial"), 
           "PLS-PM (OrdPLS)", x$Arguments$.approach_weights)
  )
  
  if(x$Arguments$.approach_weights == "PLS-PM") {
    cat2(
      col_align("\n\tInner weighting scheme", 35), "= ", 
      x$Arguments$.PLS_weight_scheme_inner
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
    col_align("\n\tType of path model", 35), "= ", x$Model$model_type,
    col_align("\n\tDisattenuated", 35), "= ", 
    ifelse(x$Arguments$.disattenuate & any(x$Model$construct_type == "Common factor"), 
           ifelse(x$Arguments$.approach_weights == "PLS-PM", "Yes (PLSc)",
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
      col_align("\n\tResample methode", 35), "= ", xx$Information_resample$Method,
      col_align("\n\tNumber of resamples", 35), "= ", xx$Information_resample$Number_of_runs
    )
    if(xx$Information_resample$Method2 %in% c("bootstrap", "jackknife")) {
      cat2(
        col_align("\n\tResample of resample methode", 35), "= ", xx$Information_resample$Method2,
        col_align("\n\tNumber of resamples per resample", 35), "= ", xx$Information_resample$Number_of_runs2
      ) 
    }
    cat2(
      col_align("\n\tNumber of admissible results ", 35), "= ", xx$Information_resample$Number_of_admissibles,
      col_align("\n\tApproach to handle inadmissibles ", 35), "= ", xx$Information_resample$Handle_inadmissibles,
      col_align("\n\tSign change option", 35), "= ", xx$Information_resample$Sign_change_option
    )
    if(!isFALSE(xx$Information_resample$Seed)) {
      cat2(
        col_align("\n\tRandom seed", 35), "= ", xx$Information_resample$Seed
      )
    }
  }
}

#' Helper for print.summarize_default and print.summerize_2ndorder
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
    c_2nd_order           <- grep("_temp", rownames(x2$Model$structural), 
                                  value = TRUE, invert = TRUE)
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
    cat2(col_align("Mode", 5))
  }
  cat("\n")
  
  if(inherits(.summarize_object, "cSEMSummarize_2ndorder")) {
    
    # First stage
    for(i in names(x1$Model$construct_type)) {
      cat("\n\t", 
          col_align(i, max(l, nchar("Name")) + 2), 
          col_align(x1$Model$construct_type[i], 13 + 2), 
          col_align("First order", 12 + 2), sep = "")
      if(x1$Arguments$.approach_weights == "PLS-PM") {
        cat(col_align(x1$Weight_info$Modes[i], 5), sep = "")
      }
    }
    # Second stage
    for(i in names(x2$Model$construct_type[c_2nd_order])) {
      cat("\n\t", 
          col_align(i, max(l, nchar("Name")) + 2), 
          col_align(x2$Model$construct_type[i], 13 + 2), 
          col_align("Second order", 12 + 2), sep = "")
      if(x2$Arguments$.approach_weights == "PLS-PM") {
        cat(col_align(x2$Weight_info$Modes[i], 5), sep = "")
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
        cat2(col_align(x2$Weight_info$Modes[i], 5))
      }
    }
  }
}

#' Helper for print.summarize_default and print.summerize_2ndorder
#' @noRd
#' 
printSummarizePath <- function(.summarize_object, .ci_colnames, .what = "Path") {
  
  ## Check the class
  x <- if(inherits(.summarize_object, "cSEMSummarize_2ndorder")) {
    switch (.what,
      "Path" = {x <- .summarize_object$Second_stage$Estimates$Path_estimates},
      "Total effect" = {.summarize_object$Second_stage$Estimates$Effect_estimates$Total_effect},
      "Indirect effect" = {.summarize_object$Second_stage$Estimates$Effect_estimates$Indirect_effect}
    )

  } else {
    switch (.what,
      "Path" = {.summarize_object$Estimates$Path_estimates},
      "Total effect" = {.summarize_object$Estimates$Effect_estimates$Total_effect},
      "Indirect effect" = {.summarize_object$Estimates$Effect_estimates$Indirect_effect}
    )
  }
  
  l <- max(nchar(x[, "Name"]))
  
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
    col_align(.what, max(l, nchar(.what)) + 2), 
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
                   sprintf("%7.4f", x[i, j+1]), "]"), 20, align = "center")
        )
      } 
    }
  }
}

#' Helper for print.summarize_default and print.summerize_2ndorder
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
                     sprintf("%7.4f", x[i, j+1]), "]"), 20, align = "center")
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

  cat2("\n\nEstimated Loadings:\n===================")
  
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

  cat2("\n\nEstimated Weights:\n==================")
  
  if(length(.ci_colnames) != 0) {
    cat2("\n  ",  col_align("", width = max(l, nchar("Weights")) + 44))
    for(i in interval_names) {
      cat2(col_align(i, width = 20*length(sig_level_names), align = "center"))
    }
  }
  
  cat2(
    "\n  ", 
    col_align("Weights", max(l, nchar("Weights")) + 2), 
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