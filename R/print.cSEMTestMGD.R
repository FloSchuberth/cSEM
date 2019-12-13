#' `cSEMTestMGD` method for `print()`
#'
#' The `cSEMTestMGD` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cSEMResults], [testMGD()]
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