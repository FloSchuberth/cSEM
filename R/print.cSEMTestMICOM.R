#' `cSEMTestMICOM` method for `print()`
#'
#' The `cSEMTestMICOM` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cSEMResults], [testMICOM()]
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