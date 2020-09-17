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
  
  ## AVE, R2, R2_adj -----------------------------------------------------------
  nn <- intersect(names(x), c("AVE", "R2", "R2_adj"))
  
  if(length(nn) > 0) {
    # Select all construct names that are either endogenous (LHS) variables (for 
    # the R2 and the R2_adj) or constructs for which the AVE was computed. 
    c_names <- unique(c(names(x[["AVE"]]), names(x[["R2"]]), names(x[["R2_adj"]])))
    if(length(c_names) > 0) {
      l <- max(nchar(c_names)) 
      nn <- list(nn)
      
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
  
  ## Reliability ---------------------------------------------------------------
  
  if(any(names(x) == "Reliability")) {
    rel <- x$Reliability
    c_names <- names(rel[[1]])
    if(length(c_names) > 0) {
      l <- max(nchar(c_names))
      cat2("\n\n", rule2("Common (internal consistency) reliability estimates"), "\n")
      
      cat2(
        "\n\t", 
        col_align("Construct", l + 2, align = "left"), 
        col_align(names(rel)[1], 17, align = "center"), 
        col_align(names(rel)[2], 18, align = "center"),
        col_align(names(rel)[3], 26, align = "center")
      )
      
      for(i in c_names) {
        cat2(
          "\n\t",
          col_align(i, l + 2, align = "left"),
          col_align(sprintf("%.4f", rel[[1]][i]), 17, align = "center"),
          col_align(sprintf("%.4f", rel[[2]][i]), 18, align = "center"),
          col_align(sprintf("%.4f", rel[[3]][i]), 26, align = "center")
        )
        
      } 
    }
  }
  
  ## Alternative reliability measures ------------------------------------------
  
  nn <- intersect(names(x), c("RhoC", "RhoC_mm", "RhoC_weighted", "RhoC_weighted_mm",
                              "RhoT", "RhoT_weighted"))
  
  if(length(nn) > 0) {
    # Max name length 
    c_names <- names(x[[nn[1]]])
    if(length(c_names) > 0) {
      cat2("\n\n", rule2("Alternative (internal consistency) reliability estimates"))
      l <- max(nchar(c_names)) 
      
      ## If more than 4 quality criteria are to be printed: open a second block
      ## otherwise print one block
      if(length(nn) > 3) {
        nn1 <- nn[1:3]
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
  
  ## Distance and fit measures -------------------------------------------------
  
  if(any(names(x) %in% c("Chi_square", "Chi_square_df", "CFI", "CN", "GFI", "IFI", 
                         "NFI", "NNFI", "RMSEA", "Df",
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
    
    if(any(names(x) == "Chi_square")) {
      cat2(col_align("\n\tChi_square", 17), "= ", x$Chi_square)
    }
    if(any(names(x) == "Chi_square_df")) {
      cat2(col_align("\n\tChi_square_df", 17), "= ", x$Chi_square_df)
    }
    if(any(names(x) == "CFI")) {
      cat2(col_align("\n\tCFI", 17), "= ", x$CFI)
    }
    if(any(names(x) == "CN")) {
      cat2(col_align("\n\tCN", 17), "= ")
      cat(x$CN, sep = ", ")
    }
    if(any(names(x) == "GFI")) {
      cat2(col_align("\n\tGFI", 17), "= ", x$GFI)
    }
    if(any(names(x) == "IFI")) {
      cat2(col_align("\n\tIFI", 17), "= ", x$IFI)
    }
    if(any(names(x) == "NFI")) {
      cat2(col_align("\n\tNFI", 17), "= ", x$NFI)
    }
    if(any(names(x) == "NNFI")) {
      cat2(col_align("\n\tNNFI", 17), "= ", x$NNFI)
    }
    if(any(names(x) == "RMSEA")) {
      cat2(col_align("\n\tRMSEA", 17), "= ", x$RMSEA)
    }
    if(any(names(x) == "RMS_theta")) {
      cat2(col_align("\n\tRMS_theta", 17), "= ", x$RMS_theta)
    }
    if(any(names(x) == "RMS_theta_mi")) {
      cat2(col_align("\n\tRMS_theta_mi", 17), "= ", x$RMS_theta)
    }
    if(any(names(x) == "SRMR")) {
      cat2(col_align("\n\tSRMR", 17), "= ", x$SRMR)
    }
    
    if(any(names(x) == "Df")) {
      cat2(col_align("\n\n\tDegrees of freedom", 25), "= ", x$Df)
    }
  }
  
  ## Model selection criteria --------------------------------------------------
  nn <- intersect(names(x), c("AIC", "AICc", "AICu", "BIC", "FPE", "GM", "HQ", 
                              "HQc", "Mallows_Cp"))
  if(length(nn) > 0) {
    # Get names of the endogenous variables
    c_names <- names(x[[nn[1]]])
    if(length(c_names) > 0) {
      cat2("\n\n", rule2("Model selection criteria"))
      l <- max(nchar(c_names)) 
      
      ## If more than 4 selection criteria are to be printed: open a second/third
      ## block otherwise print one block
      if(length(nn) > 3) {
        nn1 <- nn[1:3]
        nn2 <- setdiff(nn, nn1)
        if(length(nn2) > 3) {
          nn2 <- nn2[1:3]
          nn3 <- setdiff(nn, c(nn1, nn2))
        }
        nn <- list(nn1, nn2, nn3)
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

  ## Variance inflation factors ------------------------------------------------
  
  if(any(names(x) == "VIF")) {
    cat2("\n\n", rule2("Variance inflation factors (VIFs)"))
    for(i in rownames(x$VIF)) {
      cat2("\n\n  Dependent construct: '", i, "'\n")
      cat2(
        "\n\t", 
        col_align("Independent construct", max(nchar(names(x$VIF[i, ])), nchar("Independent construct")) + 2), 
        col_align("VIF value", 12, align = "center")
      )
      for(j in names(which(x$VIF[i, ] != 0))) {
        cat2(
          "\n\t", 
          col_align(j, max(nchar(names(x$VIF[i, ])), nchar("Independent construct")) + 2), 
          col_align(sprintf("%.4f", x$VIF[i, j]), 12, align = "center")
        )  
      }
    }
  }
  
  ## Variance inflation factors for modeB constructs ---------------------------
  
  if(any(names(x) == "VIF_modeB") && !is.na(x$VIF_modeB)) {
    cat2("\n\n", rule2("Variance inflation factors (VIFs) for modeB constructs"))
    for(i in rownames(x$VIF_modeB)) {
      cat2("\n\n  Construct: '", i, "'\n")
      cat2(
        "\n\t", 
        col_align("Weight", max(nchar(names(x$VIF_modeB[i, ])), nchar("Weight")) + 2), 
        col_align("VIF value", 12, align = "center")
      )
      for(j in names(which(x$VIF_modeB[i, ] != 0))) {
        cat2(
          "\n\t", 
          col_align(j, max(nchar(names(x$VIF_modeB[i, ])), nchar("Weight")) + 2), 
          col_align(sprintf("%.4f", x$VIF_modeB[i, j]), 12, align = "center")
        )  
      }
    }
  }
  
  ## Effect size analysis ------------------------------------------------------
  
  if(any(names(x) == "F2")) {
    cat2("\n\n", rule2("Effect sizes (Cohen's f^2)"))
    for(i in rownames(x$F2)) {
      ll <- nchar(colnames(x$F2[i, x$F2[i, ] != 0, drop = FALSE]))
      cat2("\n\n  Dependent construct: '", i, "'\n")
      cat2(
        "\n\t", 
        col_align("Independent construct", max(ll, nchar("Independent construct")) + 2), 
        col_align("f^2", 12, align = "center")
      )
      for(j in colnames(x$F2[i, x$F2[i, ] != 0, drop = FALSE])) {
        cat2(
          "\n\t", 
          col_align(j, max(ll, nchar("Independent construct")) + 2), 
          col_align(sprintf("%.4f", x$F2[i, j]), 12, align = "center")
        )  
      }
    }
  }
  
  ## Validity assessment -------------------------------------------------------
  
  if(any(names(x) %in% c("Fornell-Larcker", "HTMT"))) {
    cat2("\n\n", rule2("Validity assessment"))
    if(any(names(x) == "HTMT") && !is.null(x$HTMT)) {
      cat2("\n\n\tHeterotrait-monotrait ratio of correlations matrix (HTMT matrix)\n\n")
      if(x$Information$.inference) {
        cat2("\tValues in the upper triangular part are the ", 
             paste0(100*(1 - x$Information$.alpha), "%-quantile of the\n", 
            "\tbootstrap distribution (using .ci = '", x$Information$.ci, "').\n\n")) 
      }
      print(x$HTMT)
    }
    
    if(any(names(x) == "Fornell-Larcker")) {
      cat2("\n\n\tFornell-Larcker matrix\n\n")
      print(x$`Fornell-Larcker`)
    }
  }
  
  ## Effects -------------------------------------------------------------------
  
  if(any(names(x) == "Effects")) {
    cat2("\n\n", rule2("Effects"), "\n\n")
    
    ## Confidence intervals
    # Get the column names of the columns containing confidence intervals
    ci_colnames <- colnames(x$Effects$Total_effects)[-c(1:6)]
    
    # Are there more confidence intervals than the default (the 95% percentile CI)
    # Inform the user to use xxx instead.
    if(length(ci_colnames) > 2) {
      cat2(
        "By default, only one confidence interval is printed."
      )
      ci_colnames <- ci_colnames[1:2]
      cat("\n\n")
    }
    
    if(any(names(x$Effects) == "Total_effect")) {
      ## Total effects
      cat2("Estimated total effects:\n========================")
      printEffects(x$Effects$Total_effect, .ci_colnames = ci_colnames, .what = "Total effect")
    }
    if(any(names(x$Effects) == "Indirect_effect")) {
      ## Indirect effects
      cat2("\n\nEstimated indirect effects:\n===========================")
      printEffects(x$Effects$Indirect_effect, .ci_colnames = ci_colnames, .what = "Indirect effect")
    }
    if(any(names(x$effects) == "Variance_accounted_for")) {
      ### Variance accounted for -------------------------------------------------
      cat2("\n\nVariance accounted for (VAF):\n=============================")
      printEffects(x$Effects$Variance_accounted_for, .ci_colnames = ci_colnames, .what = "Effects")
    }
  }
  
  cat2("\n", rule2(type = 2), "\n")
}