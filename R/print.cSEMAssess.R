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
  nn <- intersect(names(x), c("AVE", "R2", "R2_adj"))
  
  if(length(nn) > 0) {
    # Max name length 
    c_names <- names(x[[nn[1]]])
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
  
  if(any(names(x) == "reliability")) {
    rel <- x$reliability
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
  
  
  if(any(names(x) %in% c("Chi_square", "Chi_square_df", "CFI", "GFI", "IFI", 
                         "NFI", "NNFI", "RMSEA", "Df",
                         "RMS_theta", "RMS_theta_mi", "SRMR", "DG", "DL", "DML"))) {
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
  
  if(any(names(x) == "VIF")) {
    cat2("\n\n", rule2("Variance inflation factors (VIFs)"))
    for(i in rownames(x$VIF)) {
      cat2("\n\n  Dependent construct: '", i, "'\n")
      cat2(
        "\n\t", 
        col_align("Independent construct", max(nchar(names(x$VIF[i, ])), nchar("Independent construct")) + 2), 
        col_align("VIF value", 12, align = "center")
      )
      for(j in names(x$VIF[i, ])) {
        cat2(
          "\n\t", 
          col_align(j, max(nchar(names(x$VIF[i, ])), nchar("Independent construct")) + 2), 
          col_align(sprintf("%.4f", x$VIF[i, j]), 12, align = "center")
        )  
      }
    }
  }
  
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
  
  if(any(names(x) %in% c("Fornell-Larcker", "HTMT"))) {
    cat2("\n\n", rule2("Validity assessment"))
    if(any(names(x) == "HTMT") && !is.null(x$HTMT)) {
      cat2("\n\n\tHeterotrait-monotrait ratio of correlations matrix (HTMT matrix)\n\n")
      if(x$Information$.inference) {
        cat2("Values in the upper triangular part are the ", 
             paste0(100*x$Information$.p, "%-Quantile of the bootstrap distribution.\n\n")) 
        x$HTMT <- round(x$HTMT, 4)
        diag(x$HTMT) <- "."
      } else {
        x$HTMT <- round(x$HTMT, 4)
        x$HTMT[upper.tri(x$HTMT, diag = TRUE)] <- "."
      }
      print(x$HTMT)
    }
    
    if(any(names(x) == "Fornell-Larcker")) {
      cat2("\n\n\tFornell-Larcker matrix\n\n")
      print(round(x$`Fornell-Larcker`, 4))
    }
  }
  
  if(any(names(x) == "Effects")) {
    ### Effects ----------------------------------------------------------------
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
    
    ## Total effects
    cat2("Estimated total effects:\n========================")
    printEffects(x$Effects$Total_effect, .ci_colnames = ci_colnames, .what = "Total effect")
    
    ## Indirect effects
    cat2("\n\nEstimated indirect effects:\n===========================")
    printEffects(x$Effects$Indirect_effect, .ci_colnames = ci_colnames, .what = "Indirect effect")
    
    ### Variance accounted for -------------------------------------------------
    cat2("\n\nVariance accounted for (VAF):\n=============================")
    printEffects(x$Effects$Variance_accounted_for, .ci_colnames = ci_colnames, .what = "Effects")
  }
  
  cat2("\n", rule2(type = 2))
}