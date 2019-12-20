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
  nn <- intersect(names(x), c("RhoC", "RhoC_weighted", "RhoT", "RhoT_weighted",
                              "R2", "R2_adj", "reliability"))
  
  if(any(names(x) %in% c("CFI", "GFI", "IFI", "NFI", "NNFI", "RMSEA", 
                         "RMS_theta", "SRMR", "DG", "DL", "DML"))) {
    cat2("\n\n", rule2("Reliability"), "\n")
    # # Reliability
    # nn <- c("AVE", "RhoC", "RhoC_weighted", "RhoT", "RhoT_weighted", 
    #   "R2", "R2_adj", "reliability")
    # ## If more than 4 quality criteria are to be printed: open a second block
    # ## otherwise print one block
    # if(length(nn) > 4) {
    #   nn1 <- nn[1:4]
    #   nn2 <- setdiff(nn, nn1)
    #   nn <- list(nn1, nn2)
    # } else {
    #   nn <- list(nn)
    # }
    
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
    cat2("\n\n", rule2("Effect sizes (Cohen's f^2)"))
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
  
  if(any(names(x) %in% c("Fornell-Larcker", "HTMT"))) {
    cat2("\n\n", rule2("Validity assessment"))
    if(any(names(x) == "HTMT") && !anyNA(x$HTMT)) {
      cat2("\n\n\tHeterotrait-montrait ratio of correlation matrix (HTMT matrix)\n\n")
      print(x$HTMT)
    }
    
    if(any(names(x) == "Fornell-Larcker")) {
      cat2("\n\n\tFornell-Larcker matrix\n\n")
      print(x$`Fornell-Larcker`)
    }
    
    # if(any(names(x) == "RA") && !anyNA(x$RA)) {
    #   cat2("\n\n\tRedundancy analysis")
    #   cat2(
    #     "\n\n\t", 
    #     col_align("Construct", max(l, nchar("Construct")) + 2),
    #     col_align("Value", 14, align = "center")
    #   )
    #   for(i in names(x$RA)) {
    #     cat2(
    #       "\n\t", 
    #       col_align(i, max(l, nchar("Construct")) + 2),
    #       col_align(sprintf("%.4f",x$RA[i]), 14, align = "center")
    #     ) 
    #   }
    # }
  }
  
  cat2("\n", rule2(type = 2))
}