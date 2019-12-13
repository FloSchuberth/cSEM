#' `cSEMVerify` method for `print()`
#'
#' The `cSEMVerify` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cSEMResults], [verify()]
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