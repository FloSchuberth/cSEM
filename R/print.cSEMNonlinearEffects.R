#' `cSEMNonlinearEffectsAnalysis` method for `print()`
#'
#' The `cSEMNonlinearEffectsAnalysis` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cSEMResults], [doNonlinearEffectsAnalysis()], [plot.cSEMNonlinearEffects()]
#'
#' @export
#' @keywords internal
print.cSEMNonlinearEffects <- function(x, ...) {
  
  if(!inherits(x, "cSEMNonlinearEffects")) {
    stop2("x must be of class `cSEMNonlinearEffects`")
  }
  cat2(
    rule2(type = 2), "\n",
    rule2("Simple slope analysis", type = 1, align = "center")
  )

  info  <- x$Information_print
  dep   <- x$Information$dependent
  indep <- x$Information$independent
  mod   <- x$Information$moderator
  JN    <- x$out_floodlight$Johnson_Neyman_points
  
  cat2("\n\nFloodlight analysis is based on a value of ", 
       x$Information$value_independent , " for the independent variable `", indep, "`.")
  
  for(i in 1:nrow(info)) {
    cat2("\n\nSlope of `", indep, "` when `", mod, "` is at ", 
         switch (i,
                 "1" = {paste0(info[i, 2], " SDs from mean")},
                 "2" = {paste0(info[i, 2], " SD from mean")},
                 "3" = {"the mean (= 0)"},
                 "4" = {paste0("+", info[i, 2], " SD from mean")},
                 "5" = {paste0("+", info[i, 2], " SDs from mean")}
         ), "\n")
    
    cat2(
      "\n  ", 
      col_align("Effect", 16, align = "right"), 
      col_align(colnames(info)[3], 12, align = "right"),
      col_align(colnames(info)[4], 12, align = "right"),
      "\n  ", 
      col_align(sprintf("%.4f", info[i, 1]), 16, align = "right"),
      col_align(sprintf("%.4f", info[i, 3]), 12, align = "right"),
      col_align(sprintf("%.4f", info[i, 4]), 12, align = "right")
    )
  }
  
  ### Johnson-Neyman points ----------------------------------------------------
  # (18.05.2020) is JN always a matrix?; make sure it is
  if(nrow(JN) > 0) {
    cat2("\n\n", rule2("Johnson-Neyman points"))
    
    cat2(
      "\n  ", 
      col_align("x", 12, align = "right"),
      col_align("y", 12, align = "right"))
    
    for(i in 1:nrow(JN)) {
      cat2(
        "\n  ",
        col_align(sprintf("%.4f", JN[i, 1]), 12, align = "right"),
        col_align(sprintf("%.4f", JN[i, 2]), 12, align = "right")
      )
    } 
  } 
  
  cat2("\n", rule2(type = 2), "\n")
}  