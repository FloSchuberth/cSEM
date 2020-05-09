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
  
  x <- neffects
  
  info  <- x$Information_print
  dep   <- x$Information$dependent
  indep <- x$Information$independent
  mod   <- x$Information$moderator
  # Johnson Neyman-point
  # Es kann folgende Faelle geben kein, ein oder zwei Punkt,
  # da muss man sich was überlegen
  x$out_floodlight$Johnson_Neyman_points
  
  
  # Slope when moderator = mean(=0)-2SD
  cat2("\n\nSlope when", mod, " = ? (-2 SD)")
  x$Information_print[x$Information_print[,'value_z']==-2,]['direct_effect']
  # CI
  # lb:
  x$Information_print[x$Information_print[,'value_z']==-2,]['lb']
  # ub:
  x$Information_print[x$Information_print[,'value_z']==-2,]['ub']
  
  # Slope when moderator = mean(=0)-1SD
  x$Information_print[x$Information_print[,'value_z']==-1,]['direct_effect']
  # CI
  # lb:
  x$Information_print[x$Information_print[,'value_z']==-1,]['lb']
  # ub:
  x$Information_print[x$Information_print[,'value_z']==-1,]['ub']
  
  # Slope when moderator at mean (=0)
  x$Information_print[x$Information_print[,'value_z']==0,]['direct_effect']
  # CI
  # lb:
  x$Information_print[x$Information_print[,'value_z']==0,]['lb']
  # ub:
  x$Information_print[x$Information_print[,'value_z']==0,]['ub']
  
  # Slope when moderator = mean(=0)+1SD
  x$Information_print[x$Information_print[,'value_z']==1,]['direct_effect']
  # CI
  # lb:
  x$Information_print[x$Information_print[,'value_z']==1,]['lb']
  # ub:
  x$Information_print[x$Information_print[,'value_z']==1,]['ub']
  
  # Slope when moderator = mean(=0)+2SD
  x$Information_print[x$Information_print[,'value_z']==2,]['direct_effect']
  # CI
  # lb:
  x$Information_print[x$Information_print[,'value_z']==2,]['lb']
  # ub:
  x$Information_print[x$Information_print[,'value_z']==2,]['ub']
  cat2("\n", rule2(type = 2), "\n")
}  