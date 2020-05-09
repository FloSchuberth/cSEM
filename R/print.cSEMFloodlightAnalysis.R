#' `cSEMFloodlight` method for `print()`
#'
#' The `cSEMFloodlight` method for the generic function [print()]. 
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [cSEMResults], [doFloodlightAnalysis()], [plot.cSEMFloodlight()]
#'
#' @export
#' @keywords internal
print.cSEMFloodlight <- function(x, ...) {
 
  if(!inherits(x, "cSEMFloodlight")) {
    stop2("x must be of class `cSEMFloodlight`.")
  }
  
  # Johnson Neyman-point
  # Es kann folgende Faelle geben kein, ein oder zwei Punkt,
  # da muss man sich was überlegen
 x$Johnson_Neyman_points
  
   
  # Slope when moderator = mean(=0)-2SD
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

  
  }  