#' `cSEMIPA` method for `plot()`
#'
#' Plot the importance-performance matrix
#' 
#' @param x An R object of class `cSEMSurface`.
#' @param .plot_type A character vector indicating the plot package used. Options are
#' "*plotly*", "*rsm*", and "*persp*". Defaults to "*plotly*". 
#' @param ... Additional parameters that can be passed to 
#' \code{\link[graphics:persp]{graphics::persp}}, e.g., to rotate the plot.
#' 
#' @seealso [doIPA()]
#' @export

plot.cSEMIPA <- function(x,.dependent,.level){
  # select relevant variables
  rel_imp_value<-x$Importance[.dependent,,drop=FALSE]
  
  # Remove 0s from the row
  rel_imp_value<-rel_imp_value[rel_imp_value!=0]
  
  # data.frame for dots to be plotted
  data_plot
  
}
