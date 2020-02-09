#' `cSEMSurface` method for `plot()`
#'
#' Plot the predicted values of an independent variable (z) 
#' The values are predicted based on a certain moderator and a certain 
#' independent variable including all their higher-order terms.
#' 
#' @param x An R object of class `cSEMSurface`.
#' @param .plot_type A character vector indicating the plot package used.
#' @param ... Currently ignored.
#' @export
plot.cSEMSurface <- function(x,
                             .plot_type = c('plotly','rsm'),
                             ...) {
  if(!(.plot_type %in% c('plotly','rsm'))){
    stop2("Currenlty only the plotly and the rsm package are supported.")
  }
  
  if(.plot_type == 'plotly'){
    
  
  ## Install ggplot2 if not already installed
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop2(
      "Package `plotly` required. Use `install.packages(\"plotly\")` and rerun.")
  }

  plot1 <- plotly::plot_ly( x = x$out$x, y = x$out$z, z = x$out$y, type = "surface",...)
  plot1

  
  }

  if(.plot_type == 'rsm'){
  # see rsm package for other plots

  }
}