#' `cSEMSurface` method for `plot()`
#'
#' Plot the predicted values of an independent variable (z) 
#' The values are predicted based on a certain moderator and a certain 
#' independent variable including all their higher-order terms.
#' 
#' @param x An R object of class `cSEMSurface`.
#' @param ... Currently ignored.
#' @export
plot.cSEMFloodlight <- function(x, ...) {
  
  ## Install ggplot2 if not already installed
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop2(
      "Package `plotly` required. Use `install.packages(\"plotly\")` and rerun.")
  }

  plot1 <- plotly::plot_ly( x = out$x, y = out$z, z = out$y, type = "surface")  
  options
  plot1
}