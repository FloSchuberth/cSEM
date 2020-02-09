#' `cSEMSurface` method for `plot()`
#'
#' Plot the predicted values of an independent variable (z) 
#' The values are predicted based on a certain moderator and a certain 
#' independent variable including all their higher-order terms.
#' 
#' @param x An R object of class `cSEMSurface`.
#' @param ... Currently ignored.
#' @export
# plot.cSEMSurface <- function(x, ...) {
# 
#   ## Install ggplot2 if not already installed
#   if (!requireNamespace("plotly", quietly = TRUE)) {
#     stop2(
#       "Package `plotly` required. Use `install.packages(\"plotly\")` and rerun.")
#   }
# 
#   plot1 <- plotly::plot_ly( x = out$x, y = out$z, z = out$y, type = "surface")
#   options
#   plot1
#   
#   
#   # see rsm package for other plots
#   
#   # Needs to be added to dependencies
#   library(plotly)
#   # For plotting 3D figures options(viewer=NULL) must be set, perhaps there is a more elegant way
#   
#   plot_ly( x = ret$x, y = ret$z, z = ret$y, type = "scatter3d",mode='lines')
#   
#   scatter3d, mode='lines'.
#   
#   p1 <- plot_ly(x= ret$x, y=ret$y,z = ret$z, scene='scene1', lighting = list(ambient = 0.2)) %>%
#     +   add_surface(showscale=FALSE)
#   
#   p1 = plot_ly(ret, x = ~x, y = ~y, z = ~z) 
#   %>%
#     +   add_surface(p1,showscale=FALSE)
#   
#   
# }