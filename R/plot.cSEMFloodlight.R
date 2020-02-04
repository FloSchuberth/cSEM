#' `cSEMFloodlight` method for `plot()`
#'
#' Plot the direct effect of an independent variable (z) on a dependent variable
#' (y) conditional on the values of a moderator variable (x), including
#' the confidence interval and the Johnson-Neyman points. 
#' 
#' @param x An R object of class `cSEMFloodlight`.
#' @param ... Currently ignored.
#' @export
plot.cSEMFloodlight <- function(x, ...) {
  
  ## Install ggplot2 if not already installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop2(
      "Package `ggplot2` required. Use `install.packages(\"ggplot2\")` and rerun.")
  }
  
  plot1 <- ggplot2::ggplot(as.data.frame(x$out), ggplot2::aes(x = x$out[, 'value_z'], 
                                                              y = x$out[, 'direct_effect'])) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = x$out[,'lb'], ymax = x$out[, 'ub']), alpha = 0.2) +
    ggplot2::labs(
      x = paste('Level of ', x$Information['moderator']), 
      y = paste('Effect of', x$Information['independent'], 'on \n', x$Information['dependent']),
      caption = paste0("Created using cSEM version: ", packageVersion("cSEM"), 
                       " (", Sys.Date(), ")")) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  
  # Add Johnson-Neyman points, if they exist in the considered range
  if(length(x$Johnson_Neyman_points$JNlb) == 2){
    JN <- x$Johnson_Neyman_points$JNlb
    plot1 <- plot1 +
      ggplot2::geom_point(x = JN['x'], y = JN['y'], size = 2)  
  }
  
  if(length(x$Johnson_Neyman_points$JNub) == 2){
    JN = x$Johnson_Neyman_points$JNub
    plot1 = plot1 +
      ggplot2::geom_point(x = JN['x'], y = JN['y'], size = 2)  
  }
  
  # Plot
  plot1
}