#' `cSEMSimpleEffects` method for `plot()`
#'
#' Creates a simple effects plot. The expected value of the dependent variable is plotted
#' for different values of the independent variable and the moderator
#' 
#' @param x An R object of class `cSEMSimpleEffects`.
#' @param ... Currently ignored.
#' @export
plot.cSEMSimpleEffects <- function(x, ...) {
  
  if(!inherits(x, "cSEMSimpleEffects")) {
    stop2("x must be of class `cSEMSimpleEffects`.")
  }
  ## Install ggplot2 if not already installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop2(
      "Package `ggplot2` required. Use `install.packages(\"ggplot2\")` and rerun.")
  }
  browser()
  plot1 <- ggplot2::ggplot(x$out, ggplot2::aes(x = values_ind, 
                                               y = values_dep,
                                               color = values_mod))+
    ggplot2::geom_point() 
  # +
  #   ggplot2::geom_ribbon(ggplot2::aes(ymin = x$out[,'lb'], ymax = x$out[, 'ub']), alpha = 0.2) +
  #   ggplot2::labs(
  #     x = paste('Level of ', x$Information['moderator']), 
  #     y = paste('Effect of', x$Information['independent'], 'on \n', x$Information['dependent']),
  #     caption = paste0("Created using cSEM version: ", packageVersion("cSEM"), 
  #                      " (", Sys.Date(), ")")) +
  #   ggplot2::theme_bw() +
  #   ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  
  # Plot
  plot1
}