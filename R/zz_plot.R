#' `cSEMPredict` method for `plot()`
#'
#' The `cSEMPredict` method for the generic function [plot()].
#'
#' @inheritParams csem_arguments
#'
#' @seealso [csem()], [predict()]
#'
#' @export
#' @keywords internal
#' 
plot.cSEMPredict <- function(x, ...) {
  
  ## Install ggplot2 if not already installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop2(
      "Package `ggplot2` required. Use `install.packages(\"ggplot2\")` and rerun.")
  }
  
  ## Install tidyr and dplyr if not already installed
  if (!requireNamespace("dplyr", quietly = TRUE) | !requireNamespace("tidyr", quietly = TRUE)) {
    stop2(
      "Package `dplyr` and `tidyr` required. Use `install.packages(\"dplyr\")` and `install.packages(\"tidyr\")` and rerun.")
  }
  x <- lapply(x[1:4], as.data.frame)
  xx <- x %>% 
    lapply(function(x) {
      x$Obs <- as.numeric(rownames(x))
      x
    }) %>% 
    dplyr::bind_rows(.id = "Type") %>% 
    tidyr::pivot_longer(cols = colnames(x[[1]]), names_to = "Indicator") %>% 
    tidyr::pivot_wider(names_from = "Type", values_from = "value") %>% 
    as.data.frame()
  
  ## Actual vs. predicted
  plot_actual_vs_predicted <- 
    ggplot2::ggplot(xx, ggplot2::aes(x = xx[, "Actual_target"], y = xx[, "Predictions_target"])) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(. ~ Indicator, ncol = 4) + 
    ggplot2::geom_smooth(method = "lm", se = FALSE) +
    ggplot2::labs(
      title = "Actual vs. predicted values",
      x = "Actual values",
      y = "Predicted values",
      caption = "Created using cSEM"
    )
  
  ## Density predicted
  plot_density_predicted <- ggplot2::ggplot(as.data.frame(xx),  ggplot2::aes(x = xx[, "Predictions_target"])) +
    ggplot2::geom_density() +
    ggplot2::facet_wrap(. ~ Indicator, ncol = 4) + 
    ggplot2::labs(
      title = "Estimated density function - predicted values",
      x = "Predicted values",
      y = "Density",
      caption = "Created using cSEM"
    )
  
  list("Actual vs. predicted" = plot_actual_vs_predicted, 
       "Density predicted" = plot_density_predicted)
}

#' `cSEMFloodlight` method for `plot()`
#'
#' Plot the direct effect of an independent variable (z) on a dependent variable
#' (y) conditional on the values of a moderator variable (x), including
#' the confidence interval and the Johnson-Neyman points. 
#' 
#' @param x An R object of class `cSEMFloodlight` resulting from a call to [doFloodlightAnalysis].
#' @param ... ignored.
#' 
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
      "x" = paste('Level of ', x$Information['moderator']), 
      "y" = paste('Effect of', x$Information['independent'], 'on \n', x$Information['dependent'])) +
    ggplot2::theme_bw() +
    # scale_x_continuous(breaks=seq(-3,3,0.5))+
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