#' `cSEMPredict` method for `plot()`
#'
#' The `cSEMPredict` method for the generic function [plot()].
#'
#' @param x An R object of class `cSEMPredict`.
#' @param ... Currently ignored.
#'
#' @seealso [predict()]
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
    ggplot2::ggplot(xx, ggplot2::aes(x = xx[, "Actual"], y = xx[, "Predictions_target"])) +
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