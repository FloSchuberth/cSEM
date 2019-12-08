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
  xx <- lapply(x[1:4], as.data.frame)
  xx1 <- xx %>% 
    lapply(function(x) {
      x$Obs <- as.numeric(rownames(x))
      x
    }) %>% 
    dplyr::bind_rows(.id = "Type") %>% 
    tidyr::pivot_longer(cols = colnames(xx[[1]]), names_to = "Indicator") %>% 
    tidyr::pivot_wider(names_from = "Type", values_from = "value")
  
  ## Actual vs. predicted
  plot_actual_vs_predicted <- ggplot2::ggplot(xx1, ggplot2::aes(x = Actual_target, y = Predictions_target)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(. ~ Indicator, ncol = 4) + 
    ggplot2::geom_smooth(method = "lm", se = FALSE) +
    ggplot2::labs(
      title = "Actual vs. predicted values",
      x = "Actual values",
      y = "Predicted values",
      caption = "Created using cSEM"
    )
  
  plot_actual_vs_predicted
}