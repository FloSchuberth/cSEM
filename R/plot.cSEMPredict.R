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
  
  if(!inherits(x, "cSEMPredict")) {
    stop2("x must be of class `cSEMPredict`.")
  }
  
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
      caption = paste0("Created using cSEM version: ", packageVersion("cSEM"), 
                       " (", Sys.Date(), ")")
    )
  
  ## Density residuals
  plot_density_residuals <- ggplot2::ggplot(as.data.frame(xx),  ggplot2::aes(x = xx[, "Residuals_target"])) +
    ggplot2::geom_density() +
    ggplot2::facet_wrap(. ~ Indicator, ncol = 4) + 
    ggplot2::labs(
      title = "Estimated density function - residuals",
      x = "Residuals",
      y = "Density",
      caption = paste0("Created using cSEM version: ", packageVersion("cSEM"), 
                       " (", Sys.Date(), ")")
    )
  
  ## Density predicted
  plot_density_predicted <- ggplot2::ggplot(as.data.frame(xx),  ggplot2::aes(x = xx[, "Predictions_target"])) +
    ggplot2::geom_density() +
    ggplot2::facet_wrap(. ~ Indicator, ncol = 4) + 
    ggplot2::labs(
      title = "Estimated density function - predicted values",
      x = "Predicted values",
      y = "Density",
      caption = paste0("Created using cSEM version: ", packageVersion("cSEM"), 
                       " (", Sys.Date(), ")")
    )
  
  out <- list(
    "Actual vs. predicted" = plot_actual_vs_predicted, 
    "Density residuals" = plot_density_residuals,
    "Density predicted" = plot_density_predicted
    )
  class(out) <- "cSEMPlotPredict"
  out
}

#' `cSEMPlotPredict` method for `print()`
#'
#' The `cSEMPlotPredict` method for the generic function [print()].
#'
#' @param x An R object of class `cSEMPlotPredict`.
#' @param ... Currently ignored.
#'
#' @seealso [predict()], [plot.cSEMPredict()]
#'
#' @export
#' @keywords internal

print.cSEMPlotPredict <- function(x, ...) {
  cat2(rule2(type = 2))
  cat2("\n\tEnter the number of the plot you want to show:\n\n",
       "\t\t1 - Actual values vs. predicted values\n",
       "\t\t2 - Estimated density function - residuals\n",
       "\t\t3 - Estimated density function - predicted values"
       # \n\n",
       # "\tTo get the code that created the plots enter:\n\n",
       # "\t\t`code_x`\n\n", 
       # "\twhere `x` is the number of the plot you want the code for."
       )
  cat2("\n", rule2(type = 2))
  choice <- readline("Your choice:")
  if(!(choice %in% c(1:3, paste0("code_", 1:3)))) {
    stop2("`", choice, "`", " is not a valid choice.")
  } else if(choice == 1) {
    print(x$`Actual vs. predicted`)
  } else if(choice == 2) {
    print(x$`Density residuals`)
  } else if(choice == 3) {
    print(x$`Density predicted`)
  } 
  # else if(choice == "code_1") {
  #   cat2("Here is the code:")
  # }
}
