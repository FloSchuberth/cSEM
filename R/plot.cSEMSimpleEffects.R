#' `cSEMSimpleEffects` method for `plot()`
#'
#' Creates a simple effects plot. The expected value of the dependent variable is plotted
#' for different values of the independent variable and the moderator. As levels for the moderator
#' -2, -1, 0, 1, and 2 are used. Since the constructs are standardized the values of the moderator
#' equals the deviation from its mean measured in standard deviations.   
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

  plot1 <- ggplot2::ggplot(x$out, ggplot2::aes(x = x$out[,'values_ind'], 
                                               y = x$out[,'values_dep'],
                                               group = x$out[,'values_mod']))+
    ggplot2::geom_line(ggplot2::aes(linetype=x$out[,'values_mod']))  +
    ggplot2::labs(
      x = paste('Level of ', x$Information['independent']),
      y = paste('Expected value of', x$Information['dependent']),
      # linetype=paste("Level of",x$Information['moderator']),
      caption = paste0("Created using cSEM version: ", packageVersion("cSEM"),
                       " (", Sys.Date(), ")")) +
    ggplot2::scale_linetype_discrete(name = paste("Level of",x$Information['moderator'],'\n',
                                                  'SDs from the mean'))+
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  
  # Plot
  plot1
}
