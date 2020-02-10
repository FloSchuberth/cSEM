#' `cSEMSurface` method for `plot()`
#'
#' Plot the predicted values of an independent variable (z) 
#' The values are predicted based on a certain moderator and a certain 
#' independent variable including all their higher-order terms.
#' 
#' @param x An R object of class `cSEMSurface`.
#' @param .plot_type A character vector indicating the plot package used. Options are
#' plotly, and rsm. Defaults to '*plotly*'. 
#' @param ... arguments that are passed to the used plotting function.
#' @export
plot.cSEMSurface <- function(x,
                             .plot_type = c('plotly'),
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

    # plot_ly(x = x$out$x, y = x$out$z, z = x$out$y, type = 'scatter3d', mode = 'lines')
    
  plot1 <- plotly::plot_ly( x = x$out$x, y = x$out$z, z = x$out$y, type = "surface",
                            , colors = 'Greys') %>%
   add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor="#ff0000",
          project=list(z=TRUE)
        )
      )) %>%
    layout( # Axis labeling
      title = paste("Surface Analysis:",x$Information$dependent),
      scene = list(
        xaxis = list(title = x$Information$independent),
        yaxis = list(title = x$Information$moderator),
        zaxis = list(title = x$Information$dependent)
      ))%>% hide_colorbar()
  op=options()
  viewer_old <- op$viewer
  options(viewer = NULL)
  return(plot1)
options(viewer = viewer_old)
  
  }

  if(.plot_type == 'rsm'){
  # see rsm package for other plots

  }
}