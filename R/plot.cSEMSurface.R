#' `cSEMSurface` method for `plot()`
#'
#' Plot the predicted values of an independent variable (z) 
#' The values are predicted based on a certain moderator and a certain 
#' independent variable including all their higher-order terms.
#' 
#' @param x An R object of class `cSEMSurface`.
#' @param .plot_type A character vector indicating the plot package used. Options are
#' "*plotly*", and "*persp*". Defaults to "*plotly*". 
#' @param ... Additional parameters that can be passed to 
#' \code{\link[graphics:persp]{graphics::persp}}, e.g., to rotate the plot.
#' 
#' @seealso [doSurfaceAnalysis()]
#' @export

plot.cSEMSurface <- function(
  x,
  .plot_type = c('plotly'),
  ...) {
  
  if(!inherits(x, "cSEMSurface")) {
    stop2("x must be of class `cSEMSurface`")
  }
  
  if(!(.plot_type %in% c('plotly','rsm','persp'))){
    stop2("Currenlty only plotly and rsm are supported for plotting.")
  }
  
  if(.plot_type == 'plotly'){
    
  ## Install plotly if not already installed
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop2(
      "Package `plotly` required. Use `install.packages(\"plotly\")` and rerun.")
  }
    
  plot1 <- plotly::plot_ly( x = x$out$values_ind1,
                            y = x$out$values_ind2,
                            z = x$out$values_dep,
                            type = "surface",
                            colors = 'Greys') %>%

    plotly::add_surface(
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=TRUE,
          highlightcolor="#ff0000",
          project=list(z=TRUE)
        ))) %>%
    plotly::layout( # Axis labeling
      title = paste("Surface analysis:",x$Information$dependent),
      scene = list(
        xaxis = list(title = paste('Standardized ',x$Information$independent_1)),
        yaxis = list(title = paste('Standardized ',x$Information$independent_2)),
        zaxis = list(title = x$Information$dependent)
      ))%>% plotly::hide_colorbar()

  # save current viewer settings. Set to null that figure is opened in browser
   op=options()
  viewer_old <- op$viewer
  options(viewer = NULL)
  # Create plot
  return(plot1)
  # Restore viewer settings
  options(viewer = viewer_old)
  }

  if(.plot_type == 'persp'){
    # using the persp function
    # phi and theta are used for rotation
    graphics::persp(x = x$out$values_ind1, y = x$out$values_ind2, z = t(x$out$values_dep),
          xlab = paste('Standardized ',x$Information$independent_1),
          ylab = paste('Standardized ',x$Information$independent_2),
          zlab = x$Information$dependent,
          main = paste('Surface analysis: ',x$Information$dependent ),...)
    
  }
  # 
  # if(.plot_type == 'rsm'){
  # # see rsm package for other plots
  # 
  # }
}