#' `cSEMSurface` method for `plot()`
#'
#' Plot the predicted values of an independent variable (z) 
#' The values are predicted based on a certain moderator and a certain 
#' independent variable including all their higher-order terms.
#' 
#' @param x An R object of class `cSEMSurface`.
#' @param .plot_type A character vector indicating the plot package used. Options are
#' plotly, persp, and rsm. Defaults to '*plotly*'. 
#' @param ... arguments that are passed to the used plotting function.
#' @export
plot.cSEMSurface <- function(x,
                             .plot_type = c('plotly'),
                             ...) {
  if(!(.plot_type %in% c('plotly','rsm','persp'))){
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
                            colors = 'Greys') %>%
   # plotly::add_surface(
   #    contours = list(
   #      z = list(
   #        show=TRUE,
   #        usecolormap=TRUE,
   #        highlightcolor="#ff0000",
   #        project=list(z=TRUE)
   #      ))) %>%
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
        xaxis = list(title = paste('Standardized ',x$Information$independent)),
        yaxis = list(title = paste('Standardized ',x$Information$moderator)),
        zaxis = list(title = x$Information$dependent)
      ))%>% plotly::hide_colorbar()

  # save viewer settings. Set to null that figure is opened in browser
   op=options()
  viewer_old <- op$viewer
  options(viewer = NULL)
  return(plot1)
options(viewer = viewer_old)
  }

  if(.plot_type == 'persp'){
    # # Using the rgl package
    # ylim <- range(x$out$z)
    # ylen <- ylim[2] - ylim[1] + 1
    # 
    # colorlut <- terrain.colors(ylen)
    # col <- colorlut[ x$out$z - ylim[1] + 1 ]
    # rgl.surface(x = x$out$x, z = x$out$z, y = x$out$y,color = col, back="lines")
    # 
    # using the persp function
    # phi and theta are used for rotation
    persp(x = x$out$x, y = x$out$z, z = x$out$y, phi = 40, theta = 140,
          xlab = paste('Standardized ',x$Information$independent),
          ylab = paste('Standardized ',x$Information$moderator),
          zlab = x$Information$dependent,
          main = paste('Surface analysis: ',x$Information$dependent ))
    
  }
  
  
  if(.plot_type == 'rsm'){
  # see rsm package for other plots

  }
}