#' `cSEMNonlinearEffects` method for `plot()`
#'
#' This plot method can be used to create plots to analyze non-linear models in more
#' depth. In doing so the following plot types can be selected:
#' \describe{
#' \item{`.plot_type = "simpleeffects"`:}{
#' The plot of a simple effects analysis displays the predicted value of the 
#' dependent variable for different values of the independent variable and the moderator.
#' As levels for the moderator the levels provided to the `doNonlinearEffectsAnalysis()` function
#' are used. Since the constructs are standardized the values of the moderator
#' equals the deviation from its mean measured in standard deviations.  
#' }
#' \item{`.plot_type = "surface"`:}{
#' The plot of a surface analyis displays the predicted values of an 
#' independent variable (z). The values are predicted based on the values of the moderator 
#' and the independent variable including all their higher-order terms. 
#' For the values of the moderator and the indepedent variable steps between
#'  their minimum and maximum values  are used.
#' }
#' \item{`.plot_type = "floodlight"`:}{
#' The plot of a floodlight analysis displays the direct effect of an continuous 
#' independent variable (z) on a dependent variable (y) conditional on the values
#' of a continuous moderator variable (x), including
#' the confidence interval and the Johnson-Neyman points. It is 
#' noted that in the floodlight plot only moderation is taken into account and higher
#' order terms are ignored. For more details, see \insertCite{Spiller2013;textual}{cSEM}. 
#' }
#' }
#' Plot the predicted values of an independent variable (z) 
#' The values are predicted based on a certain moderator and a certain 
#' independent variable including all their higher-order terms.
#' 
#' @param x An R object of class `cSEMNonlinearEffects`.
#' @param .plot_type A character string indicating the type of plot that should be produced.
#' Options are "*simpleeffects*", "*surface*", and "*floodlight*". Defaults to "*simpleeffects*".
#' @param .plot_package A character vector indicating the plot package used. Options are
#' "*plotly*", and "*persp*". Defaults to "*plotly*". 
#' @param ... Additional parameters that can be passed to 
#' \code{\link[graphics:persp]{graphics::persp}}, e.g., to rotate the plot.
#' 
#' @seealso [doNonlinearEffectsAnalysis()]
#' @export

plot.cSEMNonlinearEffects <- function(
  x,
  .plot_type = 'simpleeffects',
  .plot_package = 'plotly',
  ...) {
  
  if(!inherits(x, "cSEMNonlinearEffects")) {
    stop2("x must be of class `cSEMNonlinearEffects`")
  }
  
  if(!(.plot_type %in% c('simpleeffects','surface','floodlight'))){
    stop2("Currenlty only the plot types `simpleeffects`, `surface`, and `floodlight` are supported.")
  }
  
  ## Create Simple Effects plot -----------------------
  if(.plot_type == 'simpleeffects'){
    ## Install ggplot2 if not already installed
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop2(
        "Package `ggplot2` required. Use `install.packages(\"ggplot2\")` and rerun.")
    }
    
    # places_mod_ordered <- levels(with(x$out, reorder(values_mod,values_ind, values_dep)))
    plot1 <- ggplot2::ggplot(x$out_simpleeffects, 
                             ggplot2::aes(x = x$out_simpleeffects[,'values_ind'], 
                                          y = x$out_simpleeffects[,'values_dep'],
                                          group = x$out_simpleeffects[,'values_mod']))+
      ggplot2::geom_line(ggplot2::aes(linetype=x$out_simpleeffects[,'values_mod']))  +
      ggplot2::labs(
        x = paste('Level of ', x$Information['independent']),
        y = paste('Expected value of', x$Information['dependent']),
        # linetype=paste("Level of",x$Information['moderator']),
        caption = paste0("Created using cSEM version: ", packageVersion("cSEM"),
                         " (", Sys.Date(), ")")) +
      ggplot2::scale_linetype_discrete(name = paste("Level of",x$Information['moderator'],'\n',
                                                    'SDs from the mean'),
                                       breaks = sort(as.numeric(levels(x$out_simpleeffects$values_mod))))+
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    
    # Plot
    return(plot1)
  }# end of simple effects plot
  
  ## Create Surface plot -----------------------
  if(.plot_type == 'surface'){
  if(!(.plot_package %in% c('plotly','rsm','persp'))){
    stop2("Currenlty only plotly and persp are supported for plotting.")
  }
  
  if(.plot_package == 'plotly'){
    
    ## Install plotly if not already installed
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop2(
        "Package `plotly` required. Use `install.packages(\"plotly\")` and rerun.")
    }
  
    
    
    plot1 <- plotly::plot_ly( x = x$out_surface$values_ind1, 
                              y = x$out_surface$values_ind2, 
                              z = x$out_surface$values_dep, 
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
          xaxis = list(title = paste('Standardized ',x$Information$independent)),
          yaxis = list(title = paste('Standardized ',x$Information$moderator)),
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
  
  if(.plot_package == 'persp'){
    # using the persp function
    # phi and theta are used for rotation
    graphics::persp(x = x$out_surface$values_ind1, y = x$out_surface$values_ind2, z = t(x$out_surface$values_dep),
                    xlab = paste('Standardized ',x$Information$independent),
                    ylab = paste('Standardized ',x$Information$moderator),
                    zlab = x$Information$dependent,
                    main = paste('Surface analysis: ',x$Information$dependent ),...)
    
  }
  # 
  # if(.plot_package == 'rsm'){
  # # see rsm package for other plots
  # 
  # }
  } #end surface plot
  
  ## Create floodlight plot -----------------------
  if(.plot_type == 'floodlight'){
    ## Install ggplot2 if not already installed
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop2(
        "Package `ggplot2` required. Use `install.packages(\"ggplot2\")` and rerun.")
    }
    
    plot1 <- ggplot2::ggplot(as.data.frame(x$out_floodlight$out), 
                             ggplot2::aes(x = x$out_floodlight$out[, 'value_z'],
                                          y = x$out_floodlight$out[, 'direct_effect'])) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = x$out_floodlight$out[,'lb'], 
                                        ymax = x$out_floodlight$out[, 'ub'],
                                        fill=paste0(100*(1-x$Information$alpha),'% CI')), 
                           alpha = 0.2) +
      ggplot2::scale_fill_manual('',values="grey12")+
      ggplot2::labs(
        x = paste('Level of ', x$Information['moderator']), 
        y = paste('Effect of', x$Information['independent'], 'on \n', x$Information['dependent']),
        caption = paste0("Created using cSEM version: ", packageVersion("cSEM"), 
                         " (", Sys.Date(), ")")) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    
    # Add Johnson-Neyman points, if they exist in the considered range
    if(nrow(x$out_floodlight$Johnson_Neyman_points)>0){
    for(row in 1:nrow(x$out_floodlight$Johnson_Neyman_points)){
      plot1 <- plot1 +
            ggplot2::geom_point(x = x$out_floodlight$Johnson_Neyman_points[row,'x'], 
                                y = x$out_floodlight$Johnson_Neyman_points[row,'y'], size = 2)
    }
}
    # Plot
    return(plot1)
  }
}