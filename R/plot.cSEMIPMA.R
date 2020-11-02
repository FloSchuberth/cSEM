#' `cSEMIPMA` method for `plot()`
#'
#' Plot the importance-performance matrix.
#' 
#' @param x An R object of class `cSEMIPMA`.
#' @param .dependent Character string. Name of the target construct for which the 
#'   importance-performance matrix should be created.  
#' @param .attributes Character string. A vector containing indicator/construct 
#'   names that should be plotted in the importance-performance matrix. 
#'   It must be at least of length 2.
#' @param .level Character string. Indicates the level for which the 
#'   importance-performance matrix should be plotted. One of `"construct"` or 
#'   `"indicator"`. Defaults to `"construct"`.
#' @param ... Currently ignored.
#' 
#' @seealso [doIPMA()]
#' @export

plot.cSEMIPMA <- function(
  x           = NULL,
  .dependent  = NULL,
  .attributes = NULL,
  .level      = c("construct", "indicator"),
  ...
  ){
  
  .level <- match.arg(.level)
  

  # Check wether the input is of class cSEMIPMA
  if(!inherits(x, "cSEMIPMA")) {
    stop2("x must be of class `cSEMIPMA`.")
  }

  
  # check whether .dependent is supplied:
  if(is.null(.dependent)){
    stop2("Please provide the name of the dependent construct.")
  }
  
  # check whether .dependent is supplied:
  if(is.null(.attributes)){
    stop2("Please provide the names of the constructs/indicators that should be plotted in the IPM.")
  }
    
  # Check whether specified dependent variable is valid
  if(!.dependent %in% x$Construct_names | !is.character(.dependent)){
    stop2(".dependent is not a valid construct name.")
  }
  
  # Check whether specified attributes are valid
  if(.level == "construct"&!all(.attributes %in% x$Construct_names)){
    stop2("If `.level == 'construct'`, .attributes must only contain valid construct names.")
  }
  
  if(.level=="indicator"&!all(.attributes %in% x$Indicator_names)){
    stop2("If `.level == 'indicator'`, .attributes must only contain valid indicator names.")
  }
  
  # If only one attribute is provided the axis cannot be properly be scaled 
  # as I currently use the sd for this purpose. IPMA should be only conducted 
  # if more than one attribute is considered. 
  if(length(.attributes)<2){
    stop2("You need to provide at least construct/indicators names", 
          " to argument .attributes.")
  }
  
  ## Install ggplot2 if not already installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop2(
      "Package `ggplot2` required. Use `install.packages(\"ggplot2\")` and rerun.")
  }
  
  # Distinguihs between IPMA for constructs and for indicators
  if(.level == 'construct'){
    Importance <- x$Construct$Importance[.dependent,.attributes]
    # indep_vars <- colnames(Importance_all)[Importance_all!=0]
    # Importance <- Importance_all[,attributes]
    Performance <- x$Construct$Performance[.attributes]
    
  } else if(.level == 'indicator'){
    Importance <- x$Indicator$Importance[.dependent,.attributes] 
    # indep_vars <- colnames(Importance_all)[Importance_all!=0]
    # Importance <- Importance_all[,indep_vars]
    Performance <- x$Indicator$Performance[.attributes]
  }
  
  # prepare data for plotting
  data_plot <- data.frame("Importance"  = Importance,
                          "Performance" = Performance,
                          "Name"        = .attributes)
  
  # Calculate coordinates of the lines to seperate the quadrants
 # browser()
  horizontal <- mean(Performance) 
  vertical   <- mean(Importance)
  
  # Make quadrants symmetric
  # Determine the start and endpoint for the y-axis (Performance)
  devfrommeanPerf<-abs(c(min(Performance),max(Performance))-horizontal)
  if(devfrommeanPerf[1]>=devfrommeanPerf[2]){
    starty<-min(Performance)-sd(Performance)/2
    endy<-horizontal+devfrommeanPerf[1]+sd(Performance)/2
  }else{
    starty<-horizontal-devfrommeanPerf[2]-sd(Performance)/2
    endy<-max(Performance)+sd(Performance)/2
  }
  
  # Determine the start and endpoint for the x-axis (Importance)
  devfrommeanImp <-abs(c(min(Importance),max(Importance))-vertical)
  if(devfrommeanImp[1]>=devfrommeanImp[2]){
    startx<-min(Importance)-sd(Importance)/2
    endx<-vertical+devfrommeanImp[1]+sd(Importance)/2
  }else{
    startx<-vertical-devfrommeanImp[2]-sd(Importance)/2
    endx<-max(Importance)+sd(Importance)/2
  }
  
  plot1 <- ggplot2::ggplot(data_plot, ggplot2::aes(x = data_plot[, "Importance"], 
                                                   y = data_plot[, "Performance"],
                                                   label = data_plot[, "Name"])) + 
    ggplot2::geom_text(ggplot2::aes(label = data_plot[, "Name"]), hjust = 0, vjust = 1.2) +
    ggplot2::geom_point() + 
    ggplot2::geom_vline(xintercept = vertical) +
    ggplot2::geom_hline(yintercept = horizontal) +
    ggplot2::labs(
      x = "Importance", 
      y = "Performance",
      caption = paste0("Created using cSEM version: ", packageVersion("cSEM"), 
                       " (", Sys.Date(), ")")) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    # ggplot2::coord_fixed(ratio = 1, xlim = c(startx,endx), ylim = c(starty,endy), expand = TRUE, clip = "on")
    ggplot2::xlim(startx,endx) +
    ggplot2::ylim(starty,endy)
  
  # Return plot
  plot1
}
