#' `cSEMIPMA` method for `plot()`
#'
#' Plot the importance-performance matrix
#' 
#' @param x An R object of class `cSEMIPMA`.
#' @param .dependent Character. Name of the target construct for which the 
#' importance performance matrix should be created.  
#' @param .attributes Character string. A vector containing indicator/construct names that 
#' should be plotted in the importance performance matrix. It must be at least of length 2.
#' @param .level Character. It indicates whether the importance performance matrix 
#' should be created on indicator or construct level. Options are `'construct'` and `'indicator'`.
#' Defaults to `'construct'`.
#' 
#' @seealso [doIPMA()]
#' @export

plot.cSEMIPMA <- function(x=NULL,
                          .dependent=NULL,
                          .attributes=NULL,
                          .level='construct'){
  
  # Check wether the input is of class cSEMIPMA
  if(!inherits(x, "cSEMIPMA")) {
    stop2("x must be of class `cSEMIPMA`.")
  }
  
  # Check whether specified dependent variable is valid
  if(!.dependent %in% x$Construct_names | !is.character(.dependent)){
    stop2(".dependent is not a valid construct name.")
  }
  
  # Check whether specified attributes are valid
  if(.level=="construct"&!all(.attributes %in% x$Construct_names)){
    stop2("If .level == `construct`, .attributes must only contain valid construct names.")
  }
  
  if(.level=="indicator"&!all(.attributes %in% x$Indicator_names)){
    stop2("If .level == `indicator`, .attributes must only contain valid indicator names.")
  }
  
  # If only one attribute is provided the axis cannot be properly scaled as I currently use the sd 
  # for this purpose. IPMA should be only conducted if more than one attribute is considered. 
  if(length(.attributes)<=2){
    stop2("You need to provide at least construct/indicators names to argument .attributes.")
  }
  
  ## Install ggplot2 if not already installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop2(
      "Package `ggplot2` required. Use `install.packages(\"ggplot2\")` and rerun.")
  }
  
  # Distinguihs between IPMA for constructs and for indicators
  if(.level=='construct'){
    Importance <- x$Construct$Importance[.dependent,.attributes]
    # indep_vars <- colnames(Importance_all)[Importance_all!=0]
    # Importance <- Importance_all[,attributes]
    Performance <- x$Construct$Performance[.attributes]
    
  }else if(.level=='indicator'){
    Importance <- x$Indicator$Importance[.dependent,.attributes] 
    # indep_vars <- colnames(Importance_all)[Importance_all!=0]
    # Importance <- Importance_all[,indep_vars]
    Performance <- x$Indicator$Performance[.attributes]
  }else{
    stop2("Allowed choices for .level are only 'construct' and 'indicator'.")
  }
  
  # prepare data for plotting
  data_plot <- data.frame("Importance"=Importance,
                          "Performance"=Performance,
                          "Name"= .attributes)
  
  # Calculate coordinates of the lines to seperate the quadrants
  horizontal <- mean(Performance) 
  vertical <- mean(Importance)
  
  plot1 <- ggplot2::ggplot(data_plot, ggplot2::aes(x = Importance, 
                                                   y = Performance,
                                                   label=Name))+
    ggplot2::geom_text(ggplot2::aes(label=Name),hjust=0, vjust=1.1)+
    ggplot2::geom_point()+
     ggplot2::geom_vline(xintercept = vertical)+
     ggplot2::geom_hline(yintercept = horizontal)+
    ggplot2::labs(
      x = "Importance", 
      y = "Performance",
      caption = paste0("Created using cSEM version: ", packageVersion("cSEM"), 
                       " (", Sys.Date(), ")"))+
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())+
    ggplot2::xlim(min(Importance)-sd(Importance)/2,max(Importance)+sd(Importance)/2)+
  ggplot2::ylim(min(Performance)-sd(Performance)/2,max(Performance)+sd(Performance)/2)
  
  # Return plot
  plot1
  
}
