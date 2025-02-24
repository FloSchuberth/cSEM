#' savePlot
#'
#' This function saves a given plot of a cSEMResults object to a specified file format.
#'
#' @usage savePlot(
#'  .plot_object,
#'  .filename,
#'  .path = NULL) 
#'
#' @param .plot_object Object returned by one of the following functions [plot.cSEMResults_default()], [plot.cSEMResults_multi()], or
#' [plot.cSEMResults_2ndorder()].
#' @param .filename Character string. The name of the file to save the plot to (supports 'pdf', 'png', 'svg', and 'dot' formats).
#' @param .path Character string. Path of the directory to save the file to. Defaults
#'   to the current working directory.
#'
#' @return NULL
#'
#' @seealso [plot.cSEMResults_default()] [plot.cSEMResults_multi()] [plot.cSEMResults_2ndorder()]
#' 
#' @export
savePlot <- function(
    .plot_object,
    .filename,
    .path = NULL
    ){   
  
  if (missing(.plot_object)) {
    stop2(".plot_object must be provided to savePlot().")
  }
  
  # Ensure .filename is provided
  if (missing(.filename)) {
    stop2("Filename must be specified.")
  }
  
  ## Install DiagrammeR if not already installed
  if (!requireNamespace("DiagrammeRsvg", quietly = TRUE)) {
    stop2(
      "Package `DiagrammeRsvg` required. Use `install.packages(\"DiagrammeRsvg\")` and rerun.")
  }
  
  # Check if the provided plot object is a multi-plot (list of plots)
  is_multi <- inherits(.plot_object, "cSEMPlot_multi") 
  
  
  # Handle saving for multi-plot and single-plot cases
  if (is_multi) {
    # --- MULTI-PLOT CASE ---
    # Save each sub-plot in the list with an appended identifier in the filename
    for (i in seq_along(.plot_object)) {
      this_plot  <- .plot_object[[i]]
      this_name  <- names(.plot_object)[i]  # Extract the name of the sub-plot
      
      # Generate a unique filename for each sub-plot
      base_name <- gsub("\\.(pdf|png|svg|dot)$", "", .filename, ignore.case = TRUE)
      ext       <- tools::file_ext(.filename)
      file_i    <- paste0(base_name, "_", this_name, ".", ext)
      
      # Save the sub-plot
      save_single_plot(this_plot, file_i, .path) #added .path here
    }
  } else {
    # --- SINGLE-PLOT CASE ---
    save_single_plot(.plot_object, .filename, .path) #added .path here
  }
  
  invisible(NULL)  # Return NULL invisibly
}