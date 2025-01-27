#' Internal: get significance stars
#'
#' Transforms a p-value into stars. 
#' 
#' `.pvalue` Numeric. A p-value that is transformed into a star.
#' 
#' 
#' @usage get_significance_stars(
#' .pvalue
#' )
#' 
#' @inheritParams csem_arguments
#' 
#' @return Character string. A p-value transformed into a star.
#' 
#' @keywords internal
#'
get_significance_stars <- function(.pvalue) {
  if (is.na(.pvalue)) return("")
  else if (.pvalue < 0.001) return("***")
  else if (.pvalue < 0.01) return("**")
  else if (.pvalue < 0.05) return("*")
  else return("")
}


#' Internal: save_single_plot
#' Helper function to save a single DiagrammeR plot based on the file extension
#' 
#' `diagrammer_obj` DiagrammeR plot object to be saved.
#' `out_file` The name of the file to save the plot to (supports 'pdf', 'png', 'svg', and 'dot' formats).
#' 
#' 
#' @usage save_single_plot(
#' diagrammer_obj, 
#' out_file
#' )
#' 
#' @inheritParams csem_arguments
#' 
#' @return NULL.
#' 
#' @keywords internal
#' 
save_single_plot <- function(diagrammer_obj, out_file) {
  extension <- tolower(tools::file_ext(out_file))  # Get the file extension
  if (extension == "pdf") {
    export_svg(diagrammer_obj) %>% charToRaw() %>% rsvg_pdf(out_file)
    message("Plot saved to ", out_file)
  } else if (extension == "png") {
    export_svg(diagrammer_obj) %>% charToRaw() %>% rsvg_png(out_file)
    message("Plot saved to ", out_file)
  } else if (extension == "svg") {
    svg_code <- export_svg(diagrammer_obj)
    writeLines(svg_code, con = out_file)
    message("Plot saved to ", out_file)
  } else if (extension == "dot") {
    dot_code <- diagrammer_obj$x$diagram  # Extract DOT code
    writeLines(dot_code, con = out_file)
    message("Plot saved to ", out_file)
  } else {
    stop("Unsupported file format! Please use 'pdf', 'png', 'svg', or 'dot'.")
  }
}

#' Save plotMOdel object to a file
#'
#' This function saves a given plot of a cSEM model to a specified file format.
#'
#'  @usage plotModelSave(
#'  plot_object,
#'  .file) 
#'
#' @param plot_object The DiagrammeR plot object to be saved. If not provided, the function will try to retrieve `latest_plot` from the global environment.
#' @param .file Character string. The name of the file to save the plot to (supports 'pdf', 'png', 'svg', and 'dot' formats).
#'
#' @return NULL (the function saves the plot to a file and provides messages on completion)
#'
#' @examples
#' \dontrun{
#' # Example 1: Automatically retrieve the latest plot and save as DOT format
#' plotModel(.object)
#' plotModelSave(.file = "sem_plot.dot")
#'
#' # Example 2: Save a plot object to a pdf file
#' plot_object <- plotModel(.object)
#' plotModelSave(plot_object, .file = "sem_plot.pdf")
#' }
#' @export
plotModelSave <- function(plot_object, .file) {
  
  # Check if plot_object is provided, otherwise try retrieving the latest plot from the global environment
  if (missing(plot_object)) {
    if (!exists("latest_plot", envir = .GlobalEnv)) {
      stop("No plot object provided and no plot available to save. Please run plotModel() first or provide a plot object.")
    }
    plot_object <- get("latest_plot", envir = .GlobalEnv) # Retrieve the latest plot if available
  }
  
  # Ensure .file is provided
  if (missing(.file)) {
    stop("Filename must be specified.")
  }
  
  # Check if the provided plot object is a multi-plot (list of plots)
  is_multi <- inherits(plot_object, "cSEMPlot_multi") 
  

  # Handle saving for multi-plot and single-plot cases
  if (is_multi) {
    # --- MULTI-PLOT CASE ---
    # Save each sub-plot in the list with an appended identifier in the filename
    for (i in seq_along(plot_object)) {
      this_plot  <- plot_object[[i]]
      this_name  <- names(plot_object)[i]  # Extract the name of the sub-plot
      
      # Generate a unique filename for each sub-plot
      base_name <- gsub("\\.(pdf|png|svg|dot)$", "", .file, ignore.case = TRUE)
      ext       <- tools::file_ext(.file)
      file_i    <- paste0(base_name, "_", this_name, ".", ext)
      
      # Save the sub-plot
      save_single_plot(this_plot, file_i)
    }
  } else {
    # --- SINGLE-PLOT CASE ---
    save_single_plot(plot_object, .file)
  }
  
  invisible(NULL)  # Return NULL invisibly
}