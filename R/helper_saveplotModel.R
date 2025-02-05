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
save_single_plot <- function(
    diagrammer_obj, 
    out_file
    ){
  extension <- tolower(tools::file_ext(out_file))  # Get the file extension
  if (extension == "pdf") {
    DiagrammeRsvg::export_svg(diagrammer_obj) %>% charToRaw() %>% rsvg::rsvg_pdf(out_file)
    message("Plot saved to ", out_file)
  } else if (extension == "png") {
    DiagrammeRsvg::export_svg(diagrammer_obj) %>% charToRaw() %>% rsvg::rsvg_png(out_file)
    message("Plot saved to ", out_file)
  } else if (extension == "svg") {
    svg_code <- DiagrammeRsvg::export_svg(diagrammer_obj)
    writeLines(svg_code, con = out_file)
    message("Plot saved to ", out_file)
  } else if (extension == "dot") {
    dot_code <- diagrammer_obj$x$diagram  # Extract DOT code
    writeLines(dot_code, con = out_file)
    message("Plot saved to ", out_file)
  } else {
    stop2("Unsupported file format! Please use 'pdf', 'png', 'svg', or 'dot'.")
  }
}
