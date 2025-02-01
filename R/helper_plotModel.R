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
get_significance_stars <- function(
    .pvalue
    ){
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
save_single_plot <- function(
    diagrammer_obj, 
    out_file){
  extension <- tolower(tools::file_ext(out_file))  # Get the file extension
  if (extension == "pdf") {
    export_svg(diagrammer_obj) %>% charToRaw() %>% rsvg::rsvg_pdf(out_file)
    message("Plot saved to ", out_file)
  } else if (extension == "png") {
    export_svg(diagrammer_obj) %>% charToRaw() %>% rsvg::rsvg_png(out_file)
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
#'
#' # Example: Save a plot object to a pdf file
#' plot_object <- plotModel(.object)
#' plotModelSave(plot_object, .file = "sem_plot.pdf")
#' }
#' @export
plotModelSave <- function(
    plot_object,
    .file) {
  if (missing(plot_object)) {
    stop2("plot_object must be provided to plotModelSave().")
  }
  
  # Ensure .file is provided
  if (missing(.file)) {
    stop2("Filename must be specified.")
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




#' Internal: secondOrderMeasurementEdges
#'
#' Build measurement edges for a first–order model.
#' 
#'
#' @usage firstOrderMeasurementEdges(
#'     construct, 
#'     weights,
#'     loadings, 
#'     weight_p_values, 
#'     loading_p_values, 
#'     plot_signif, 
#'     constructTypes
#'     )
#'
#' 
#' @inheritParams csem_arguments
#' 
#' @return Character string. .
#' 
#' @keywords internal
#'
firstOrderMeasurementEdges <- function(
    construct, 
    weights,
    loadings, 
    weight_p_values, 
    loading_p_values, 
    plot_signif, 
    constructTypes){
  
  code <- paste0("subgraph cluster_", construct, " {\nlabel=\"\";\nstyle=invis;\n")
  if (!is.null(weights[construct, ])) {
    non_zero_idx <- which(weights[construct, ] != 0)
    if (length(non_zero_idx) > 0) {
      indicators <- colnames(weights)[non_zero_idx]
      for (indicator in indicators) {
        code <- paste0(code, indicator, " [shape=box];\n")
        weight <- round(weights[construct, indicator], 3)
        loading <- round(loadings[construct, indicator], 3)
        weight_name <- paste(construct, "<~", indicator)
        loading_name <- paste(construct, "=~", indicator)
        weight_stars <- if (plot_signif) get_significance_stars(weight_p_values[weight_name]) else ""
        loading_stars <- if (plot_signif) get_significance_stars(loading_p_values[loading_name]) else ""
        label <- if (constructTypes[construct] == "Common factor") {
          paste0(loading, loading_stars)
        } else {
          paste0(weight, weight_stars)
        }
        direction <- if (constructTypes[construct] == "Common factor") {
          paste0(construct, " -> ", indicator)
        } else {
          paste0(indicator, " -> ", construct)
        }
        code <- paste0(code, direction, " [label=\"", label, "\"];\n")
      }
    }
  }
  code <- paste0(code, "}\n")
  return(code)
}


#' Internal: secondOrderMeasurementEdges
#'
#' Build measurement edges for a second–order model. 
#' 
#'
#'@usage secondOrderMeasurementEdges(
#'  construct,
#'  weights_first,
#'  loadings_first,
#'  weight_p_first,
#'  loading_p_first,
#'  weights_second,
#'  weight_p_second,
#'  plot_signif, 
#'  constructTypes,
#'  only_second_stage = FALSE
#'  )
#'
#' 
#' @inheritParams csem_arguments
#' 
#' @return Character string. .
#' 
#' @keywords internal
#'
secondOrderMeasurementEdges <- function(
    construct,
    weights_first,
    loadings_first,
    weight_p_first, 
    loading_p_first,
    weights_second, 
    weight_p_second,
    plot_signif,
    constructTypes,
    only_second_stage = FALSE
    ){
  
  code <- paste0("subgraph cluster_", gsub("[^A-Za-z0-9]", "_", construct),
                 " {\nlabel=\"\";\nstyle=invis;\n")
  if (!only_second_stage) {
    # Draw first–stage measurement edges if not skipping them.
    if (!is.null(weights_first) && (construct %in% rownames(weights_first))) {
      non_zero_idx <- which(weights_first[construct, ] != 0)
      if (length(non_zero_idx) > 0) {
        indicators <- colnames(weights_first)[non_zero_idx]
        for (indicator in indicators) {
          code <- paste0(code, "\"", indicator, "\"", " [shape=box];\n")
          weight <- round(weights_first[construct, indicator], 3)
          loading <- round(loadings_first[construct, indicator], 3)
          weight_name <- paste(construct, "<~", indicator)
          loading_name <- paste(construct, "=~", indicator)
          weight_stars <- if (plot_signif) get_significance_stars(weight_p_first[weight_name]) else ""
          loading_stars <- if (plot_signif) get_significance_stars(loading_p_first[loading_name]) else ""
          label <- if (constructTypes[construct] == "Common factor") {
            paste0(loading, loading_stars)
          } else {
            paste0(weight, weight_stars)
          }
          direction <- if (constructTypes[construct] == "Common factor") {
            paste0("\"", construct, "\" -> \"", indicator, "\"")
          } else {
            paste0("\"", indicator, "\" -> \"", construct, "\"")
          }
          code <- paste0(code, direction, " [label=\"", label, "\"];\n")
        }
      }
    }
  }
  # Always draw second–stage measurement edges as dashed indicator -> construct.
  if (!is.null(weights_second) && (construct %in% rownames(weights_second))) {
    non_zero_idx <- which(weights_second[construct, ] != 0)
    if (length(non_zero_idx) > 0) {
      indicators <- colnames(weights_second)[non_zero_idx]
      for (indicator in indicators) {
        code <- paste0(code, "\"", indicator, "\"", " [shape=box];\n")
        weight <- round(weights_second[construct, indicator], 3)
        weight_name <- paste(construct, "<~", indicator)
        weight_stars <- if (plot_signif) get_significance_stars(weight_p_second[weight_name]) else ""
        label <- paste0(weight, weight_stars)
        code <- paste0(code, "\"", indicator, "\"", " -> ", "\"", construct, "\"",
                       " [label=\"", label, "\", style=dashed];\n")
      }
    }
  }
  code <- paste0(code, "}\n")
  return(code)
}


#' Internal: buildDotCode
#'
#' Build DOT code from common components.
#' 
#'
#'@usage buildDotCode(
#'title, 
#'graph_attrs,
#'constructs, 
#'r2_values,
#'measurement_edge_fun,
#'path_coefficients, 
#'path_p_values,
#'correlations,
#'plot_significances,
#'plot_indicator_correlations,
#'plot_structural_model_only,
#'is_second_order = FALSE
#')
#'
#' 
#' @inheritParams csem_arguments
#' 
#' @return Character string. .
#' 
#' @keywords internal
#'
buildDotCode <- function(title, 
                         graph_attrs,
                         constructs, 
                         r2_values,
                         measurement_edge_fun,
                         path_coefficients, 
                         path_p_values,
                         correlations,
                         plot_significances,
                         plot_indicator_correlations,
                         plot_structural_model_only,
                         is_second_order = FALSE
                         ){
  
  dot_code <- "digraph SEM {\n"
  
  # Add title and basic graph attributes
  if (title != "") {
    dot_code <- paste0(dot_code,
                       "labelloc=\"t\";\n",
                       "label=<", title, " >;\n",
                       "fontsize=20;\n",
                       "fontname=\"Helvetica\";\n")
  }
  if (!is.null(graph_attrs)) {
    for (attr in graph_attrs) {
      dot_code <- paste0(dot_code, attr, ";\n")
    }
  }
  
  # Define nodes for constructs and add measurement subgraphs.
  for (construct in names(constructs)) {
    constrType <- constructs[[construct]]
    shape <- if (constrType == "Composite") "hexagon" else "circle"
    fixed_size <- "width=1.5, height=1.5, fixedsize=true"
    if (!is.na(r2_values[construct])) {
      r2_label <- paste0("R\u00b2 = ", round(r2_values[construct], 3))
      dot_code <- paste0(dot_code, construct,
                         " [label=\"", construct, "\\n", r2_label,
                         "\", shape=", shape, ", ", fixed_size, "];\n")
    } else {
      dot_code <- paste0(dot_code, construct,
                         " [label=\"", construct,
                         "\", shape=", shape, ", ", fixed_size, "];\n")
    }
    # For second-order models, always add measurement edges (which may be drawn only from second stage)
    if (is_second_order) {
      dot_code <- paste0(dot_code, measurement_edge_fun(construct))
    } else {
      if (!plot_structural_model_only) {
        dot_code <- paste0(dot_code, measurement_edge_fun(construct))
      }
    }
  }
  
  # Add structural model (path) edges.
  for (dependent in rownames(path_coefficients)) {
    depPaths <- path_coefficients[dependent, ]
    predictors <- names(depPaths[which(depPaths != 0)])
    for (predictor in predictors) {
      coefficient <- round(depPaths[predictor], 3)
      path_name <- paste(dependent, "~", predictor)
      stars <- if (plot_significances) get_significance_stars(path_p_values[path_name]) else ""
      if (grepl("\\.", predictor)) {
        # For predictors with a dot, create a diamond node.
        predictor_clean <- gsub("_temp", "", predictor)
        dot_code <- paste0(dot_code,
                           "\"", predictor_clean, "\"",
                           " [label=\"", predictor_clean, "\", shape=diamond, width=1.5, height=1.5, fixedsize=true];\n")
        dot_code <- paste0(dot_code,
                           "\"", predictor_clean, "\"", " -> ", dependent,
                           " [label=<", coefficient, stars, ">",
                           ", style=dashed, penwidth=2.5];\n")
      } else {
        dot_code <- paste0(dot_code,
                           predictor, " -> ", dependent,
                           " [label=<", coefficient, stars, ">",
                           ", penwidth=2.5];\n")
      }
    }
  }
  
  # Validate the .plot_indicator_correlations argument.
  valid_options <- c("exo", "both", "none")
  if (!plot_indicator_correlations %in% valid_options) {
    stop("Invalid argument for .plot_indicator_correlations. Valid options are 'exo', 'both', 'none'.")
  }
  
  # Add dashed edges for exogenous construct correlations.
  if (plot_indicator_correlations %in% c("exo", "both")) {
    for (i in seq_along(correlations$exo$names)) {
      corr <- round(correlations$exo$estimates[i], 3)
      stars <- if (plot_significances) get_significance_stars(correlations$exo$p_values[i]) else ""
      names_split <- strsplit(correlations$exo$names[i], " ~~ ")[[1]]
      dot_code <- paste0(dot_code,
                         names_split[1], " -> ", names_split[2],
                         " [label=<(", corr, stars, ")>",
                         ", style=dashed, dir=none];\n")
    }
  }
  
  # Add dashed edges for indicator correlations if requested.
  if (plot_indicator_correlations == "both") {
    for (i in seq_along(correlations$ind$names)) {
      corr <- round(correlations$ind$estimates[i], 3)
      stars <- if (plot_significances) get_significance_stars(correlations$ind$p_values[i]) else ""
      names_split <- strsplit(correlations$ind$names[i], " ~~ ")[[1]]
      dot_code <- paste0(dot_code,
                         names_split[1], " -> ", names_split[2],
                         " [label=<(", corr, stars, ")>",
                         ", style=dashed, dir=none];\n")
    }
  }
  
  dot_code <- paste0(dot_code, "}")
  return(dot_code)
}