#' Internal: get significance stars
#'
#' Transforms a p-value into stars. 
#' 
#' `.pvalue` Numeric. A p-value that is transformed into significance stars.
#' 
#' @usage get_significance_stars(
#' .pvalue
#' )
#' @inheritParams csem_arguments
#' @return Character string. A p-value transformed into a star.
#' @keywords internal
get_significance_stars <- function(
    .pvalue
){
  if (is.na(.pvalue)) return("")
  else if (.pvalue < 0.001) return("***")
  else if (.pvalue < 0.01) return("**")
  else if (.pvalue < 0.05) return("*")
  else return("")
}

#' Internal: firstOrderMeasurementEdges
#'
#' Build measurement edges for a first–order model.
#'
#' @usage firstOrderMeasurementEdges(
#'     construct,
#'     weights,
#'     loadings,
#'     weight_p_values,
#'     loading_p_values,
#'     plot_signif,
#'     plot_labels,
#'     constructTypes
#'     )
#'
#' @inheritParams csem_arguments
#' @return Character string containing DOT code.
#' @keywords internal
firstOrderMeasurementEdges <- function(
    construct, 
    weights,
    loadings, 
    weight_p_values,
    loading_p_values, 
    plot_signif,
    plot_labels,
    constructTypes
){
  lines <- c()
  lines <- c(lines, paste0("subgraph cluster_", construct, " {"))
  lines <- c(lines, "label=\"\";")
  lines <- c(lines, "style=invis;")
  if (!is.null(weights[construct, ])) {
    non_zero_idx <- which(weights[construct, ] != 0)
    if (length(non_zero_idx) > 0) {
      indicators <- colnames(weights)[non_zero_idx]
      for (indicator in indicators) {
        lines <- c(lines, paste0(indicator, " [shape=box];"))
        weight <- round(weights[construct, indicator], 3)
        loading <- round(loadings[construct, indicator], 3)
        weight_name <- paste(construct, "<~", indicator)
        loading_name <- paste(construct, "=~", indicator)
        weight_stars <- if (plot_signif) get_significance_stars(weight_p_values[weight_name]) else ""
        loading_stars <- if (plot_signif) get_significance_stars(loading_p_values[loading_name]) else ""
        label <- if(plot_labels) {
          if (constructTypes[construct] == "Common factor") {
            paste0(loading, loading_stars)
          } else {
            paste0(weight, weight_stars)
          }
        } else {
          ""
        }
        lines <- c(lines, paste0(
          if(constructTypes[construct] == "Common factor") { paste0(construct, " -> ", indicator) } else { paste0(indicator, " -> ", construct) },
          " [label=\"", label, "\"];"
        ))
      }
    }
  }
  lines <- c(lines, "}")
  return(paste(lines, collapse = "\n"))
}

#' Internal: secondOrderMeasurementEdges
#'
#' Build measurement edges for a second–order model. 
#'
#'@usage secondOrderMeasurementEdges(
#'  construct,
#'  weights_first,
#'  loadings_first,
#'  weight_p_first,
#'  loading_p_first,
#'  weights_second,
#'  loadings_second,
#'  weight_p_second,
#'  loading_p_second,
#'  plot_signif,
#'  plot_labels,
#'  constructTypes,
#'  only_second_stage = FALSE
#'  )
#'
#' @inheritParams csem_arguments
#' @return Character string.
#' @keywords internal
secondOrderMeasurementEdges <- function(
    construct,
    weights_first,
    loadings_first,
    weight_p_first, 
    loading_p_first,
    weights_second,
    loadings_second,
    weight_p_second,
    loading_p_second,
    plot_signif,
    plot_labels,    
    constructTypes,
    only_second_stage = FALSE
){
  
  lines <- c()
  clust <- gsub("[^A-Za-z0-9]", "_", construct)
  lines <- c(lines, paste0("subgraph cluster_", clust, " {"))
  lines <- c(lines, "label=\"\";")
  lines <- c(lines, "style=invis;")
  
  ## First–stage measurement edges (solid)
  if (!only_second_stage) {
    if (!is.null(weights_first) && (construct %in% rownames(weights_first))) {
      non_zero_idx <- which(weights_first[construct, ] != 0)
      if (length(non_zero_idx) > 0) {
        indicators <- colnames(weights_first)[non_zero_idx]
        for (indicator in indicators) {
          if (cleanNode(construct) == cleanNode(indicator)) next
          lines <- c(lines, paste0("\"", indicator, "\" [shape=box];"))
          weight  <- round(weights_first[construct, indicator], 3)
          loading <- round(loadings_first[construct, indicator], 3)
          weight_name  <- paste(construct, "<~", indicator)
          loading_name <- paste(construct, "=~", indicator)
          weight_stars  <- if (plot_signif) get_significance_stars(weight_p_first[weight_name]) else ""
          loading_stars <- if (plot_signif) get_significance_stars(loading_p_first[loading_name]) else ""
          label <- if(plot_labels) {
            if (constructTypes[construct] == "Common factor") {
              paste0(loading, loading_stars)
            } else {
              paste0(weight, weight_stars)
            }
          } else {
            ""
          }
          lines <- c(lines, paste0(
            if (constructTypes[construct] == "Common factor") { paste0("\"", construct, "\" -> \"", indicator, "\"") } else { paste0("\"", indicator, "\" -> \"", construct, "\"") },
            " [label=\"", label, "\"];"
          ))
        }
      }
    }
  }
  
  ## Second–stage measurement edges (dashed)
  if (!only_second_stage) {
    if (constructTypes[construct] == "Common factor") {
      if (!is.null(loadings_second) && (construct %in% rownames(loadings_second)) &&
          any(loadings_second[construct, ] != 0)) {
        indicators <- colnames(loadings_second)[which(loadings_second[construct, ] != 0)]
        source_matrix <- "second"
      } else {
        indicators <- colnames(weights_first)[which(weights_first[construct, ] != 0)]
        source_matrix <- "first"
      }
      if (length(indicators) > 0) {
        for (indicator in indicators) {
          if (cleanNode(construct) == cleanNode(indicator)) next
          lines <- c(lines, paste0("\"", indicator, "\" [shape=circle, width=1.5, height=1.5, fixedsize=false];"))
          if (source_matrix == "second" && !is.na(loadings_second[construct, indicator]) && 
              loadings_second[construct, indicator] != 0) {
            loading <- round(loadings_second[construct, indicator], 3)
            loading_name <- paste(construct, "=~", indicator)
            loading_stars <- if (plot_signif) get_significance_stars(loading_p_second[loading_name]) else ""
          } else {
            loading <- round(loadings_first[construct, indicator], 3)
            loading_name <- paste(construct, "=~", indicator)
            loading_stars <- if (plot_signif) get_significance_stars(loading_p_first[loading_name]) else ""
          }
          label <- if(plot_labels) {
            paste0(loading, loading_stars)
          } else {
            ""
          }
          lines <- c(lines, paste0("\"", construct, "\" -> \"", indicator, "\" [label=\"", label, "\", penwidth=2.5, style=dashed];"))
        }
      }
    } else {
      if (!is.null(weights_second) && (construct %in% rownames(weights_second)) &&
          any(weights_second[construct, ] != 0)) {
        indicators <- colnames(weights_second)[which(weights_second[construct, ] != 0)]
        source_matrix <- "second"
      } else {
        indicators <- colnames(weights_first)[which(weights_first[construct, ] != 0)]
        source_matrix <- "first"
      }
      if (length(indicators) > 0) {
        for (indicator in indicators) {
          if (cleanNode(construct) == cleanNode(indicator)) next
          lines <- c(lines, paste0("\"", indicator, "\" [shape=circle, width=1.5, height=1.5, fixedsize=false];"))
          if (source_matrix == "second" && !is.na(weights_second[construct, indicator]) &&
              weights_second[construct, indicator] != 0) {
            weight <- round(weights_second[construct, indicator], 3)
            weight_name <- paste(construct, "<~", indicator)
            weight_stars <- if (plot_signif) get_significance_stars(weight_p_second[weight_name]) else ""
          } else {
            weight <- round(weights_first[construct, indicator], 3)
            weight_name <- paste(construct, "<~", indicator)
            weight_stars <- if (plot_signif) get_significance_stars(weight_p_first[weight_name]) else ""
          }
          label <- if(plot_labels) {
            paste0(weight, weight_stars)
          } else {
            ""
          }
          lines <- c(lines, paste0("\"", indicator, "\" -> \"", construct, "\" [label=\"", label, "\", penwidth=2.5, style=dashed];"))
        }
      }
    }
  }
  
  lines <- c(lines, "}")
  return(paste(lines, collapse = "\n"))
}

#' Internal: Clean a node name.
#'
#' Removes a trailing "_temp" from a node name.
#'
#' @param node A node name.
#'
#' @return A cleaned node name.
#' @keywords internal
cleanNode <- function(node) {
  gsub("_temp$", "", node)
}

#' Internal: Check whether two indicators belong to the same construct.
#'
#' Checks whether two indicators belong to the same construct.
#'
#' @param .indicator1 Character string. The name of the indicator 1.
#' @param .indicator2 Character string. The name of the indicator 1.
#' @param .model_measurement Matrix. The measurement matrix indicating the relationship
#' between constructs and indicators. 
#' @param .model_error_cor Matrix. The matrix indicates the error correlation structure.
#'
#' @return TRUE if both indicators belong to the same construct, FALSE otherwise.
#' @keywords internal
check_connection <- function(.indicator1, .indicator2, .model_measurement, .model_error_cor) {
  x <- .model_measurement[, .indicator1, drop = FALSE]
  construct1 <- rownames(x)[x == 1]
  x <- .model_measurement[, .indicator2, drop = FALSE]
  construct2 <- rownames(x)[x == 1]
  if (construct1 == construct2 || .model_error_cor[.indicator1, .indicator2] == 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Internal: Build DOT code for the SEM plot, including construct correlations.
#'
#' Constructs the DOT script for the SEM path diagram, now including
#' correlations between constructs (not just exogenous ones). Correctly
#' handles drawing only one edge per correlation.
#'
#' @param title The title of the plot.
#' @param graph_attrs Optional graph attributes.
#' @param constructs A vector of constructs.
#' @param r2_values Named vector of R2 values.
#' @param measurement_edge_fun Function to generate measurement edge code.
#' @param path_coefficients Matrix/data frame of path coefficients.
#' @param path_p_values Named vector of path p-values. Used for construct correlations too.
#' @param correlations List containing correlations (exogenous and indicator).
#' @param plot_significances Logical. Whether to display significance levels.
#' @param plot_correlations Option for indicator correlations ("none", "exo", or "all").
#' @param plot_structural_model_only Logical. Whether to display only the structural model.
#' @param is_second_order Logical. Whether the model is second-order.
#' @param model_measurement a matrix. The measurement matrix.
#' @param model_error_cor a matrix.
#' @param construct_correlations A matrix. The construct correlation matrix.
#' @param indicator_correlations A matrix. The indicator correlation matrix.
#'
#' @return A character string containing the complete DOT code.
#' @keywords internal
buildDotCode <- function(title, 
                         graph_attrs,
                         constructs, 
                         r2_values,
                         measurement_edge_fun,
                         path_coefficients, 
                         path_p_values,
                         correlations,
                         plot_significances,
                         plot_correlations, 
                         plot_structural_model_only,
                         plot_labels,
                         is_second_order = FALSE,
                         model_measurement = NULL,     # rows = constructs, columns = indicators
                         model_error_cor = NULL,
                         construct_correlations = NULL,  # Construct_VCV matrix (no p–values)
                         indicator_correlations = NULL  # Indicator_VCV matrix (no p–values)
){
  dot_lines <- c("digraph SEM {")
  
  # Title and graph attributes
  if (title != "") {
    dot_lines <- c(dot_lines,
                   "labelloc=\"t\";",
                   paste0("label=<", title, " >;"),
                   "fontsize=20;",
                   "fontname=\"Helvetica\";")
  }
  if (!is.null(graph_attrs)) {
    for (attr in graph_attrs) {
      dot_lines <- c(dot_lines, paste0(attr, ";"))
    }
  }
  
  # Define nodes.
  # When plot_labels is TRUE, include the R² value as in the original code.
  # When FALSE, show only the construct name.
  for (construct in names(constructs)) {
    constrType <- constructs[[construct]]
    shape <- if (constrType == "Composite") "hexagon" else "circle"
    fixed_size <- "width=1.5, height=1.5, fixedsize=false"
    node_label <- if (plot_labels) {
      if (!is.na(r2_values[construct])) {
        r2_label <- paste0("R\u00b2 = ", round(r2_values[construct], 3))
        paste0(construct, "\\n", r2_label)
      } else {
        construct
      }
    } else {
      construct
    }
    dot_lines <- c(dot_lines,
                   paste0(construct, " [label=\"", node_label,
                          "\", shape=", shape, ", ", fixed_size, "];"))
    # In a full plot we normally add measurement edges.
    # Here even if plot_structural_model_only is TRUE, if indicator correlations are being plotted with option "all"
    # we want to add the measurement (indicator) edges so that indicator correlation edges are attached.
    if (!plot_structural_model_only || (plot_structural_model_only && plot_correlations == "all")) {
      dot_lines <- c(dot_lines, measurement_edge_fun(construct))
    }
  }
  
  # ----- For hierarchical models, attach subordinate latent constructs.
  # If the measurement matrix indicates that a latent construct (row) "owns" indicators that are themselves latent nodes (present in 'constructs')
  # then add a connecting edge so that the subordinate indicator nodes aren’t isolated.
  if (!is.null(model_measurement)) {
    for (construct in names(constructs)) {
      if (construct %in% rownames(model_measurement)) {
        # Get the “indicators” according to the measurement matrix
        subs <- colnames(model_measurement)[which(as.numeric(model_measurement[construct, ]) == 1)]
        # Filter to those that are also latent nodes (i.e. in our constructs list)
        subs <- subs[subs %in% names(constructs)]
        if (length(subs) > 0) {
          for (sub in subs) {
            # Add a solid edge
            dot_lines <- c(dot_lines, paste0(construct, " -> ", sub, " [style=solid, penwidth=1.5];"))
          }
        }
      }
    }
  }
  
  # Add structural model (path) edges.
  for (dependent in rownames(path_coefficients)) {
    depPaths <- path_coefficients[dependent, ]
    predictors <- names(depPaths)[which(depPaths != 0)]
    for (predictor in predictors) {
      cleaned_dependent <- cleanNode(dependent)
      cleaned_predictor <- cleanNode(predictor)
      if (is_second_order && (cleaned_dependent == cleaned_predictor)) next
      edge_label <- if(plot_labels) {
        coefficient <- round(depPaths[predictor], 3)
        stars <- if (plot_significances) get_significance_stars(path_p_values[paste(dependent, "~", predictor)]) else ""
        paste0(coefficient, stars)
      } else {
        ""
      }
      if (grepl("\\.", predictor)) {
        cleaned_predictor <- gsub("_temp", "", predictor)
        dot_lines <- c(dot_lines,
                       paste0("\"", cleaned_predictor, "\" [label=\"", cleaned_predictor,
                              "\", shape=diamond, width=1.5, height=1.5, fixedsize=false];"),
                       paste0("\"", cleaned_predictor, "\" -> ", cleaned_dependent,
                              " [label=<", edge_label, ">, style=dashed, penwidth=2.5, weight=10, minlen=2];"))
      } else {
        dot_lines <- c(dot_lines,
                       paste0(cleaned_predictor, " -> ", cleaned_dependent,
                              " [label=<", edge_label, ">, penwidth=2.5, weight=10, minlen=2];"))
      }
    }
  }
  
  # Add correlation edges based on plot_correlations option.
  if (plot_correlations == "none") {
    # Do nothing.
  } else if (plot_correlations == "exo") {
    for (i in seq_along(correlations$exo$names)) {
      label_text <- if(plot_labels) {
        corr <- round(correlations$exo$estimates[i], 3)
        stars <- if (plot_significances) get_significance_stars(correlations$exo$p_values[i]) else ""
        paste0("(", corr, stars, ")")
      } else {
        ""
      }
      names_split <- strsplit(correlations$exo$names[i], " ~~ ")[[1]]
      if (plot_structural_model_only) {
        dot_lines <- c(dot_lines,
                       paste0(names_split[1], " -> ", names_split[2],
                              " [label=<", label_text, ">, style=dashed, arrowhead=vee, arrowtail=vee, dir=both, penwidth=2.5, weight=10];"))
      } else {
        dot_lines <- c(dot_lines,
                       paste0(names_split[1], " -> ", names_split[2],
                              " [label=<", label_text, ">, style=dashed, arrowhead=vee, arrowtail=vee, dir=both, penwidth=2.5, weight=10];"))
      }
    }
  } else if (plot_correlations == "all") {
    for (i in seq_along(correlations$exo$names)) {
      label_text <- if(plot_labels) {
        corr <- round(correlations$exo$estimates[i], 3)
        stars <- if (plot_significances) get_significance_stars(correlations$exo$p_values[i]) else ""
        paste0("(", corr, stars, ")")
      } else {
        ""
      }
      names_split <- strsplit(correlations$exo$names[i], " ~~ ")[[1]]
      if (plot_structural_model_only) {
        dot_lines <- c(dot_lines,
                       paste0(names_split[1], " -> ", names_split[2],
                              " [label=<", label_text, ">, style=dashed, arrowhead=vee, arrowtail=vee, dir=both];"))
      } else {
        dot_lines <- c(dot_lines,
                       paste0(names_split[1], " -> ", names_split[2],
                              " [label=<", label_text, ">, style=dashed, arrowhead=vee, arrowtail=vee, dir=both, penwidth=1];"))
      }
    }
    for (i in seq_along(correlations$ind$names)) {
      label_text <- if(plot_labels) {
        corr <- round(correlations$ind$estimates[i], 3)
        stars <- if (plot_significances) get_significance_stars(correlations$ind$p_values[i]) else ""
        paste0("(", corr, stars, ")")
      } else {
        ""
      }
      names_split <- strsplit(correlations$ind$names[i], " ~~ ")[[1]]
      if (check_connection(names_split[1], names_split[2], model_measurement, model_error_cor)) {
        if (plot_structural_model_only) {
          dot_lines <- c(dot_lines,
                         paste0(names_split[1], " -> ", names_split[2],
                                " [label=<", label_text, ">, style=dashed, arrowhead=vee, arrowtail=vee, dir=both];"))
        } else {
          dot_lines <- c(dot_lines,
                         paste0(names_split[1], " -> ", names_split[2],
                                " [label=<", label_text, ">, style=dashed, arrowhead=vee, arrowtail=vee, dir=both, penwidth=1];"))
        }
      }
    }
  } else {
    stop("Invalid option for plot_correlations. Valid options are 'none', 'exo', or 'all'.")
  }
  
  dot_lines <- c(dot_lines, "}")
  dot_code <- paste(dot_lines, collapse = "\n")
  return(dot_code)
}