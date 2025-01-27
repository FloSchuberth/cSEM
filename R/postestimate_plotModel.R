#' Plot model
#' 
#' \lifecycle{experimental}
#'
#' Create a plot of a cSEM model
#' 
#' @usage plotModel(
#' .object, 
#' .title = "",
#' .plot_significances = args_default()$.plot_significances,  
#' .plot_indicator_correlations = args_default()$.plot_indicator_correlations,
#' .plot_structural_model_only = args_default()$.plot_structural_model_only,
#' .plot_stage_SOC = c("second", "first"),
#' .node_ranks = NULL, 
#' .graph_attrs = NULL,
#' .remove_colors = FALSE, 
#' .edge_lengths = NULL
#' ) 
#'
#' @inheritParams csem_arguments
#' @param .title Character string. Plot title.
#' @param .node_ranks List specifying the ranks for nodes (e.g., list(rank1 = c("Node1", "Node2"))), DiagrammeR syntax.
#' @param .graph_attrs Character vector of custom graph attributes for the DOT graph (e.g., c("rankdir=LR", "ranksep=1.0")). DiagrammeR syntax.
#' @param .remove_colors Logical. Should the plot be colored? Defaults to `FALSE`. 
#' @param .edge_lengths Integer. Length of edges in the DOT graph.
#' @param .plot_stage_SOC Character string. Which stage of the two-stage approach for models containing second-order constructs should be displayed? 
#' Defaults to `second`.
#' 
#' 
#' @return A DiagrammeR graph object or a list of DiagrammeR graph objects in case of a multi analysis.
#' 
#' @seealso [csem()], [cSEMResults]
#' 
#' @example inst/examples/example_plotModel.R
#' 
#' @export
plotModel <- function(
    .object, 
    .title = "",
    .plot_significances = args_default()$.plot_significances, 
    .plot_indicator_correlations = args_default()$.plot_indicator_correlations,
    .plot_structural_model_only = args_default()$.plot_structural_model_only,
    .plot_stage_SOC = c("second", "first"),
    .node_ranks = NULL, 
    .graph_attrs = NULL,
    .remove_colors = FALSE, 
    .edge_lengths = NULL
                      ) {
  
  # Match the stage argument
  .plot_stage_SOC <- match.arg(.plot_stage_SOC)
  
  # Check if the object is a multi-group object
  if (inherits(.object, "cSEMResults_multi")) {
    plots <- lapply(names(.object), function(group_name) {
      
      group_object <- .object[[group_name]]
      
      if (.title == "") {
        group_title <- paste0("Group_", group_name)
      } else {
        group_title <- paste0(.title, " Group_", group_name)
      }
      plotModel(group_object, .plot_significances, .title = group_title,
                .plot_indicator_correlations,.plot_structural_model_only,
                .node_ranks, .graph_attrs, .remove_colors, .edge_lengths, .plot_stage_SOC)
    })
    
    # name each plot for later use
    names(plots) <- names(.object)
    
    # Assign the list of plots to the global variable latest_plot
    class(plots) <- c("cSEMPlot_multi", class(plots))
    # assign("latest_plot", plots, envir = .GlobalEnv)
    
    return(plots)
    
  } else {
    # Original code for single-group and multi-group, but also second-order
    # models is used (after minor adaptations)
    # Check if the object is a second-order model and perform plot generation accordingly
    if(inherits(.object, "cSEMResults_2ndorder")){
      if(.plot_stage_SOC == "second"){
        # Ensure that the second-stage path estimates exist
        if (is.null(.object$Second_stage$Estimates$Path_estimates)) stop2("Second-stage path estimates could not be found.")
        
        # Summarize results to get estimates and p-values
        results <- summarize(.object)$Second_stage
        
        # Extract estimates and model information, and remove "_temp" if present
        path_coefficients <- as.data.frame(.object$Second_stage$Estimates$Path_estimates)
        rownames(path_coefficients) <- gsub("_temp$", "", rownames(path_coefficients))
        colnames(path_coefficients) <- gsub("_temp$", "", colnames(path_coefficients))
        
        weights <- as.data.frame(.object$Second_stage$Estimates$Weight_estimates)
        rownames(weights) <- gsub("_temp$", "", rownames(weights))
        
        loadings <- as.data.frame(.object$Second_stage$Estimates$Loading_estimates)
        rownames(loadings) <- gsub("_temp$", "", rownames(loadings))
        
        construct_types <- results$Information$Model$construct_type
        
        # Extract R2 values
        r2_values <- results$Estimates$R2
        
        # Extract p-values for path, weight, and loading estimates, remove "_temp"
        path_p_values <- results$Estimates$Path_estimates$p_value
        names(path_p_values) <- gsub("_temp","",results$Estimates$Path_estimates$Name)
        
        weight_p_values <- results$Estimates$Weight_estimates$p_value
        names(weight_p_values) <- gsub("_temp","",results$Estimates$Weight_estimates$Name)
        
        loading_p_values <- results$Estimates$Loading_estimates$p_value
        names(loading_p_values) <- gsub("_temp","",results$Estimates$Loading_estimates$Name)
        
        # Extract correlations and estimates for indicators and exogenous constructs, remove "_temp" if present
        indicator_correlation_names <- gsub("_temp","",results$Estimates$Indicator_correlation$Name)
        indicator_correlation_estimates <- results$Estimates$Indicator_correlation$Estimate
        exo_construct_correlation_names <- gsub("_temp","",results$Estimates$Exo_construct_correlation$Name)
        exo_construct_correlation_estimates <- results$Estimates$Exo_construct_correlation$Estimate
        
        # Extract p-values for indicator and exogenous construct correlations
        indicator_correlation_p_values <- results$Estimates$Indicator_correlation$p_value
        exo_construct_correlation_p_values <- results$Estimates$Exo_construct_correlation$p_value
        
      }else{ # first stage
        # Ensure that the first-stage path estimates exist
        if (is.null(.object$First_stage$Estimates$Path_estimates)) stop("First-stage path estimates not found in the model.")
        
        # Summarize results to get estimates and p-values
        results <- summarize(.object)$First_stage
        
        # Extract estimates and model information
        path_coefficients <- as.data.frame(.object$First_stage$Estimates$Path_estimates)
        weights <- as.data.frame(.object$First_stage$Estimates$Weight_estimates)
        loadings <- as.data.frame(.object$First_stage$Estimates$Loading_estimates)
        construct_types <- .object$First_stage$Information$Model$construct_type
        
        # Extract R2
        r2_values <- results$Estimates$R2
        
        # Extract p-values for path, weight, and loading estimates
        path_p_values <- results$Estimates$Path_estimates$p_value
        names(path_p_values) <- results$Estimates$Path_estimates$Name
        
        weight_p_values <- results$Estimates$Weight_estimates$p_value
        names(weight_p_values) <- results$Estimates$Weight_estimates$Name
        
        loading_p_values <- results$Estimates$Loading_estimates$p_value
        names(loading_p_values) <- results$Estimates$Loading_estimates$Name
        
        # Extract correlations and estimates for indicators and exogenous constructs
        indicator_correlation_names <- results$Estimates$Indicator_correlation$Name
        indicator_correlation_estimates <- results$Estimates$Indicator_correlation$Estimate
        exo_construct_correlation_names <- results$Estimates$Exo_construct_correlation$Name
        exo_construct_correlation_estimates <- results$Estimates$Exo_construct_correlation$Estimate
        
        # Extract p-values for indicator and exogenous construct correlations
        indicator_correlation_p_values <- results$Estimates$Indicator_correlation$p_value
        exo_construct_correlation_p_values <- results$Estimates$Exo_construct_correlation$p_value
      }
      
    } else {
      # Ensure that the path estimates exist
      if (is.null(.object$Estimates$Path_estimates)) stop2("Path estimates could not be found.")
      
      # Summarize results to get estimates and p-values
      results <- summarize(.object)
      
      # Extract estimates and model information
      path_coefficients <- as.data.frame(.object$Estimates$Path_estimates)
      weights <- as.data.frame(.object$Estimates$Weight_estimates)
      loadings <- as.data.frame(.object$Estimates$Loading_estimates)
      construct_types <- .object$Information$Model$construct_type
      
      # Extract R2
      r2_values <- results$Estimates$R2
      
      # Extract p-values for path, weight, and loading estimates
      path_p_values <- results$Estimates$Path_estimates$p_value
      names(path_p_values) <- results$Estimates$Path_estimates$Name
      
      weight_p_values <- results$Estimates$Weight_estimates$p_value
      names(weight_p_values) <- results$Estimates$Weight_estimates$Name
      
      loading_p_values <- results$Estimates$Loading_estimates$p_value
      names(loading_p_values) <- results$Estimates$Loading_estimates$Name
      
      # Extract correlations and estimates for indicators and exogenous constructs
      indicator_correlation_names <- results$Estimates$Indicator_correlation$Name
      indicator_correlation_estimates <- results$Estimates$Indicator_correlation$Estimate
      exo_construct_correlation_names <- results$Estimates$Exo_construct_correlation$Name
      exo_construct_correlation_estimates <- results$Estimates$Exo_construct_correlation$Estimate
      
      # Extract p-values for indicator and exogenous construct correlations
      indicator_correlation_p_values <- results$Estimates$Indicator_correlation$p_value
      exo_construct_correlation_p_values <- results$Estimates$Exo_construct_correlation$p_value
    }
    

    # Initialize DOT code for Graphviz
    dot_code <- "digraph SEM {\n"
    
    # Add title to the graph if provided
    if (.title != "") {
      dot_code <- paste0(dot_code,
                         "labelloc=\"t\";\n",
                         "label=<", .title, " >;\n",
                         "fontsize=20;\nfontname=\"Helvetica\";\n")
    }
    
    # Add rank constraints for nodes if provided
    if (!is.null(.node_ranks)) {
      for (rank_name in names(.node_ranks)) {
        nodes_in_rank <- .node_ranks[[rank_name]]
        dot_code <- paste0(dot_code,
                           "{rank=same; ",
                           paste(nodes_in_rank, collapse = "; "),
                           ";}\n")
      }
    }
    
    # Add custom graph attributes (e.g., rankdir, ranksep, nodesep)
    if (!is.null(.graph_attrs)) {
      for (attr in .graph_attrs) {
        dot_code <- paste0(dot_code, attr, ";\n")
      }
    }
    
    # Define constructs with R2 values and subgraphs for their indicators
    constructs <- names(construct_types)
    
    # Iterate over constructs
    for (construct in constructs) {
      # Get construct type and assign color and shape
      construct_type <- construct_types[construct]
      color <- if (construct_type == "Composite") "orange" else "lightblue"
      shape <- if (construct_type == "Composite") "hexagon" else "circle"
      fixed_size <- "width=1.5, height=1.5, fixedsize=true"
      
      # Add R2 value to label if available
      if (!is.na(r2_values[construct])) {
        r2_label <- paste0("R\u00b2 = ", round(r2_values[construct], 3))
        dot_code <- paste0(dot_code, construct,
                           " [label=\"", construct, "\\n", r2_label,
                           "\", shape=", shape, ", style=filled, color=", color, ", ", fixed_size, "];\n")
      } else {
        dot_code <- paste0(dot_code, construct,
                           " [label=\"", construct,
                           "\", shape=", shape, ", style=filled, color=", color, ", ", fixed_size, "];\n")
      }
      
      # Create subgraph for each construct's indicators, unless showing structural model only
      if (!.plot_structural_model_only) {
        dot_code <- paste0(dot_code, "subgraph cluster_", construct, " {\n",
                           "label=\"\";\n",
                           "style=invis;\n")
        
        # Get indicators with non-zero weights
        non_zero_indices <- which(weights[construct, ] != 0)
        non_zero_indicators <- colnames(weights)[non_zero_indices]
        # Ensure that non_zero_indicators are not NULL
        if (length(non_zero_indicators) > 0) {
          for (indicator in non_zero_indicators) {
            # Define indicator node
            dot_code <- paste0(dot_code,
                               indicator,
                               " [shape=box, style=filled, color=lightgray];\n")
            
            # Get weight and loading for the indicator
            weight <- round(weights[construct, indicator], 3)
            loading <- round(loadings[construct, indicator], 3)
            weight_name <- paste(construct, "<~", indicator)
            loading_name <- paste(construct, "=~", indicator)
            weight_stars <- if (.plot_significances) get_significance_stars(weight_p_values[weight_name]) else ""
            loading_stars <- if (.plot_significances) get_significance_stars(loading_p_values[loading_name]) else ""
            
            # Create label for edge
            label <- if (construct_type == "Common factor") {
              paste0(loading, loading_stars)
            } else {
              paste0(weight, weight_stars)
            }
            
            # Determine direction based on construct type
            direction <- if (construct_type == "Common factor") {
              paste0(construct, " -> ", indicator)
            } else {
              paste0(indicator, " -> ", construct)
            }
            # Add edge length attribute if specified
            len_attr <- if (!is.null(.edge_lengths)) paste0(", len=", .edge_lengths) else ""
            
            # Add edge from construct to indicator
            dot_code <- paste0(dot_code,
                               direction,
                               " [label=\"", label, "\", color=black", len_attr, "];\n")
          }
        }
        dot_code <- paste0(dot_code, "}\n") # Close subgraph
      }
    }
    
    # Add edges for structural model paths
    for (dependent in rownames(path_coefficients)) {
      dependent_paths <- path_coefficients[dependent, ]
      non_zero_predictors <- names(dependent_paths[which(dependent_paths != 0)])
      
      for (predictor in non_zero_predictors) {
        coefficient <- round(dependent_paths[predictor], 3)
        path_name <- paste(dependent, "~", predictor)
        
        stars <- if (.plot_significances) get_significance_stars(path_p_values[path_name]) else ""
        
        if (grepl("\\.", predictor)) {
          # Create a new node for the predictor with dot in its name
          dot_code <- paste0(dot_code,
                             "\"", predictor, "\"",
                             " [label=\"", predictor, "\", shape=diamond, style=filled, color=bisque, width=1.5, height=1.5, fixedsize=true];\n")
          
          # Create edge from the new predictor node to the dependent variable
          dot_code <- paste0(dot_code,
                             "\"", predictor, "\"", " -> ", dependent,
                             " [label=<", coefficient, stars, ">",
                             ", color=blue, style=dashed, penwidth=2.5, fontcolor=brown];\n")
          
        } else {
          # Construct the path string for predictors without dots
          dot_code <- paste0(dot_code,
                             predictor, " -> ", dependent,
                             " [label=<", coefficient, stars, ">",
                             ", color=blue, penwidth=2.5, fontcolor=brown];\n")
        }
      }
    }
    
    # Add dashed lines for exogenous construct correlations
    for (i in seq_along(exo_construct_correlation_names)) {
      correlation <- round(exo_construct_correlation_estimates[i], 3)
      stars <- if (.plot_significances) get_significance_stars(exo_construct_correlation_p_values[i]) else ""
      names_split <- strsplit(exo_construct_correlation_names[i], " ~~ ")[[1]]
      dot_code <- paste0(dot_code,
                         names_split[1], " -> ", names_split[2],
                         " [label=<(", correlation, stars,
                         ")>", ", fontcolor=darkkhaki, style=dashed, dir=none, color=gray50];\n")
    }
    
    # Add dashed lines for indicator correlations if `run_correlation_lines` is TRUE
    if (.plot_indicator_correlations) {
      for (i in seq_along(indicator_correlation_names)) {
        correlation <- round(indicator_correlation_estimates[i], 3)
        stars <- if (.plot_significances) get_significance_stars(indicator_correlation_p_values[i]) else ""
        names_split <- strsplit(indicator_correlation_names[i], " ~~ ")[[1]]
        dot_code <- paste0(dot_code,
                           names_split[1], " -> ", names_split[2],
                           " [label=<(", correlation, stars,
                           ")>", ", fontcolor=darkkhaki, style=dashed, dir=none, color=gray50];\n")
      }
    }
    
    # Remove specific colors from DOT code if `remove_colors` is TRUE
    if (.remove_colors) {
      dot_code <- gsub("color=darkkhaki", "color=black", dot_code)
      dot_code <- gsub("color=blue", "color=black", dot_code)
      dot_code <- gsub("color=brown", "color=black", dot_code)
      dot_code <- gsub("color=orange", "fillcolor=white, color=black", dot_code)
      dot_code <- gsub("color=lightblue", "fillcolor=white, color=black", dot_code)
      dot_code <- gsub("color=bisque", "fillcolor=white, color=black", dot_code)
    }
    
    # Close the DOT graph
    dot_code <- paste0(dot_code, "}")
    
    # Visualize using DiagrammeR
    latest_plot <- DiagrammeR::grViz(dot_code)
    
    # Assign latest_plot to the global environment
    # assign("latest_plot", latest_plot, envir = .GlobalEnv)
    
    # Return the plot object
    latest_plot
  }
}
