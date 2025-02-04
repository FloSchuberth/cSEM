#' Plot model
#' 
#' \lifecycle{experimental}
#'
#' The `plotModel()` function creates a plot of a cSEM model using the \link[DiagrammeR]{grViz} function of the  DiagrammeR package.
#' 
#' @usage plotModel(
#' .object, 
#' .title =  args_default()$.title,
#' .plot_significances = args_default()$.plot_significances,  
#' .plot_indicator_correlations = args_default()$.plot_indicator_correlations,
#' .plot_structural_model_only = args_default()$.plot_structural_model_only,
#' .graph_attrs = args_default()$.graph_attrs
#' ) 
#'
#' @inheritParams csem_arguments
#' 
#' 
#' @return A DiagrammeR graph object or a list of DiagrammeR graph objects in case of a multigroup analysis.
#' 
#' @seealso [csem()], [cSEMResults], \link[DiagrammeR]{grViz}
#' 
#' @example inst/examples/example_plotModel.R
#' 
#' @export
plotModel <- function(
    .object, 
    .title =  args_default()$.title,
    .plot_significances = args_default()$.plot_significances, 
    .plot_indicator_correlations = args_default()$.plot_indicator_correlations,
    .plot_structural_model_only = args_default()$.plot_structural_model_only,
    .graph_attrs = args_default()$.graph_attrs
) {
  UseMethod("plotModel")
}

#'@export
plotModel.cSEMResults_default <- function(
    .object, 
    .title = args_default()$.title,
    .plot_significances = args_default()$.plot_significances, 
    .plot_indicator_correlations = args_default()$.plot_indicator_correlations,
    .plot_structural_model_only = args_default()$.plot_structural_model_only,
    .graph_attrs = args_default()$.graph_attrs
){
  
  ## Install DiagrammeR if not already installed
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop2(
      "Package `DiagrammeR` required. Use `install.packages(\"DiagrammeR\")` and rerun.")
  }
  
  if (inherits(.object, "cSEMResults_multi")) {
    plots <- lapply(names(.object), function(group_name) {
      group_object <- .object[[group_name]]
      group_title <- if (.title == "") paste0("Group_", group_name) else paste0(.title, " Group_", group_name)
      plotModel(group_object, 
                .title = group_title, 
                .plot_significances = .plot_significances, 
                .plot_indicator_correlations = .plot_indicator_correlations,
                .plot_structural_model_only = .plot_structural_model_only,
                .graph_attrs = .graph_attrs)
    })
    names(plots) <- names(.object)
    class(plots) <- c("cSEMPlot_multi", class(plots))
    return(plots)
    
  } else {
    results <- summarize(.object)
    constructs <- .object$Information$Model$construct_type  # named vector of construct types
    r2_values <- results$Estimates$R2
    weights <- as.data.frame(.object$Estimates$Weight_estimates)
    loadings <- as.data.frame(.object$Estimates$Loading_estimates)
    weight_p_values <- results$Estimates$Weight_estimates$p_value
    names(weight_p_values) <- results$Estimates$Weight_estimates$Name
    loading_p_values <- results$Estimates$Loading_estimates$p_value
    names(loading_p_values) <- results$Estimates$Loading_estimates$Name
    path_coefficients <- as.data.frame(.object$Estimates$Path_estimates)
    path_p_values <- results$Estimates$Path_estimates$p_value
    ind_corr <- list(names = results$Estimates$Indicator_correlation$Name,
                     estimates = results$Estimates$Indicator_correlation$Estimate,
                     p_values = results$Estimates$Indicator_correlation$p_value)
    exo_corr <- list(names = results$Estimates$Exo_construct_correlation$Name,
                     estimates = results$Estimates$Exo_construct_correlation$Estimate,
                     p_values = results$Estimates$Exo_construct_correlation$p_value)
    correlations <- list(ind = ind_corr, exo = exo_corr)
    
    measurement_edge_fun <- function(construct) {
      firstOrderMeasurementEdges(construct, weights, loadings, weight_p_values, loading_p_values, .plot_significances, constructs)
    }
    
    dot_code <- buildDotCode(title = .title,
                             graph_attrs = .graph_attrs,
                             constructs = constructs,
                             r2_values = r2_values,
                             measurement_edge_fun = measurement_edge_fun,
                             path_coefficients = path_coefficients,
                             path_p_values = path_p_values,
                             correlations = correlations,
                             plot_significances = .plot_significances,
                             plot_indicator_correlations = .plot_indicator_correlations,
                             plot_structural_model_only = .plot_structural_model_only,
                             is_second_order = FALSE)
    
    return(DiagrammeR::grViz(dot_code))
  }
}

#' @export
plotModel.cSEMResults <- function(
    .object, 
    .title = args_default()$.title,
    .plot_significances = args_default()$.plot_significances, 
    .plot_indicator_correlations = args_default()$.plot_indicator_correlations,
    .plot_structural_model_only = args_default()$.plot_structural_model_only,
    .graph_attrs = args_default()$.graph_attrs
) {
  
  # If the model contains second-stage information, treat it as a second‐order model;
  # otherwise call the default (first‐order) plotModel.
  if (!is.null(.object$Second_stage)) {
    plotModel.cSEMResults_2ndorder(
      .object,
      .title = .title,
      .plot_significances = .plot_significances,
      .plot_indicator_correlations = .plot_indicator_correlations,
      .plot_structural_model_only = .plot_structural_model_only,
      .graph_attrs = .graph_attrs)
  } else {
    plotModel.cSEMResults_default(
      .object,
      .title = .title,
      .plot_significances = .plot_significances,
      .plot_indicator_correlations = .plot_indicator_correlations,
      .plot_structural_model_only = .plot_structural_model_only,
      .graph_attrs = .graph_attrs)
  }
}

#' @export
plotModel.cSEMResults_2ndorder <- function(
    .object,
    .title = args_default()$.title,
    .plot_significances = args_default()$.plot_significances,
    .plot_indicator_correlations = args_default()$.plot_indicator_correlations,
    .plot_structural_model_only = args_default()$.plot_structural_model_only,
    .graph_attrs = args_default()$.graph_attrs
){
  
  ## Handle multi-group objects:
  if (inherits(.object, "cSEMResults_multi")) {
    plots <- lapply(names(.object), function(group_name) {
      group_object <- .object[[group_name]] # Extract results for each group
      group_title <- if (.title == "") paste0("Group_", group_name) else paste0(.title, " Group_", group_name)
      plotModel(group_object, # Recursively call plotModel for each group (now a single-group 2nd-order object)
                .title = group_title,
                .plot_significances = .plot_significances,
                .plot_indicator_correlations = .plot_indicator_correlations,
                .plot_structural_model_only = .plot_structural_model_only,
                .graph_attrs = .graph_attrs)
    })
    names(plots) <- names(.object)
    class(plots) <- c("cSEMPlot_multi", class(plots))
    return(plots)
    
  } else { # --- Single-group case ---
    
    # Extract first– and second–stage models and summaries.
    fs <- .object$First_stage
    ss <- .object$Second_stage
    results_fs <- summarize(.object)$First_stage
    results_ss <- summarize(.object)$Second_stage
    
    # Merge construct types from first– and second–stage.
    ct_first <- fs$Information$Model$construct_type
    ct_second <- ss$Information$Model$construct_type
    names(ct_second) <- gsub("_temp$", "", names(ct_second))
    constructs <- c(ct_second, ct_first)
    constructs <- constructs[!duplicated(names(constructs))]
    
    # R2 values from second–stage (clean names).
    r2_values <- results_ss$Estimates$R2
    names(r2_values) <- gsub("_temp$", "", names(r2_values))
    
    # First–stage measurement parameters.
    weights_fs <- as.data.frame(fs$Estimates$Weight_estimates)
    loadings_fs <- as.data.frame(fs$Estimates$Loading_estimates)
    weight_p_fs <- results_fs$Estimates$Weight_estimates$p_value
    names(weight_p_fs) <- results_fs$Estimates$Weight_estimates$Name
    loading_p_fs <- results_fs$Estimates$Loading_estimates$p_value
    names(loading_p_fs) <- results_fs$Estimates$Loading_estimates$Name
    
    # Second–stage measurement parameters (clean row and col names).
    weights_ss <- as.data.frame(ss$Estimates$Weight_estimates)
    rownames(weights_ss) <- gsub("_temp$", "", rownames(weights_ss))
    colnames(weights_ss) <- gsub("_temp$", "", colnames(weights_ss))
    weight_p_ss <- results_ss$Estimates$Weight_estimates$p_value
    
    # Structural model paths from second–stage.
    path_ss <- as.data.frame(ss$Estimates$Path_estimates)
    rownames(path_ss) <- gsub("_temp$", "", rownames(path_ss))
    colnames(path_ss) <- gsub("_temp$", "", colnames(path_ss))
    path_p_ss <- results_ss$Estimates$Path_estimates$p_value
    
    # --- Structural paths: if .plot_structural_model_only is TRUE, use only second-stage edges.
    if (.plot_structural_model_only) {
      combined_path_coeff <- path_ss
      combined_path_p <- path_p_ss
    } else {
      # Merge second-stage and non–duplicate first-stage edges.
      path_fs <- as.data.frame(fs$Estimates$Path_estimates)
      path_p_fs <- results_fs$Estimates$Path_estimates$p_value
      names(path_p_fs) <- results_fs$Estimates$Path_estimates$Name
      combined_path_coeff <- path_ss
      combined_path_p <- path_p_ss
      for (dependent in rownames(path_fs)) {
        fs_paths <- path_fs[dependent, ]
        predictors <- names(fs_paths)[which(fs_paths != 0)]
        for (predictor in predictors) {
          found <- FALSE
          if (dependent %in% rownames(combined_path_coeff)) {
            if (predictor %in% colnames(combined_path_coeff)) {
              if (isTRUE(combined_path_coeff[dependent, predictor] != 0)) found <- TRUE
            }
          }
          if (!found) {
            combined_path_coeff[dependent, predictor] <- fs_paths[predictor]
            path_name <- paste(dependent, "~", predictor)
            combined_path_p[path_name] <- path_p_fs[path_name]
          }
        }
      }
    }
    # --- End structural path merging.
    
    # Correlations: exogenous from 2nd-stage; indicator from first-stage.
    exo_corr <- list(names = gsub("_temp$", "", results_ss$Estimates$Exo_construct_correlation$Name),
                     estimates = results_ss$Estimates$Exo_construct_correlation$Estimate,
                     p_values = results_ss$Estimates$Exo_construct_correlation$p_value)
    ind_corr <- list(names = results_fs$Estimates$Indicator_correlation$Name,
                     estimates = results_fs$Estimates$Indicator_correlation$Estimate,
                     p_values = results_fs$Estimates$Indicator_correlation$p_value)
    correlations <- list(exo = exo_corr, ind = ind_corr)
    
    # Define measurement edge function for second–order.
    # When .plot_structural_model_only is TRUE, we set only_second_stage = TRUE so that only second-stage edges are drawn.
    if (.plot_structural_model_only) {
      measurement_edge_fun <- function(construct) {
        secondOrderMeasurementEdges(construct,
                                    weights_first = weights_fs,
                                    loadings_first = loadings_fs,
                                    weight_p_first = weight_p_fs,
                                    loading_p_first = loading_p_fs,
                                    weights_second = weights_ss,
                                    weight_p_second = weight_p_ss,
                                    plot_signif = .plot_significances,
                                    constructTypes = constructs,
                                    only_second_stage = TRUE)
      }
    } else {
      measurement_edge_fun <- function(construct) {
        secondOrderMeasurementEdges(construct,
                                    weights_first = weights_fs,
                                    loadings_first = loadings_fs,
                                    weight_p_first = weight_p_fs,
                                    loading_p_first = loading_p_fs,
                                    weights_second = weights_ss,
                                    weight_p_second = weight_p_ss,
                                    plot_signif = .plot_significances,
                                    constructTypes = constructs,
                                    only_second_stage = FALSE)
      }
    }
    
    # Build the DOT script.
    # Pass is_second_order = TRUE so that buildDotCode always calls the measurement_edge_fun.
    dot_code <- buildDotCode(title = .title,
                             graph_attrs = .graph_attrs,
                             constructs = constructs,
                             r2_values = r2_values,
                             measurement_edge_fun = measurement_edge_fun,
                             path_coefficients = combined_path_coeff,
                             path_p_values = combined_path_p,
                             correlations = correlations,
                             plot_significances = .plot_significances,
                             plot_indicator_correlations = .plot_indicator_correlations,
                             plot_structural_model_only = .plot_structural_model_only,
                             is_second_order = TRUE)
    
    return(DiagrammeR::grViz(dot_code))
  }
}