#' `cSEMResults` method for `plot()`
#' 
#' \lifecycle{experimental}
#'
#' Creates a plot of a `cSEMResults` object using the \link[DiagrammeR]{grViz} function.
#' 
#' @inheritParams csem_arguments
#' @param x An R object of class `cSEMResults_default` object.
#' @param ... Currently ignored.
#' @param .plot_labels Logical. Whether to display edge labels and R² values in the nodes.
#'        Defaults to TRUE (i.e. original plot).
#' 
#' @seealso [savePlot()] [csem()], [cSEMResults], \link[DiagrammeR]{grViz}
#' @example inst/examples/example_plot.cSEMResults.R
#' @export
plot.cSEMResults_default <- function(
    x = NULL, 
    .title = args_default()$.title,
    .plot_significances = args_default()$.plot_significances, 
    .plot_correlations = args_default()$.plot_correlations,
    .plot_structural_model_only = args_default()$.plot_structural_model_only,
    .plot_labels = args_default()$.plot_labels,
    .graph_attrs = args_default()$.graph_attrs,
    ...
){
  
  ## Install DiagrammeR if not already installed  
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop2(
      "Package `DiagrammeR` is required. Use `install.packages(\"DiagrammeR\")` and rerun.")
  }
  
  results <- summarize(x)
  constructs <- x$Information$Model$construct_type  # named vector of construct types
  r2_values <- results$Estimates$R2
  weights <- as.data.frame(x$Estimates$Weight_estimates)
  loadings <- as.data.frame(x$Estimates$Loading_estimates)
  weight_p_values <- results$Estimates$Weight_estimates$p_value
  names(weight_p_values) <- results$Estimates$Weight_estimates$Name
  loading_p_values <- results$Estimates$Loading_estimates$p_value
  names(loading_p_values) <- results$Estimates$Loading_estimates$Name
  path_coefficients <- as.data.frame(x$Estimates$Path_estimates)
  path_p_values <- results$Estimates$Path_estimates$p_value
  ind_corr <- list(names = results$Estimates$Indicator_correlation$Name,
                   estimates = results$Estimates$Indicator_correlation$Estimate,
                   p_values = results$Estimates$Indicator_correlation$p_value)
  exo_corr <- list(names = results$Estimates$Exo_construct_correlation$Name,
                   estimates = results$Estimates$Exo_construct_correlation$Estimate,
                   p_values = results$Estimates$Exo_construct_correlation$p_value)
  correlations <- list(ind = ind_corr, exo = exo_corr)
  
  measurement_edge_fun <- function(construct) {
    firstOrderMeasurementEdges(construct, weights, loadings, weight_p_values, loading_p_values, .plot_significances, constructs, plot_labels = .plot_labels)
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
                           plot_correlations = .plot_correlations,
                           plot_structural_model_only = .plot_structural_model_only,
                           is_second_order = FALSE,
                           model_measurement = x$Information$Model$measurement,
                           model_error_cor = x$Information$Model$error_cor,
                           construct_correlations = x$Estimates$Construct_VCV,
                           indicator_correlations = x$Estimates$Indicator_VCV,
                           plot_labels = .plot_labels)
  
  out <- DiagrammeR::grViz(dot_code)
  class(out) <- c(class(out), "cSEMPlot_single")
  return(out)
}

#' `cSEMResults` method for `plot()` for multiple groups.
#' 
#' \lifecycle{experimental}
#'
#' Creates a plot of a `cSEMResults` object using the \link[DiagrammeR]{grViz} function.
#' 
#' @inheritParams csem_arguments
#' @param x An R object of class `cSEMResults_multi` object.
#' @param ... Currently ignored.
#' @param .plot_labels Logical. Whether to display edge labels and node R² values. Defaults to TRUE.
#' @seealso [csem()], [cSEMResults], \link[DiagrammeR]{grViz}
#' @export
plot.cSEMResults_multi <- function(
    x = NULL, 
    .title =  args_default()$.title,
    .plot_significances = args_default()$.plot_significances, 
    .plot_correlations = args_default()$.plot_correlations,
    .plot_structural_model_only = args_default()$.plot_structural_model_only,
    .plot_labels = args_default()$.plot_labels,
    .graph_attrs = args_default()$.graph_attrs,
    ...
){
  plots <- Map(function(group_name, group_object) {
    group_title <- if (.title == "") paste0("Group_", group_name) else paste0(.title, " Group_", group_name)
    plot(group_object,
         .title = group_title,
         .plot_significances = .plot_significances,
         .plot_correlations = .plot_correlations,
         .plot_structural_model_only = .plot_structural_model_only,
         .plot_labels = .plot_labels,
         .graph_attrs = .graph_attrs)
  }, names(x), x)
  
  class(plots) <- c("cSEMPlot_multi", class(plots))
  return(plots)
}

#' `cSEMResults` method for `plot()` for second-order models.
#'
#' \lifecycle{experimental}
#'
#' Creates a plot of a `cSEMResults_2ndorder` object using the \link[DiagrammeR]{grViz} function.
#' 
#' @inheritParams csem_arguments
#' @param x An R object of class `cSEMResults_2ndorder` object.
#' @param ... Currently ignored.
#' @param .plot_labels Logical. Whether to display edge labels and node R² values. Defaults to TRUE.
#' @seealso [csem()], [cSEMResults], \link[DiagrammeR]{grViz}
#' @export
plot.cSEMResults_2ndorder <- function(
    x,
    .title = args_default()$.title,
    .plot_significances = args_default()$.plot_significances,
    .plot_correlations = args_default()$.plot_correlations,
    .plot_structural_model_only = args_default()$.plot_structural_model_only,
    .plot_labels = args_default()$.plot_labels,
    .graph_attrs = args_default()$.graph_attrs,
    ...
){
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop2("Package `DiagrammeR` is required. Use `install.packages(\"DiagrammeR\")` and rerun.")
  }
  # Extract first– and second–stage models and summaries.
  results <- summarize(x)
  results_fs <- results$First_stage
  results_ss <- results$Second_stage
  
  fs <- x$First_stage
  ss <- x$Second_stage
  
  # Merge construct types from first– and second–stage.
  ct_first <- fs$Information$Model$construct_type
  ct_second <- ss$Information$Model$construct_type
  names(ct_second) <- gsub("_temp$", "", names(ct_second))
  ct2_names <- names(ct_second)
  constructs <- c(ct_second, ct_first)
  constructs <- constructs[!duplicated(names(constructs))]
  
  if (.plot_structural_model_only && !(.plot_correlations %in% c("all"))) {
    constructs <- ct_second
  }
  
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
  
  loadings_ss <- as.data.frame(ss$Estimates$Loading_estimates)
  rownames(loadings_ss) <- gsub("_temp$", "", rownames(loadings_ss))
  colnames(loadings_ss) <- gsub("_temp$", "", colnames(loadings_ss))
  loading_p_ss <- results_ss$Estimates$Loading_estimates$p_value
  
  # Structural model paths from second–stage.
  path_ss <- as.data.frame(ss$Estimates$Path_estimates)
  rownames(path_ss) <- gsub("_temp$", "", rownames(path_ss))
  colnames(path_ss) <- gsub("_temp$", "", colnames(path_ss))
  path_p_ss <- results_ss$Estimates$Path_estimates$p_value
  names(path_p_ss) <- results_ss$Estimates$Path_estimates$Name
  
  combined_path_coeff <- path_ss
  combined_path_p <- path_p_ss
  
  exo_corr <- list(names = gsub("_temp$", "", results_ss$Estimates$Exo_construct_correlation$Name),
                   estimates = results_ss$Estimates$Exo_construct_correlation$Estimate,
                   p_values = results_ss$Estimates$Exo_construct_correlation$p_value)
  ind_corr <- list(names = results_fs$Estimates$Indicator_correlation$Name,
                   estimates = results_fs$Estimates$Indicator_correlation$Estimate,
                   p_values = results_fs$Estimates$Indicator_correlation$p_value)
  correlations <- list(exo = exo_corr, ind = ind_corr)
  
  # Only skip measurement edges when no indicator correlations should be added.
  measurement_edge_fun <- function(construct) {
    if (.plot_structural_model_only && (.plot_correlations != "all")) return("")
    if (construct %in% ct2_names) {
      return(secondOrderMeasurementEdges(construct,
                                         weights_first    = weights_fs,
                                         loadings_first   = loadings_fs,
                                         weight_p_first   = weight_p_fs,
                                         loading_p_first  = loading_p_fs,
                                         weights_second   = weights_ss,
                                         loadings_second  = loadings_ss,
                                         weight_p_second  = weight_p_ss,
                                         loading_p_second = loading_p_ss,
                                         plot_signif      = .plot_significances,
                                         plot_labels      = .plot_labels,
                                         constructTypes   = constructs,
                                         only_second_stage = FALSE))
    } else {
      return(firstOrderMeasurementEdges(construct,
                                        weights = weights_fs,
                                        loadings = loadings_fs,
                                        weight_p_values = weight_p_fs,
                                        loading_p_values = loading_p_fs,
                                        plot_signif = .plot_significances,
                                        plot_labels = .plot_labels,
                                        constructTypes = constructs))
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
                           plot_correlations = .plot_correlations,
                           plot_structural_model_only = .plot_structural_model_only,
                           plot_labels = .plot_labels,
                           is_second_order = TRUE,
                           model_measurement = x$First_stage$Information$Model$measurement,
                           model_error_cor = x$First_stage$Information$Model$error_cor,
                           construct_correlations = results_fs$Estimates$Construct_VCV,
                           indicator_correlations = fs$Estimates$Indicator_VCV)
  
  out <- DiagrammeR::grViz(dot_code)
  class(out) <- c(class(out), "cSEMPlot_single")
  return(out)
}