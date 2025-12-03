# Internal: Build DOT code for the SEM plot, including construct correlations.

Constructs the DOT script for the SEM path diagram, now including
correlations between constructs (not just exogenous ones). Correctly
handles drawing only one edge per correlation.

## Usage

``` r
buildDotCode(
  title,
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
  model_measurement = NULL,
  model_error_cor = NULL,
  construct_correlations = NULL,
  indicator_correlations = NULL
)
```

## Arguments

- title:

  The title of the plot.

- graph_attrs:

  Optional graph attributes.

- constructs:

  A vector of constructs.

- r2_values:

  Named vector of R2 values.

- measurement_edge_fun:

  Function to generate measurement edge code.

- path_coefficients:

  Matrix/data frame of path coefficients.

- path_p_values:

  Named vector of path p-values. Used for construct correlations too.

- correlations:

  List containing correlations (exogenous and indicator).

- plot_significances:

  Logical. Whether to display significance levels.

- plot_correlations:

  Option for indicator correlations ("none", "exo", or "all").

- plot_structural_model_only:

  Logical. Whether to display only the structural model.

- is_second_order:

  Logical. Whether the model is second-order.

- model_measurement:

  a matrix. The measurement matrix.

- model_error_cor:

  a matrix.

- construct_correlations:

  A matrix. The construct correlation matrix.

- indicator_correlations:

  A matrix. The indicator correlation matrix.

## Value

A character string containing the complete DOT code.
