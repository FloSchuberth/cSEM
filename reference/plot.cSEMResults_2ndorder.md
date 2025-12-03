# `cSEMResults` method for `plot()` for second-order models.

**\[experimental\]**

## Usage

``` r
# S3 method for class 'cSEMResults_2ndorder'
plot(
  x,
  .title = args_default()$.title,
  .plot_significances = args_default()$.plot_significances,
  .plot_correlations = args_default()$.plot_correlations,
  .plot_structural_model_only = args_default()$.plot_structural_model_only,
  .plot_labels = args_default()$.plot_labels,
  .graph_attrs = args_default()$.graph_attrs,
  ...
)
```

## Arguments

- x:

  An R object of class `cSEMResults_2ndorder` object.

- .title:

  Character string. Title of an object. Defaults to *""*.

- .plot_significances:

  Logical. Should p-values in the form of stars be plotted? Defaults to
  `TRUE`.

- .plot_correlations:

  Character string. Specify which correlations should be plotted, i.e.,
  between the exogenous constructs (`exo`), between the exogenous
  constructs and the indicators (`all`), or not at all (`none`).
  Defaults to `exo`.

- .plot_structural_model_only:

  Logical. Should only the structural model, i.e., the constructs and
  their relationships be plotted? Defaults to `FALSE`.

- .plot_labels:

  Logical. Whether to display edge labels and node R² values. Defaults
  to TRUE.

- .graph_attrs:

  Character string. Additional attributes that should be passed to the
  DiagrammeR syntax, e.g., c("rankdir=LR", "ranksep=1.0"). Defaults to
  *c("rankdir=LR")*.

- ...:

  Currently ignored.

## Details

Creates a plot of a `cSEMResults_2ndorder` object using the
[grViz](https://rich-iannone.github.io/DiagrammeR/reference/grViz.html)
function. For more details on customizing plot, see
<https://rpubs.com/nguyen_mot/1275413>.

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md),
[grViz](https://rich-iannone.github.io/DiagrammeR/reference/grViz.html)
