# `cSEMResults` method for `plot()`

**\[experimental\]**

## Usage

``` r
# S3 method for class 'cSEMResults_default'
plot(
  x = NULL,
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

  An R object of class `cSEMResults_default` object.

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

  Logical. Whether to display edge labels and R² values in the nodes.
  Defaults to TRUE (i.e. original plot).

- .graph_attrs:

  Character string. Additional attributes that should be passed to the
  DiagrammeR syntax, e.g., c("rankdir=LR", "ranksep=1.0"). Defaults to
  *c("rankdir=LR")*.

- ...:

  Currently ignored.

## Details

Creates a plot of a `cSEMResults` object using the
[grViz](https://rich-iannone.github.io/DiagrammeR/reference/grViz.html)
function. For more details on customizing plot, see
<https://rpubs.com/nguyen_mot/1275413>.

## See also

[`savePlot()`](https://floschuberth.github.io/cSEM/reference/savePlot.md)
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md),
[grViz](https://rich-iannone.github.io/DiagrammeR/reference/grViz.html)

## Examples

``` r
if (FALSE) { # \dontrun{
model_Bergami_int="
  # Common factor and composite models
  OrgPres <~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8 
  OrgIden =~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
  AffJoy =~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
  AffLove  =~ orgcmt5 + orgcmt6 + orgcmt8

  # Structural model 
  OrgIden ~ OrgPres 
  AffLove ~ OrgPres+OrgIden+OrgPres.OrgIden
  AffJoy  ~ OrgPres+OrgIden
  "
  
  outBergamiInt <- csem(.data = BergamiBagozzi2000,.model = model_Bergami_int,
                        .disattenuate = T,
                        .PLS_weight_scheme_inner = 'factorial',
                        .tolerance = 1e-6,
                        .resample_method = 'none')
  
  outPlot <- plot(outBergamiInt)
  outPlot
  savePlot(outPlot,.file='plot.pdf')
  savePlot(outPlot,.file='plot.png')
  savePlot(outPlot,.file='plot.svg')
  savePlot(outPlot,.file='plot.dot')
} # }
```
