# `cSEMNonlinearEffects` method for `plot()`

This plot method can be used to create plots to analyze non-linear
models in more depth. In doing so the following plot types can be
selected:

- `.plot_type = "simpleeffects"`::

  The plot of a simple effects analysis displays the predicted value of
  the dependent variable for different values of the independent
  variable and the moderator. As levels for the moderator the levels
  provided to the
  [`doNonlinearEffectsAnalysis()`](https://floschuberth.github.io/cSEM/reference/doNonlinearEffectsAnalysis.md)
  function are used. Since the constructs are standardized the values of
  the moderator equals the deviation from its mean measured in standard
  deviations.

- `.plot_type = "surface"`::

  The plot of a surface analysis displays the predicted values of an
  independent variable (z). The values are predicted based on the values
  of the moderator and the independent variable including all their
  higher-order terms. For the values of the moderator and the
  independent variable steps between their minimum and maximum values
  are used.

- `.plot_type = "floodlight"`::

  The plot of a floodlight analysis displays the direct effect of an
  continuous independent variable (z) on a dependent variable (y)
  conditional on the values of a continuous moderator variable (x),
  including the confidence interval and the Johnson-Neyman points. It is
  noted that in the floodlight plot only moderation is taken into
  account and higher order terms are ignored. For more details, see
  Spiller et al. (2013) .

Plot the predicted values of an independent variable (z) The values are
predicted based on a certain moderator and a certain independent
variable including all their higher-order terms.

## Usage

``` r
# S3 method for class 'cSEMNonlinearEffects'
plot(x, .plot_type = "simpleeffects", .plot_package = "plotly", ...)
```

## Arguments

- x:

  An R object of class `cSEMNonlinearEffects`.

- .plot_type:

  A character string indicating the type of plot that should be
  produced. Options are "*simpleeffects*", "*surface*", and
  "*floodlight*". Defaults to "*simpleeffects*".

- .plot_package:

  A character vector indicating the plot package used. Options are
  "*plotly*", and "*persp*". Defaults to "*plotly*".

- ...:

  Additional parameters that can be passed to
  [`graphics::persp`](https://rdrr.io/r/graphics/persp.html), e.g., to
  rotate the plot.

## See also

[`doNonlinearEffectsAnalysis()`](https://floschuberth.github.io/cSEM/reference/doNonlinearEffectsAnalysis.md)
