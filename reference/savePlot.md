# savePlot

This function saves a given plot of a cSEMResults object to a specified
file format.

## Usage

``` r
savePlot(
 .plot_object,
 .filename,
 .path = NULL)
```

## Arguments

- .plot_object:

  Object returned by one of the following functions
  [`plot.cSEMResults_default()`](https://floschuberth.github.io/cSEM/reference/plot.cSEMResults_default.md),
  [`plot.cSEMResults_multi()`](https://floschuberth.github.io/cSEM/reference/plot.cSEMResults_multi.md),
  or
  [`plot.cSEMResults_2ndorder()`](https://floschuberth.github.io/cSEM/reference/plot.cSEMResults_2ndorder.md).

- .filename:

  Character string. The name of the file to save the plot to (supports
  'pdf', 'png', 'svg', and 'dot' formats).

- .path:

  Character string. Path of the directory to save the file to. Defaults
  to the current working directory.

## See also

[`plot.cSEMResults_default()`](https://floschuberth.github.io/cSEM/reference/plot.cSEMResults_default.md)
[`plot.cSEMResults_multi()`](https://floschuberth.github.io/cSEM/reference/plot.cSEMResults_multi.md)
[`plot.cSEMResults_2ndorder()`](https://floschuberth.github.io/cSEM/reference/plot.cSEMResults_2ndorder.md)
