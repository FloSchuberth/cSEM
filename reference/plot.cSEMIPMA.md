# `cSEMIPMA` method for `plot()`

Plot the importance-performance matrix.

## Usage

``` r
# S3 method for class 'cSEMIPMA'
plot(
  x = NULL,
  .dependent = NULL,
  .attributes = NULL,
  .level = c("construct", "indicator"),
  ...
)
```

## Arguments

- x:

  An R object of class `cSEMIPMA`.

- .dependent:

  Character string. Name of the target construct for which the
  importance-performance matrix should be created.

- .attributes:

  Character string. A vector containing indicator/construct names that
  should be plotted in the importance-performance matrix. It must be at
  least of length 2.

- .level:

  Character string. Indicates the level for which the
  importance-performance matrix should be plotted. One of `"construct"`
  or `"indicator"`. Defaults to `"construct"`.

- ...:

  Currently ignored.

## See also

[`doIPMA()`](https://floschuberth.github.io/cSEM/reference/doIPMA.md)
