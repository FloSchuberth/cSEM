# Do an importance-performance matrix analysis

**\[maturing\]**

## Usage

``` r
doIPMA(.object)
```

## Arguments

- .object:

  A `cSEMResults` object.\`

## Value

A list of class `cSEMIPA` with a corresponding method for
[`plot()`](https://rdrr.io/r/graphics/plot.default.html). See:
[`plot.cSEMIPMA()`](https://floschuberth.github.io/cSEM/reference/plot.cSEMIPMA.md).

## Details

Performs an importance-performance matrix analysis (IPMA).

To calculate the performance and importance, the weights of the
indicators are unstandardized using the standard deviation of the
original indicators but normed to have a length of 1. Normed construct
scores are calculated based on the original indicators and the
unstandardized weights.

The importance is calculated as the mean of the original indicators or
the unstandardized construct scores, respectively. The performance is
calculated as the unstandardized total effect if `.level == "construct"`
and as the normed weight times the unstandardized total effect if
`.level == "indicator"`. The literature recommends to use an estimation
as input for \``doIPMA()` that is based on normed indicators, e.g., by
scaling all indicators to 0 to 100, see e.g., Henseler (2021); Ringle
and Sarstedt (2016) .

Note, indicators are not normed internally, as theoretical
maximum/minimum can differ from the empirical maximum/minimum which
would lead to an incorrect normalization.

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md),
[`plot.cSEMIPMA()`](https://floschuberth.github.io/cSEM/reference/plot.cSEMIPMA.md)
