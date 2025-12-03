# Internal: Calculation of the CDF used in Henseler et al. (2009)

Calculates the probability that theta^1 is smaller than or equal to
theta^2. See Equation (6) in Sarstedt et al. (2011) .

## Usage

``` r
calculatePr(.resample_centered = NULL, .parameters_to_compare = NULL)
```

## Arguments

- .parameters_to_compare:

  A model in [lavaan model
  syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html) indicating
  which parameters (i.e, path (`~`), loadings (`=~`), weights (`<~`), or
  correlations (`~~`)) should be compared across groups. Defaults to
  `NULL` in which case all weights, loadings and path coefficients of
  the originally specified model are compared.

## Value

A named vector

## References

Sarstedt M, Henseler J, Ringle CM (2011). “Multigroup Analysis in
Partial Least Squares (PLS) Path Modeling: Alternative Methods and
Empirical Results.” In *Advances in International Marketing*, 195–218.
Emerald Group Publishing Limited.
[doi:10.1108/s1474-7979(2011)0000022012](https://doi.org/10.1108/s1474-7979%282011%290000022012)
.
