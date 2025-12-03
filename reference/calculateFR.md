# Internal: ANOVA F-test statistic

Calculate the ANOVA F-test statistic suggested by Sarstedt et al. (2011)
in the OTG testing procedure.

## Usage

``` r
calculateFR(.resample_sarstedt)
```

## Arguments

- .resample_sarstedt:

  A matrix containing the parameter estimates that could potentially be
  compared and an id column indicating the group adherence of each row.

## Value

A named scalar, the test statistic of the ANOVA F-test

## References

Sarstedt M, Henseler J, Ringle CM (2011). “Multigroup Analysis in
Partial Least Squares (PLS) Path Modeling: Alternative Methods and
Empirical Results.” In *Advances in International Marketing*, 195–218.
Emerald Group Publishing Limited.
[doi:10.1108/s1474-7979(2011)0000022012](https://doi.org/10.1108/s1474-7979%282011%290000022012)
.
