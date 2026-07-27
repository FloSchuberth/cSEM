# Internal: Fast polychoric correlation

Estimate the polychoric correlation between two ordinal (categorical)
variables `x` and `y`, i.e., the correlation between the two continuous,
bivariate normal latent variables assumed to underlie `x` and `y`
(Drasgow 1988) . Implemented by Kjell S. Slupphaug.

## Usage

``` r
polychor(
  x,
  y,
  control = list(),
  maxrho = 0.999,
  start = NULL,
  thresholds = FALSE
)
```

## References

Drasgow F (1988). “Polychoric and polyserial correlations.” In
*Encyclopedia of Statistical Sciences*, volume 7, 68-74. John Wiley &
Sons Inc, Hoboken.
