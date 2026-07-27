# Internal: Fast polyserial correlation

Estimate the polyserial correlation between a continuous variable `x`
and an ordinal (categorical) variable `y`, i.e., the correlation between
`x` and the continuous latent variable assumed to underlie `y` (Drasgow
1988) . Implemented by Kjell S. Slupphaug.

## Usage

``` r
polyserial(
  x,
  y,
  control = list(),
  maxrho = 0.999,
  start = rawcor(x, y),
  thresholds = FALSE
)
```

## References

Drasgow F (1988). “Polychoric and polyserial correlations.” In
*Encyclopedia of Statistical Sciences*, volume 7, 68-74. John Wiley &
Sons Inc, Hoboken.
