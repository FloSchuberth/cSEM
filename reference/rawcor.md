# Internal: Raw starting value for [`polychor()`](https://floschuberth.github.io/cSEM/reference/polychor.md)/[`polyserial()`](https://floschuberth.github.io/cSEM/reference/polyserial.md)

Compute the ordinary Bravais-Pearson correlation between `x` and `y`
after coercing both to numeric. Used as the default starting value for
the `rho` optimization in
[`polychor()`](https://floschuberth.github.io/cSEM/reference/polychor.md)
and
[`polyserial()`](https://floschuberth.github.io/cSEM/reference/polyserial.md);
for
[`polychor()`](https://floschuberth.github.io/cSEM/reference/polychor.md),
`x` and `y` are integer category codes rather than the underlying latent
variables, so this is only a rough approximation of the polychoric
correlation. Implemented by Kjell S. Slupphaug.

## Usage

``` r
rawcor(x, y)
```

## References

There are no references for Rd macro `\insertAllCites` on this help
page.
