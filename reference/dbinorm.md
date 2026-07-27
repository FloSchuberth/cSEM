# Internal: Bivariate standard normal density

Evaluate the density of the bivariate standard normal distribution with
correlation `rho` at `(u, v)`. Used internally by
[`polychor()`](https://floschuberth.github.io/cSEM/reference/polychor.md)
to compute the gradient of the polychoric log-likelihood. Implemented by
Kjell S. Slupphaug.

## Usage

``` r
dbinorm(u, v, rho, force.zero = FALSE, rho.lim = 0.9999)
```

## References

There are no references for Rd macro `\insertAllCites` on this help
page.
