# Internal: Safer alternative to [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)

A safer alternative to
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) which doesn't
hard error. In the case of an error in
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) (e.g.,. a
NA/NaN gradient) the error is caught, and a list structured like the
return value of [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)
is returned. Used in
[`polychor()`](https://floschuberth.github.io/cSEM/reference/polychor.md)
and
[`polyserial()`](https://floschuberth.github.io/cSEM/reference/polyserial.md).
Implemented by Kjell S. Slupphaug.

## Usage

``` r
snlminb(
  start,
  objective,
  gradient = NULL,
  hessian = NULL,
  ...,
  scale = 1,
  control = list(),
  lower = -Inf,
  upper = Inf,
  warn.on.failure = TRUE
)
```

## References

There are no references for Rd macro `\insertAllCites` on this help
page.
