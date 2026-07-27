# Internal: Fast (cross-)tabulation of integer-coded categorical variables

A faster, more specialized alternative to
[`base::table()`](https://rdrr.io/r/base/table.html) for
(cross-)tabulating one or two variables that are (or can be coerced to)
small positive integer codes, as used by
[`polychor()`](https://floschuberth.github.io/cSEM/reference/polychor.md)
and
[`polyserial()`](https://floschuberth.github.io/cSEM/reference/polyserial.md).
Rows with a missing value in either `x` or `y` are dropped before
tabulating.Implemented by Kjell S. Slupphaug.

## Usage

``` r
fastIntTab(x, y = NULL)
```

## References

There are no references for Rd macro `\insertAllCites` on this help
page.
