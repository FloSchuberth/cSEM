# Internal: Calculate indicator correlation matrix

Calculate the indicator correlation matrix using conventional or robust
methods.

## Usage

``` r
calculateIndicatorCor(
  .X_cleaned           = NULL, 
  .approach_cor_robust = "none"
 )
```

## Arguments

- .X_cleaned:

  A data.frame of processed data (cleaned and ordered). Note:
  `X_cleaned` may not be scaled!

- .approach_cor_robust:

  Character string. Approach used to obtain a robust indicator
  correlation matrix. One of: "*none*" in which case the standard
  Bravais-Pearson correlation is used, "*spearman*" for the Spearman
  rank correlation, or "*mcd*" via
  [`MASS::cov.rob()`](https://rdrr.io/pkg/MASS/man/cov.rob.html) for a
  robust correlation matrix. Defaults to "*none*". Note that many
  postestimation procedures (such as
  [`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md)
  or [`fit()`](https://floschuberth.github.io/cSEM/reference/fit.md)
  implicitly assume a continuous indicator correlation matrix (e.g.
  Bravais-Pearson correlation matrix). Only use if you know what you are
  doing.

## Value

A list with elements:

- `$S`:

  The (K x K) indicator correlation matrix

- `$cor_type`:

  The type(s) of indicator correlation computed ( "Pearson",
  "Polyserial", "Polychoric")

- `$thre_est`:

  Currently ignored (NULL)

## Details

If `.approach_cor_robust = "none"` (the default) the type of correlation
computed depends on the types of the columns of `.X_cleaned` (i.e., the
indicators) involved in the computation.

- `Numeric-numeric`:

  If both columns (indicators) involved are numeric, the Bravais-Pearson
  product-moment correlation is computed (via
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html)).

- `Numeric-factor`:

  If any of the columns is a factor variable, the polyserial correlation
  (Drasgow 1988) is computed (via
  [`polycor::polyserial()`](https://rdrr.io/pkg/polycor/man/polyserial.html)).

- `Factor-factor`:

  If both columns are factor variables, the polychoric correlation
  (Drasgow 1988) is computed (via
  [`polycor::polychor()`](https://rdrr.io/pkg/polycor/man/polychor.html)).

Note: logical input is treated as a 0-1 factor variable.

If `"mcd"` (= minimum covariance determinant), the MCD estimator
(Rousseeuw and Driessen 1999) , a robust covariance estimator, is
applied (via
[`MASS::cov.rob()`](https://rdrr.io/pkg/MASS/man/cov.rob.html)).

If `"spearman"`, the Spearman rank correlation is used (via
[`stats::cor()`](https://rdrr.io/r/stats/cor.html)).

## References

Drasgow F (1988). “Polychoric and polyserial correlations.” In
*Encyclopedia of Statistical Sciences*, volume 7, 68-74. John Wiley &
Sons Inc, Hoboken.  
  
Rousseeuw PJ, Driessen KV (1999). “A Fast Algorithm for the Minimum
Covariance Determinant Estimator.” *Technometrics*, **41**(3), 212–223.
[doi:10.1080/00401706.1999.10485670](https://doi.org/10.1080/00401706.1999.10485670)
.
