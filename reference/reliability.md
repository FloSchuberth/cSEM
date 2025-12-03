# Reliability

Compute several reliability estimates. See the
[Reliability](https://floschuberth.github.io/cSEM/articles/Using-assess.html#reliability)
section of the [cSEM
website](https://floschuberth.github.io/cSEM/index.html) for details.

## Usage

``` r
calculateRhoC(
  .object = NULL,
  .model_implied = TRUE,
  .only_common_factors = TRUE,
  .weighted = FALSE
)

calculateRhoT(
  .object = NULL,
  .alpha = 0.05,
  .closed_form_ci = FALSE,
  .only_common_factors = TRUE,
  .output_type = c("vector", "data.frame"),
  .weighted = FALSE,
  ...
)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .model_implied:

  Logical. Should weights be scaled using the model-implied indicator
  correlation matrix? Defaults to `TRUE`.

- .only_common_factors:

  Logical. Should only concepts modeled as common factors be included
  when calculating one of the following quality criteria: AVE, the
  Fornell-Larcker criterion, HTMT, and all reliability estimates.
  Defaults to `TRUE`.

- .weighted:

  Logical. Should estimation be based on a score that uses the weights
  of the weight approach used to obtain `.object`?. Defaults to `FALSE`.

- .alpha:

  An integer or a numeric vector of significance levels. Defaults to
  `0.05`.

- .closed_form_ci:

  Logical. Should a closed-form confidence interval be computed?
  Defaults to `FALSE`.

- .output_type:

  Character string. The type of output. One of "vector" or "data.frame".
  Defaults to "vector".

- ...:

  Ignored.

## Value

For `calculateRhoC()` and `calculateRhoT()` (if
`.output_type = "vector"`) a named numeric vector containing the
reliability estimates. If `.output_type = "data.frame"`
`calculateRhoT()` returns a `data.frame` with as many rows as there are
constructs modeled as common factors in the model (unless
`.only_common_factors = FALSE` in which case the number of rows equals
the total number of constructs in the model). The first column contains
the name of the construct. The second column the reliability estimate.
If `.closed_form_ci = TRUE` the remaining columns contain lower and
upper bounds for the (1 - `.alpha`) confidence interval(s).

## Details

Since reliability is defined with respect to a classical true score
measurement model only concepts modeled as common factors are considered
by default. For concepts modeled as composites reliability may be
estimated by setting `.only_common_factors = FALSE`, however, it is
unclear how to interpret reliability in this case.

Reliability is traditionally based on a test score (proxy) based on unit
weights. To compute congeneric and tau-equivalent reliability based on a
score that uses the weights of the weight approach used to obtain
`.object` use `.weighted = TRUE` instead.

For the tau-equivalent reliability ("`rho_T`" or "`cronbachs_alpha`") a
closed-form confidence interval may be computed (Trinchera et al. 2018)
by setting `.closed_form_ci = TRUE` (default is `FALSE`). If `.alpha` is
a vector several CIs are returned.

## Functions

- `calculateRhoC()`: Calculate the congeneric reliability

- `calculateRhoT()`: Calculate the tau-equivalent reliability

## References

Trinchera L, Marie N, Marcoulides GA (2018). “A Distribution Free
Interval Estimate for Coefficient Alpha.” *Structural Equation Modeling:
A Multidisciplinary Journal*, **25**(6), 876–887.
[doi:10.1080/10705511.2018.1431544](https://doi.org/10.1080/10705511.2018.1431544)
.

## See also

[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
