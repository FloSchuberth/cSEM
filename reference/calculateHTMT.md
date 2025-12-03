# HTMT

Computes either the heterotrait-monotrait ratio of correlations (HTMT)
based on Henseler et al. (2015) or the HTMT2 proposed by Roemer et al.
(2021) . While the HTMT is a consistent estimator for the construct
correlation in case of tau-equivalent measurement models, the HTMT2 is a
consistent estimator for congeneric measurement models. In general, they
are used to assess discriminant validity.

## Usage

``` r
calculateHTMT(
 .object               = NULL,
 .type_htmt            = c('htmt','htmt2'),
 .absolute             = TRUE,
 .alpha                = 0.05,
 .ci                   = c("CI_percentile", "CI_standard_z", "CI_standard_t", 
                           "CI_basic", "CI_bc", "CI_bca", "CI_t_interval"),
 .inference            = FALSE,
 .only_common_factors  = TRUE,
 .R                    = 499,
 .seed                 = NULL,
 ...
)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .type_htmt:

  Character string indicating the type of HTMT that should be
  calculated, i.e., the original HTMT ("*htmt*") or the HTMT2
  ("*htmt2*"). Defaults to "*htmt*"

- .absolute:

  Logical. Should the absolute HTMT values be returned? Defaults to
  `TRUE` .

- .alpha:

  A numeric value giving the significance level. Defaults to `0.05`.

- .ci:

  A character strings naming the type of confidence interval to use to
  compute the 1-alpha% quantile of the bootstrap HTMT values. For
  possible choices see
  [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md).
  Ignored if `.inference = FALSE`. Defaults to "*CI_percentile*".

- .inference:

  Logical. Should critical values be computed? Defaults to `FALSE`.

- .only_common_factors:

  Logical. Should only concepts modeled as common factors be included
  when calculating one of the following quality criteria: AVE, the
  Fornell-Larcker criterion, HTMT, and all reliability estimates.
  Defaults to `TRUE`.

- .R:

  Integer. The number of bootstrap replications. Defaults to `499`.

- .seed:

  Integer or `NULL`. The random seed to use. Defaults to `NULL` in which
  case an arbitrary seed is chosen. Note that the scope of the seed is
  limited to the body of the function it is used in. Hence, the global
  seed will not be altered!

- ...:

  Ignored.

## Value

A named list containing:

- the values of the HTMT/HTMT2, i.e., a matrix with the HTMT/HTMT2
  values at its lower triangular and if `.inference = TRUE` the upper
  triangular contains the upper limit of the 1-2\*.alpha% bootstrap
  confidence interval if the HTMT/HTMT2 is positive and the lower limit
  if the HTMT/HTMT2 is negative.

- the lower and upper limits of the 1-2\*.alpha% bootstrap confidence
  interval if `.inference = TRUE`; otherwise it is `NULL`.

- the number of admissible bootstrap runs, i.e., the number of
  HTMT/HTMT2 values calculated during bootstrap if `.inference = TRUE`;
  otherwise it is `NULL`. Note, the HTMT2 is based on the geometric and
  thus cannot always be calculated.

## Details

Computation of the HTMT/HTMT2 assumes that all intra-block and
inter-block correlations between indicators are either all-positive or
all-negative. A warning is given if this is not the case.

To obtain bootstrap confidence intervals for the HTMT/HTMT2 values, set
`.inference = TRUE`. To choose the type of confidence interval, use
`.ci`. To control the bootstrap process, arguments `.R` and `.seed` are
available. Note, that `.alpha` is multiplied by two because typically
researchers are interested in one-sided bootstrap confidence intervals
for the HTMT/HTMT2.

Since the HTMT and the HTMT2 both assume a reflective measurement model
only concepts modeled as common factors are considered by default. For
concepts modeled as composites the HTMT may be computed by setting
`.only_common_factors = FALSE`, however, it is unclear how to interpret
values in this case.

## References

Henseler J, Ringle CM, Sarstedt M (2015). “A New Criterion for Assessing
Discriminant Validity in Variance-based Structural Equation Modeling.”
*Journal of the Academy of Marketing Science*, **43**(1), 115–135.
[doi:10.1007/s11747-014-0403-8](https://doi.org/10.1007/s11747-014-0403-8)
.  
  
Roemer E, Schuberth F, Henseler J (2021). “HTMT2 – an improved criterion
for assessing discriminant validity in structural equation modeling.”
*Industrial Management & Data Systems*, **121**(12), 2637–2650.

## See also

[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[csem](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
