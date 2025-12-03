# Complete list of assess()'s ... arguments

A complete alphabetical list of all possible arguments accepted by
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)'s
`...` (dotdotdot) argument.

## Arguments

- .absolute:

  Logical. Should the absolute HTMT values be returned? Defaults to
  `TRUE` .

- .alpha:

  An integer or a numeric vector of significance levels. Defaults to
  `0.05`.

- .ci:

  A vector of character strings naming the confidence interval to
  compute. For possible choices see
  [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md).

- .closed_form_ci:

  Logical. Should a closed-form confidence interval be computed?
  Defaults to `FALSE`.

- .handle_inadmissibles:

  Character string. How should inadmissible results be treated? One of
  "*drop*", "*ignore*", or "*replace*". If "*drop*", all
  replications/resamples yielding an inadmissible result will be dropped
  (i.e. the number of results returned will potentially be less than
  `.R`). For "*ignore*" all results are returned even if all or some of
  the replications yielded inadmissible results (i.e. number of results
  returned is equal to `.R`). For "*replace*" resampling continues until
  there are exactly `.R` admissible solutions. Depending on the
  frequency of inadmissible solutions this may significantly increase
  computing time. Defaults to "*drop*".

- .inference:

  Logical. Should critical values be computed? Defaults to `FALSE`.

- .null_model:

  Logical. Should the degrees of freedom for the null model be computed?
  Defaults to `FALSE`.

- .R:

  Integer. The number of bootstrap replications. Defaults to `499`.

- .saturated:

  Logical. Should a saturated structural model be used? Defaults to
  `FALSE`.

- .seed:

  Integer or `NULL`. The random seed to use. Defaults to `NULL` in which
  case an arbitrary seed is chosen. Note that the scope of the seed is
  limited to the body of the function it is used in. Hence, the global
  seed will not be altered!

- .type_gfi:

  Character string. Which fitting function should the GFI be based on?
  One of *"ML"* for the maximum likelihood fitting function, *"GLS"* for
  the generalized least squares fitting function or *"ULS"* for the
  unweighted least squares fitting function (same as the squared
  Euclidean distance). Defaults to *"ML"*.

- .type_vcv:

  Character string. Which model-implied correlation matrix should be
  calculated? One of "*indicator*" or "*construct*". Defaults to
  "*indicator*".

## Details

Most arguments supplied to the `...` argument of
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
are only accepted by a subset of the functions called by
[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md).
The following list shows which argument is passed to which function:

- .absolute:

  Accepted by/Passed down to:
  [`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md)

- .alpha:

  Accepted by/Passed down to:
  [`calculateRhoT()`](https://floschuberth.github.io/cSEM/reference/reliability.md),
  [`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md),
  [`calculateCN()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md)

- .ci:

  Accepted by/Passed down to:
  [`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md)

- .closed_form_ci:

  Accepted by/Passed down to:
  [`calculateRhoT()`](https://floschuberth.github.io/cSEM/reference/reliability.md)

- .handle_inadmissibles:

  Accepted by/Passed down to:
  [`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md)

- .inference:

  Accepted by/Passed down to:
  [calculateHTMT](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md)

- .null_model:

  Accepted by/Passed down to:
  [`calculateDf()`](https://floschuberth.github.io/cSEM/reference/calculateDf.md)

- .R:

  Accepted by/Passed down to:
  [`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md)

- .saturated:

  Accepted by/Passed down to:
  [`calculateSRMR()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
  [`calculateDG()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md),
  [`calculateDL()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md),
  [`calculateDML()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md)and
  subsequently
  [`fit()`](https://floschuberth.github.io/cSEM/reference/fit.md).

- .seed:

  Accepted by/Passed down to:
  [`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md)

- .type_gfi:

  Accepted by/Passed down to:
  [`calculateGFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md)

- .type_vcv:

  Accepted by/Passed down to:
  [`calculateSRMR()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md),
  [`calculateDG()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md),
  [`calculateDL()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md),
  [`calculateDML()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md)
  and subsequently
  [`fit()`](https://floschuberth.github.io/cSEM/reference/fit.md).
