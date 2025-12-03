# Internal: Helper for doNonlinearEffectsAnalysis()

Function that calculates the values required for the floodlight
analysis, namely 1) partial effect of the independent variable on the
dependent variable for each bootstrap run and for the original
estimation for each step of the moderator 2) alpha/2 and 1-alpha/2
quantile of the bootstrap estimates.

## Usage

``` r
getValuesFloodlight(
  .model = NULL,
  .path_coefficients = args_default()$.path_coefficients,
  .dependent = args_default()$.dependent,
  .independent = args_default()$.independent,
  .moderator = args_default()$.moderator,
  .steps_mod = args_default()$.steps_mod,
  .value_independent = args_default()$.value_independent,
  .alpha = args_default()$.alpha
)
```

## Arguments

- .model:

  A model in [lavaan model
  syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html) or a
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)
  list.

- .path_coefficients:

  List. A list that contains the resampled and the original path
  coefficient estimates. Typically a part of a `cSEMResults_resampled`
  object. Defaults to `NULL`.

- .dependent:

  Character string. The name of the dependent variable.

- .independent:

  Character string. The name of the independent variable.

- .moderator:

  Character string. The name of the moderator variable.

- .steps_mod:

  A numeric vector. Steps used for the moderator variable in calculating
  the simple effects of an independent variable on the dependent
  variable. Defaults to `NULL`.

- .value_independent:

  Integer. Only required for floodlight analysis; The value of the
  independent variable in case that it appears as a higher-order term.

- .alpha:

  An integer or a numeric vector of significance levels. Defaults to
  `0.05`.

## Details

Only variables that comprise the independent variable are taken into
account. If it contains a variable other than the independent variable
and the moderator the effect is set to zero as the other variables are
evaluated at their means (=0), hence the effect is zero.
