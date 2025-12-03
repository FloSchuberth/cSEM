# Internal: Composite-based SEM

The central hub of the cSEM package. It acts like a foreman by
collecting all (estimation) tasks, distributing them to lower level
package functions, and eventually recollecting all of their results. It
is called by
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md) to
manage the actual calculations. It may be called directly by the user,
however, in most cases it will likely be more convenient to use
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
instead.

## Usage

``` r
foreman(
  .data                        = args_default()$.data,
  .model                       = args_default()$.model,
  .approach_cor_robust         = args_default()$.approach_cor_robust,
  .approach_nl                 = args_default()$.approach_nl,
  .approach_paths              = args_default()$.approach_paths,
  .approach_weights            = args_default()$.approach_weights,
  .conv_criterion              = args_default()$.conv_criterion,
  .disattenuate                = args_default()$.disattenuate,
  .dominant_indicators         = args_default()$.dominant_indicators,
  .estimate_structural         = args_default()$.estimate_structural,
  .id                          = args_default()$.id,
  .instruments                 = args_default()$.instruments,
  .iter_max                    = args_default()$.iter_max,
  .normality                   = args_default()$.normality,
  .PLS_approach_cf             = args_default()$.PLS_approach_cf,
  .PLS_ignore_structural_model = args_default()$.PLS_ignore_structural_model,
  .PLS_modes                   = args_default()$.PLS_modes,
  .PLS_weight_scheme_inner     = args_default()$.PLS_weight_scheme_inner,
  .reliabilities               = args_default()$.reliabilities,
  .starting_values             = args_default()$.starting_values,
  .tolerance                   = args_default()$.tolerance
  )
```

## Arguments

- .data:

  A `data.frame` or a `matrix` of standardized or unstandardized data
  (indicators/items/manifest variables). Possible column types or
  classes of the data provided are: "`logical`", "`numeric`" ("`double`"
  or "`integer`"), "`factor`" ("`ordered`" and/or "`unordered`"),
  "`character`" (converted to factor), or a mix of several types.

- .model:

  A model in [lavaan model
  syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html) or a
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)
  list.

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

- .approach_nl:

  Character string. Approach used to estimate nonlinear structural
  relationships. One of: "*sequential*" or "*replace*". Defaults to
  "*sequential*".

- .approach_paths:

  Character string. Approach used to estimate the structural
  coefficients. One of: "*OLS*" or "*2SLS*". If "*2SLS*", instruments
  need to be supplied to `.instruments`. Defaults to "*OLS*".

- .approach_weights:

  Character string. Approach used to obtain composite weights. One of:
  "*PLS-PM*", "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*",
  "*GENVAR*", "*GSCA*", "*PCA*", "*unit*", "*bartlett*", or
  "*regression*". Defaults to "*PLS-PM*".

- .conv_criterion:

  Character string. The criterion to use for the convergence check. One
  of: "*diff_absolute*", "*diff_squared*", or "*diff_relative*".
  Defaults to "*diff_absolute*".

- .disattenuate:

  Logical. Should composite/proxy correlations be disattenuated to yield
  consistent loadings and path estimates if at least one of the
  construct is modeled as a common factor? Defaults to `TRUE`.

- .dominant_indicators:

  A character vector of `"construct_name" = "indicator_name"` pairs,
  where `"indicator_name"` is a character string giving the name of the
  dominant indicator and `"construct_name"` a character string of the
  corresponding construct name. Dominant indicators may be specified for
  a subset of the constructs. Default to `NULL`.

- .estimate_structural:

  Logical. Should the structural coefficients be estimated? Defaults to
  `TRUE`.

- .id:

  Character string or integer. A character string giving the name or an
  integer of the position of the column of `.data` whose levels are used
  to split `.data` into groups. Defaults to `NULL`.

- .instruments:

  A named list of vectors of instruments. The names of the list elements
  are the names of the dependent (LHS) constructs of the structural
  equation whose explanatory variables are endogenous. The vectors
  contain the names of the instruments corresponding to each equation.
  Note that exogenous variables of a given equation **must** be supplied
  as instruments for themselves. Defaults to `NULL`.

- .iter_max:

  Integer. The maximum number of iterations allowed. If `iter_max = 1`
  and `.approach_weights = "PLS-PM"` one-step weights are returned. If
  the algorithm exceeds the specified number, weights of iteration step
  `.iter_max - 1` will be returned with a warning. Defaults to `100`.

- .normality:

  Logical. Should joint normality of \\\[\eta\_{1:p}; \zeta;
  \epsilon\]\\ be assumed in the nonlinear model? See (Dijkstra and
  Schermelleh-Engel 2014) for details. Defaults to `FALSE`. Ignored if
  the model is not nonlinear.

- .PLS_approach_cf:

  Character string. Approach used to obtain the correction factors for
  PLSc. One of: "*dist_squared_euclid*", "*dist_euclid_weighted*",
  "*fisher_transformed*", "*mean_arithmetic*", "*mean_geometric*",
  "*mean_harmonic*", "*geo_of_harmonic*". Defaults to
  "*dist_squared_euclid*". Ignored if `.disattenuate = FALSE` or if
  `.approach_weights` is not PLS-PM.

- .PLS_ignore_structural_model:

  Logical. Should the structural model be ignored when calculating the
  inner weights of the PLS-PM algorithm? Defaults to `FALSE`. Ignored if
  `.approach_weights` is not PLS-PM.

- .PLS_modes:

  Either a named list specifying the mode that should be used for each
  construct in the form `"construct_name" = mode`, a single character
  string giving the mode that should be used for all constructs, or
  `NULL`. Possible choices for `mode` are: "*modeA*", "*modeB*",
  "*modeBNNLS*", "*unit*", "*PCA*", a single integer or a vector of
  fixed weights of the same length as there are indicators for the
  construct given by `"construct_name"`. If only a single number is
  provided this is identical to using unit weights, as weights are
  rescaled such that the related composite has unit variance. Defaults
  to `NULL`. If `NULL` the appropriate mode according to the type of
  construct used is chosen. Ignored if `.approach_weight` is not PLS-PM.

- .PLS_weight_scheme_inner:

  Character string. The inner weighting scheme used by PLS-PM. One of:
  "*centroid*", "*factorial*", or "*path*". Defaults to "*path*".
  Ignored if `.approach_weight` is not PLS-PM.

- .reliabilities:

  A character vector of `"name" = value` pairs, where `value` is a
  number between 0 and 1 and `"name"` a character string of the
  corresponding construct name, or `NULL`. Reliabilities may be given
  for a subset of the constructs. Defaults to `NULL` in which case
  reliabilities are estimated by
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).
  Currently, only supported for `.approach_weights = "PLS-PM"`.

- .starting_values:

  A named list of vectors where the list names are the construct names
  whose indicator weights the user wishes to set. The vectors must be
  named vectors of `"indicator_name" = value` pairs, where `value` is
  the (scaled or unscaled) starting weight. Defaults to `NULL`.

- .tolerance:

  Double. The tolerance criterion for convergence. Defaults to `1e-05`.

## See also

[csem](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
