# Calculate composite weights using PLS-PM

Calculate composite weights using the partial least squares path
modeling (PLS-PM) algorithm (Wold 1975) .

## Usage

``` r
calculateWeightsPLS(
  .data                        = args_default()$.data,
  .S                           = args_default()$.S,
  .csem_model                  = args_default()$.csem_model,
  .conv_criterion              = args_default()$.conv_criterion,
  .iter_max                    = args_default()$.iter_max,
  .PLS_ignore_structural_model = args_default()$.PLS_ignore_structural_model,
  .PLS_modes                   = args_default()$.PLS_modes,
  .PLS_weight_scheme_inner     = args_default()$.PLS_weight_scheme_inner,
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

- .S:

  The (K x K) empirical indicator correlation matrix.

- .csem_model:

  A (possibly incomplete)
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)-list.

- .conv_criterion:

  Character string. The criterion to use for the convergence check. One
  of: "*diff_absolute*", "*diff_squared*", or "*diff_relative*".
  Defaults to "*diff_absolute*".

- .iter_max:

  Integer. The maximum number of iterations allowed. If `iter_max = 1`
  and `.approach_weights = "PLS-PM"` one-step weights are returned. If
  the algorithm exceeds the specified number, weights of iteration step
  `.iter_max - 1` will be returned with a warning. Defaults to `100`.

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

- .starting_values:

  A named list of vectors where the list names are the construct names
  whose indicator weights the user wishes to set. The vectors must be
  named vectors of `"indicator_name" = value` pairs, where `value` is
  the (scaled or unscaled) starting weight. Defaults to `NULL`.

- .tolerance:

  Double. The tolerance criterion for convergence. Defaults to `1e-05`.

## Value

A named list. J stands for the number of constructs and K for the number
of indicators.

- `$W`:

  A (J x K) matrix of estimated weights.

- `$E`:

  A (J x J) matrix of inner weights.

- `$Modes`:

  A named vector of modes used for the outer estimation.

- `$Conv_status`:

  The convergence status. `TRUE` if the algorithm has converged and
  `FALSE` otherwise. If one-step weights are used via `.iter_max = 1` or
  a non-iterative procedure was used, the convergence status is set to
  `NULL`.

- `$Iterations`:

  The number of iterations required.

## References

Wold H (1975). “Path models with latent variables: The NIPALS approach.”
In Blalock HM, Aganbegian A, Borodkin FM, Boudon R, Capecchi V (eds.),
*Quantitative Sociology*, International Perspectives on Mathematical and
Statistical Modeling, 307–357. Academic Press, New York.
