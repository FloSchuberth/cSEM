# Internal: Calculate Reliabilities

Internal: Calculate Reliabilities

## Usage

``` r
calculateReliabilities(
  .X = args_default()$.X,
  .S = args_default()$.S,
  .W = args_default()$.W,
  .approach_weights = args_default()$.approach_weights,
  .csem_model = args_default()$.csem_model,
  .disattenuate = args_default()$.disattenuate,
  .PLS_approach_cf = args_default()$.PLS_approach_cf,
  .reliabilities = args_default()$.reliabilities
)
```

## Arguments

- .X:

  A matrix of processed data (scaled, cleaned and ordered).

- .S:

  The (K x K) empirical indicator correlation matrix.

- .W:

  A (J x K) matrix of weights.

- .approach_weights:

  Character string. Approach used to obtain composite weights. One of:
  "*PLS-PM*", "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*",
  "*GENVAR*", "*GSCA*", "*PCA*", "*unit*", "*bartlett*", or
  "*regression*". Defaults to "*PLS-PM*".

- .csem_model:

  A (possibly incomplete)
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)-list.

- .disattenuate:

  Logical. Should composite/proxy correlations be disattenuated to yield
  consistent loadings and path estimates if at least one of the
  construct is modeled as a common factor? Defaults to `TRUE`.

- .PLS_approach_cf:

  Character string. Approach used to obtain the correction factors for
  PLSc. One of: "*dist_squared_euclid*", "*dist_euclid_weighted*",
  "*fisher_transformed*", "*mean_arithmetic*", "*mean_geometric*",
  "*mean_harmonic*", "*geo_of_harmonic*". Defaults to
  "*dist_squared_euclid*". Ignored if `.disattenuate = FALSE` or if
  `.approach_weights` is not PLS-PM.

- .reliabilities:

  A character vector of `"name" = value` pairs, where `value` is a
  number between 0 and 1 and `"name"` a character string of the
  corresponding construct name, or `NULL`. Reliabilities may be given
  for a subset of the constructs. Defaults to `NULL` in which case
  reliabilities are estimated by
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).
  Currently, only supported for `.approach_weights = "PLS-PM"`.
