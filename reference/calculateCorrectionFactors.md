# Internal: Calculate PLSc correction factors

Calculates the correction factor used by PLSc.

## Usage

``` r
calculateCorrectionFactors(
 .S               = args_default()$.S,
 .W               = args_default()$.W,
 .modes           = args_default()$.modes,
 .csem_model      = args_default()$.csem_model,
 .PLS_approach_cf = args_default()$.PLS_approach_cf
 )
```

## Arguments

- .S:

  The (K x K) empirical indicator correlation matrix.

- .W:

  A (J x K) matrix of weights.

- .modes:

  A vector giving the mode for each construct in the form
  `"name" = "mode"`. Only used internally.

- .csem_model:

  A (possibly incomplete)
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)-list.

- .PLS_approach_cf:

  Character string. Approach used to obtain the correction factors for
  PLSc. One of: "*dist_squared_euclid*", "*dist_euclid_weighted*",
  "*fisher_transformed*", "*mean_arithmetic*", "*mean_geometric*",
  "*mean_harmonic*", "*geo_of_harmonic*". Defaults to
  "*dist_squared_euclid*". Ignored if `.disattenuate = FALSE` or if
  `.approach_weights` is not PLS-PM.

## Value

A numeric vector of correction factors with element names equal to the
names of the J constructs used in the measurement model.

## Details

Currently, seven approaches are available:

- "dist_squared_euclid" (default)

- "dist_euclid_weighted"

- "fisher_transformed"

- "mean_geometric"

- "mean_harmonic"

- "mean_arithmetic"

- "geo_of_harmonic" (not yet implemented)

See (Dijkstra 2013) for details.

## References

Dijkstra TK (2013). “A Note on How to Make Partial Least Squares
Consistent.” *Working Paper*.
[doi:10.13140/RG.2.1.4547.5688](https://doi.org/10.13140/RG.2.1.4547.5688)
.
