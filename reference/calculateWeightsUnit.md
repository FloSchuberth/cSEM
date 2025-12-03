# Calculate composite weights using unit weights

Calculate unit weights for all blocks, i.e., each indicator of a block
is equally weighted.

## Usage

``` r
calculateWeightsUnit(
 .S                 = args_default()$.S,
 .csem_model        = args_default()$.csem_model,
 .starting_values   = args_default()$.starting_values
  )
```

## Arguments

- .S:

  The (K x K) empirical indicator correlation matrix.

- .csem_model:

  A (possibly incomplete)
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)-list.

- .starting_values:

  A named list of vectors where the list names are the construct names
  whose indicator weights the user wishes to set. The vectors must be
  named vectors of `"indicator_name" = value` pairs, where `value` is
  the (scaled or unscaled) starting weight. Defaults to `NULL`.

## Value

A named list. J stands for the number of constructs and K for the number
of indicators.

- `$W`:

  A (J x K) matrix of estimated weights.

- `$E`:

  `NULL`

- `$Modes`:

  The mode used. Always "unit".

- `$Conv_status`:

  `NULL` as there are no iterations

- `$Iterations`:

  0 as there are no iterations
