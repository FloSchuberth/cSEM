# Calculate composite weights using principal component analysis (PCA)

Calculate weights for each block by extracting the first principal
component of the indicator correlation matrix S_jj for each blocks,
i.e., weights are the simply the first eigenvector of S_jj.

## Usage

``` r
calculateWeightsPCA(
 .S                 = args_default()$.S,
 .csem_model        = args_default()$.csem_model
  )
```

## Arguments

- .S:

  The (K x K) empirical indicator correlation matrix.

- .csem_model:

  A (possibly incomplete)
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)-list.

## Value

A named list. J stands for the number of constructs and K for the number
of indicators.

- `$W`:

  A (J x K) matrix of estimated weights.

- `$E`:

  `NULL`

- `$Modes`:

  The mode used. Always "PCA".

- `$Conv_status`:

  `NULL` as there are no iterations

- `$Iterations`:

  0 as there are no iterations
