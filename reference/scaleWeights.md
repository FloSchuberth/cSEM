# Internal: Scale weights

Scale weights such that the formed composite has unit variance.

## Usage

``` r
scaleWeights(
  .S = args_default()$.S, 
  .W = args_default()$.W
  )
```

## Arguments

- .S:

  The (K x K) empirical indicator correlation matrix.

- .W:

  A (J x K) matrix of weights.

## Value

The (J x K) matrix of scaled weights.
