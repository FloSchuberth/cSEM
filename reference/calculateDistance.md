# Internal: Matrix difference

Calculates the average of the differences between all possible pairs of
(symmetric) matrices in a list using a given distance measure.

## Usage

``` r
calculateDistance(
  .matrices = NULL, 
  .distance = args_default()$.distance
  )
```

## Arguments

- .matrices:

  A list of at least two matrices.

- .distance:

  Character string. A distance measure. One of: "*geodesic*" or
  "*squared_euclidean*". Defaults to "*geodesic*".

## Value

A numeric vector of length one containing the (arithmetic) mean of the
differences between all possible pairs of matrices supplied via
`.matrices`.

## Details

`.matrices` must be a list of at least two matrices. If more than two
matrices are supplied the arithmetic mean of the differences between all
possible pairs of (symmetric) matrices in a list is computed.
Mathematically this is n chose 2. Hence, supplying a large number of
matrices will become computationally challenging.

Currently two distance measures are supported:

- `geodesic`:

  (Default) The geodesic distance.

- `squared_euclidean`:

  The squared Euclidean distance
