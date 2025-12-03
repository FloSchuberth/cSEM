# Internal: Calculate composite variance-covariance matrix

Calculate the sample variance-covariance (VCV) matrix of the
composites/proxies.

## Usage

``` r
calculateCompositeVCV(
 .S  = args_default()$.S,
 .W  = args_default()$.W
 )
```

## Arguments

- .S:

  The (K x K) empirical indicator correlation matrix.

- .W:

  A (J x K) matrix of weights.

## Value

A (J x J) composite VCV matrix.
