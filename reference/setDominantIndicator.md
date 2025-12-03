# Internal: Set the dominant indicator

Set the dominant indicator for each construct. Since the sign of the
weights, and thus the loadings is often not determined, a dominant
indicator can be chosen per block. The sign of the weights are chosen
that the correlation between the dominant indicator and the composite is
positive.

## Usage

``` r
setDominantIndicator(
 .W                   = args_default()$.W,
 .dominant_indicators = args_default()$.dominant_indicators, 
 .S                   = args_default()$.S
 )
```

## Arguments

- .W:

  A (J x K) matrix of weights.

- .dominant_indicators:

  A character vector of `"construct_name" = "indicator_name"` pairs,
  where `"indicator_name"` is a character string giving the name of the
  dominant indicator and `"construct_name"` a character string of the
  corresponding construct name. Dominant indicators may be specified for
  a subset of the constructs. Default to `NULL`.

- .S:

  The (K x K) empirical indicator correlation matrix.

## Value

The (J x K) matrix of weights with the dominant indicator set.
