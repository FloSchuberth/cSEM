# Internal: Calculate the outer weights for PLS-PM

Calculates outer weights in PLS-PM. Currently, the originally suggested
mode A and mode B are suggested. Additionally, non-negative least
squares (modeBNNLS) and weights of principal component analysis (PCA)
are implemented.

## Usage

``` r
calculateOuterWeightsPLS(
   .data   = args_default()$.data,  
   .S      = args_default()$.S,
   .W      = args_default()$.W,
   .E      = args_default()$.E,
   .modes  = args_default()$.modes
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

- .W:

  A (J x K) matrix of weights.

- .E:

  A (J x J) matrix of inner weights.

- .modes:

  A vector giving the mode for each construct in the form
  `"name" = "mode"`. Only used internally.

## Value

A (J x K) matrix of outer weights.
