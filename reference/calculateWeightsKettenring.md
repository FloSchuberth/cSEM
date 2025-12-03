# Calculate composite weights using GCCA

Calculates composite weights according to one of the the five criteria
"*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*", and "*GENVAR*"
suggested by Kettenring (1971) .

## Usage

``` r
calculateWeightsKettenring(
  .S              = args_default()$.S, 
  .csem_model     = args_default()$.csem_model,   
  .approach_gcca  = args_default()$.approach_gcca
  )
```

## Arguments

- .S:

  The (K x K) empirical indicator correlation matrix.

- .csem_model:

  A (possibly incomplete)
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)-list.

- .approach_gcca:

  Character string. The Kettenring approach to use for GCCA. One of
  "*SUMCORR*", "*MAXVAR*", "*SSQCORR*", "*MINVAR*" or "*GENVAR*".
  Defaults to "*SUMCORR*".

## Value

A named list. J stands for the number of constructs and K for the number
of indicators.

- `$W`:

  A (J x K) matrix of estimated weights.

- `$E`:

  `NULL`

- `$Modes`:

  The GCCA mode used for the estimation.

- `$Conv_status`:

  The convergence status. `TRUE` if the algorithm has converged and
  `FALSE` otherwise. For `.approach_gcca = "MINVAR"` or
  `.approach_gcca = "MAXVAR"` the convergence status is `NULL` since
  both are closed-form estimators.

- `$Iterations`:

  The number of iterations required. 0 for `.approach_gcca = "MINVAR"`
  or `.approach_gcca = "MAXVAR"`

## References

Kettenring JR (1971). “Canonical Analysis of Several Sets of Variables.”
*Biometrika*, **58**(3), 433–451.
