# Calculate composite weights using GSCA

Calculate composite weights using generalized structure component
analysis (GSCA). The first version of this approach was presented in
Hwang and Takane (2004) . Since then, several advancements have been
proposed. The latest version of GSCA can been found in Hwang and Takane
(2014) . This is the version cSEMs implementation is based on.

## Usage

``` r
calculateWeightsGSCA(
  .X                           = args_default()$.X,
  .S                           = args_default()$.S,
  .csem_model                  = args_default()$.csem_model,
  .conv_criterion              = args_default()$.conv_criterion,
  .iter_max                    = args_default()$.iter_max,
  .starting_values             = args_default()$.starting_values,
  .tolerance                   = args_default()$.tolerance
   )
```

## Arguments

- .X:

  A matrix of processed data (scaled, cleaned and ordered).

- .S:

  The (K x K) empirical indicator correlation matrix.

- .csem_model:

  A (possibly incomplete)
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)-list.

- .conv_criterion:

  Character string. The criterion to use for the convergence check. One
  of: "*diff_absolute*", "*diff_squared*", or "*diff_relative*".
  Defaults to "*diff_absolute*".

- .iter_max:

  Integer. The maximum number of iterations allowed. If `iter_max = 1`
  and `.approach_weights = "PLS-PM"` one-step weights are returned. If
  the algorithm exceeds the specified number, weights of iteration step
  `.iter_max - 1` will be returned with a warning. Defaults to `100`.

- .starting_values:

  A named list of vectors where the list names are the construct names
  whose indicator weights the user wishes to set. The vectors must be
  named vectors of `"indicator_name" = value` pairs, where `value` is
  the (scaled or unscaled) starting weight. Defaults to `NULL`.

- .tolerance:

  Double. The tolerance criterion for convergence. Defaults to `1e-05`.

## Value

A named list. J stands for the number of constructs and K for the number
of indicators.

- `$W`:

  A (J x K) matrix of estimated weights.

- `$E`:

  `NULL`

- `$Modes`:

  A named vector of Modes used for the outer estimation, for GSCA the
  mode is automatically set to "gsca".

- `$Conv_status`:

  The convergence status. `TRUE` if the algorithm has converged and
  `FALSE` otherwise.

- `$Iterations`:

  The number of iterations required.

## References

Hwang H, Takane Y (2004). “Generalized Structured Component Analysis.”
*Psychometrika*, **69**(1), 81–99.  
  
Hwang H, Takane Y (2014). *Generalized Structured Component Analysis: A
Component-Based Approach to Structural Equation Modeling*, Chapman &
Hall/CRC Statistics in the Social and Behavioral Sciences. Chapman and
Hall/CRC.
