# Calculate weights using GSCAm

Calculate composite weights using generalized structured component
analysis with uniqueness terms (GSCAm) proposed by Hwang et al. (2017) .

## Usage

``` r
calculateWeightsGSCAm(
  .X                           = args_default()$.X,
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

A list with the elements

- `$W`:

  A (J x K) matrix of estimated weights.

- `$C`:

  The (J x K) matrix of estimated loadings.

- `$B`:

  The (J x J) matrix of estimated path coefficients.

- `$E`:

  `NULL`

- `$Modes`:

  A named vector of Modes used for the outer estimation, for GSCA the
  mode is automatically set to 'gsca'.

- `$Conv_status`:

  The convergence status. `TRUE` if the algorithm has converged and
  `FALSE` otherwise.

- `$Iterations`:

  The number of iterations required.

## Details

If there are only constructs modeled as common factors calling
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md) with
`.appraoch_weights = "GSCA"` will automatically call
`calculateWeightsGSCAm()` unless `.disattenuate = FALSE`. GSCAm
currently only works for pure common factor models. The reason is that
the implementation in cSEM is based on (the appendix) of Hwang et al.
(2017) . Following the appendix, GSCAm fails if there is at least one
construct modeled as a composite because calculating weight estimates
with GSCAm leads to a product involving the measurement matrix. This
matrix does not have full rank if a construct modeled as a composite is
present. The reason is that the measurement matrix has a zero row for
every construct which is a pure composite (i.e. all related loadings are
zero) and, therefore, leads to a non-invertible matrix when multiplying
it with its transposed.

## References

Hwang H, Takane Y, Jung K (2017). “Generalized structured component
analysis with uniqueness terms for accommodating measurement error.”
*Frontiers in Psychology*, **8**(2137), 1–12.
