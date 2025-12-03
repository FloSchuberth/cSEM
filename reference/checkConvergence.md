# Internal: Check convergence

Check convergence of an algorithm using one of the following criteria:

- `diff_absolute`:

  Checks if the largest elementwise absolute difference between two
  matrices `.W_new` and `W.old` is smaller than a given tolerance.

- `diff_squared`:

  Checks if the largest elementwise squared difference between two
  matrices `.W_new` and `W.old` is smaller than a given tolerance.

- `diff_relative`:

  Checks if the largest elementwise absolute rate of change (new - old /
  new) for two matrices `.W_new` and `W.old` is smaller than a given
  tolerance.

## Usage

``` r
checkConvergence(
  .W_new          = args_default()$.W_new,
  .W_old          = args_default()$.W_old,
  .conv_criterion = args_default()$.conv_criterion,
  .tolerance      = args_default()$.tolerance
  )
```

## Arguments

- .W_new:

  A (J x K) matrix of weights.

- .W_old:

  A (J x K) matrix of weights.

- .conv_criterion:

  Character string. The criterion to use for the convergence check. One
  of: "*diff_absolute*", "*diff_squared*", or "*diff_relative*".
  Defaults to "*diff_absolute*".

- .tolerance:

  Double. The tolerance criterion for convergence. Defaults to `1e-05`.

## Value

`TRUE` if converged; `FALSE` otherwise.
