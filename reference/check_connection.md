# Internal: Check whether two indicators belong to the same construct.

Checks whether two indicators belong to the same construct.

## Usage

``` r
check_connection(
  .indicator1,
  .indicator2,
  .model_measurement,
  .model_error_cor
)
```

## Arguments

- .indicator1:

  Character string. The name of the indicator 1.

- .indicator2:

  Character string. The name of the indicator 1.

- .model_measurement:

  Matrix. The measurement matrix indicating the relationship between
  constructs and indicators.

- .model_error_cor:

  Matrix. The matrix indicates the error correlation structure.

## Value

TRUE if both indicators belong to the same construct, FALSE otherwise.
