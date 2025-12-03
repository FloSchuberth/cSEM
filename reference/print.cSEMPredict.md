# `cSEMPredict` method for `print()`

The `cSEMPredict` method for the generic function
[`print()`](https://rdrr.io/r/base/print.html).

## Usage

``` r
# S3 method for class 'cSEMPredict'
print(x, .metrics = c("MAE", "RMSE", "Q2"), ...)
```

## Arguments

- .metrics:

  Character string or a vector of character strings. Which prediction
  metrics should be displayed? One of: "*MAE*", "*RMSE*", "*Q2*",
  "*MER*", "*MAPE*, "*MSE2*", "*U1*", "*U2*", "*UM*", "*UR*", or "*UD*".
  Default to c("*MAE*", "*RMSE*", "*Q2*").

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md),
[`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md)
