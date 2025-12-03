# `cSEMTestMGD` method for `print()`

The `cSEMTestMGD` method for the generic function
[`print()`](https://rdrr.io/r/base/print.html).

## Usage

``` r
# S3 method for class 'cSEMTestMGD'
print(
  x,
  .approach_mgd = c("none", "Klesel", "Chin", "Sarstedt", "Keil", "Nitzl", "Henseler",
    "CI_para", "CI_overlap"),
  ...
)
```

## Arguments

- .approach_mgd:

  Character string or a vector of character strings. For which approach
  should details be displayed? One of: "*none*", "*Klesel*", "*Chin*",
  "*Sarstedt*", "*Keil*, "*Nitzl*", "*Henseler*", "*CI_para*", or
  "*CI_overlap*". Default to "*none*" in which case no details are
  displayed.

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md),
[`testMGD()`](https://floschuberth.github.io/cSEM/reference/testMGD.md)
