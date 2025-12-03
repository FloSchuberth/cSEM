# Internal: Convert second order cSEMModel

Uses a
[cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)
containing second order constructs and turns it into an estimable model
using either the "2stage" approach or the "mixed" approach.

## Usage

``` r
convertModel(
 .csem_model        = NULL, 
 .approach_2ndorder = "2stage",
 .stage             = "first"
 )
```

## Arguments

- .csem_model:

  A (possibly incomplete)
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)-list.

- .approach_2ndorder:

  Character string. Approach used for models containing second-order
  constructs. One of: "*2stage*", or "*mixed*". Defaults to "*2stage*".

- .stage:

  Character string. The stage the model is needed for. One of "*first*"
  or "*second*". Defaults to "*first*".

## Value

A
[cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)
list that may be passed to any function requiring `.csem_model` as a
mandatory argument.
