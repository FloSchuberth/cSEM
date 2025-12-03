# Internal: Second/Third stage of the two-stage approach for second order constructs

Performs the second and third stage for a model containing second order
constructs.

## Usage

``` r
calculate2ndStage(
 .csem_model          = args_default()$.csem_model,
 .first_stage_results = args_default()$.first_stage_results,
 .original_arguments  = args_default()$.original_arguments,
 .approach_2ndorder   = args_default()$.approach_2ndorder
  )
```

## Arguments

- .csem_model:

  A (possibly incomplete)
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)-list.

- .original_arguments:

  The list of arguments used within
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .approach_2ndorder:

  Character string. Approach used for models containing second-order
  constructs. One of: "*2stage*", or "*mixed*". Defaults to "*2stage*".

## Value

A cSEMResults object.
