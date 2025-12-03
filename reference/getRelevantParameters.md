# Internal: Extract relevant parameters from several cSEMResults_multi

Extract the relevant parameters from a cSEMResult_multi object in
`.object`.

## Usage

``` r
getRelevantParameters(
  .object     = args_default()$.object,
  .model      = args_default()$.model
)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .model:

  A model in [lavaan model
  syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html) indicating
  which parameters (i.e., path (`~`), loadings (`=~`), or weights
  (`<~`)) should be compared across groups. Defaults to `NULL` in which
  case all parameters of the model are compared.

## Value

A list of length equal to the number of groups in `.object`. Each list
element is itself a list of three. The first list element contains the
relevant parameter estimates of the structural model, the second list
element the relevant estimated loadings, and the third the relevant
estimated weights.
