# Internal: Parameter differences across groups

Calculate the difference between one or more parameter estimates across
all possible pairs of groups (data sets) in `.object`.

## Usage

``` r
calculateParameterDifference(
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

A list of length equal to the number of possible pairs of groups in
`.object` (mathematically, this is n choose 2, i.e., 3 if there are
three groups and 6 if there are 4 groups). Each list elements is itself
a list of three. The first list element contains the difference between
parameter estimates of the structural model, the second list element the
difference between estimated loadings, and the third the difference
between estimated weights.
