# Internal: Parameter names

Based on a model in [lavaan model
syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html), returns the
names of the parameters of the structural model, the
measurement/composite model and the weight relationship. Used by
[`testMGD()`](https://floschuberth.github.io/cSEM/reference/testMGD.md)
to extract the names of the parameters to compare across groups
according to the test proposed by Chin and Dibbern (2010) .

## Usage

``` r
getParameterNames(
           .object  = args_default()$.object,
           .model   = args_default()$.model
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
  which parameters (i.e, path (`~`), loadings (`=~`), or weights (`<~`))
  should be compared across groups. Defaults to `NULL` in which case all
  parameters of the model are compared.

## Value

A list with elements `names_path`, `names_loadings`, and `names_weights`
containing the names of the structural parameters, the loadings, and the
weight to compare across groups.

## References

Chin WW, Dibbern J (2010). “An Introduction to a Permutation Based
Procedure for Multi-Group PLS Analysis: Results of Tests of Differences
on Simulated Data and a Cross Cultural Analysis of the Sourcing of
Information System Services Between Germany and the USA.” In *Handbook
of Partial Least Squares*, 171–193. Springer Berlin Heidelberg.
[doi:10.1007/978-3-540-32827-8_8](https://doi.org/10.1007/978-3-540-32827-8_8)
.
