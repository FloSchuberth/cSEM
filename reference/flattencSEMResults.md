# Internal: Flatten a cSEMResults object

Recursively traverses a
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object (or any nested named list / S4 object) and returns a flat, named
list whose names are the full access paths to each leaf element. This is
primarily a helper for regression testing: comparing two flattened
objects element-by-element makes it possible to trace exactly *which*
result has changed between two versions of cSEM instead of only learning
*that* the objects differ.

## Arguments

- .object:

  An object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  or any (nested) named list or S4 object to be flattened.

- .path:

  Character string. The access path of the current element. Used
  internally for the recursion; users should rely on the default.

## Value

A named `list` with one element per leaf of `.object`. The names give
the full `"$"`-separated access path to each leaf.

## Details

A "leaf" is any `NULL`, atomic vector/matrix, factor, function, or
`formula`. All other list-like or S4 objects are descended into. Path
components are separated by `"$"` and follow the same notation used to
address elements of a
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object, e.g. `"Estimates$Path_estimates"`.
