# Internal: Process data

Prepare, standardize, check, and clean data provided via the `.data`
argument.

## Usage

``` r
processData(
  .data        = NULL, 
  .model       = NULL, 
  .instruments = NULL
  )
```

## Arguments

- .data:

  A `data.frame` or a `matrix` of standardized or unstandardized data
  (indicators/items/manifest variables). Possible column types or
  classes of the data provided are: "`logical`", "`numeric`" ("`double`"
  or "`integer`"), "`factor`" ("`ordered`" and/or "`unordered`"),
  "`character`" (converted to factor), or a mix of several types.

- .model:

  A model in [lavaan model
  syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html) or a
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)
  list.

- .instruments:

  A named list of vectors of instruments. The names of the list elements
  are the names of the dependent (LHS) constructs of the structural
  equation whose explanatory variables are endogenous. The vectors
  contain the names of the instruments corresponding to each equation.
  Note that exogenous variables of a given equation **must** be supplied
  as instruments for themselves. Defaults to `NULL`.

## Value

A (N x K) data.frame containing the standardized data with columns
ordered according to the order they appear in the measurement model
equations provided via the `.model` argument.
