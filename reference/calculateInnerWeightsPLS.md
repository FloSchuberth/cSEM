# Internal: Calculate the inner weights for PLS-PM

PLS-PM forms "inner" composites as a weighted sum of its *I* related
composites. These inner weights are obtained using one of the following
schemes (Lohmöller 1989) :

- `centroid`:

  According to the centroid weighting scheme each inner weight used to
  form composite *j* is either 1 if the correlation between composite
  *j* and its via the structural model related composite *i = 1, ..., I*
  is positive and -1 if it is negative.

- `factorial`:

  According to the factorial weighting scheme each inner weight used to
  form inner composite *j* is equal to the correlation between composite
  *j* and its via the structural model related composite *i = 1, ...,
  I*.

- `path`:

  Lets call all construct that have an arrow pointing to construct *j*
  **predecessors of j** and all arrows going from j to other constructs
  **followers of j**. According the path weighting scheme, inner weights
  are computed as follows. Take construct *j*:

  - For all predecessors of *j* set the inner weight of predecessor *i*
    to the correlation of *i* with *j*.

  - For all followers of *j* set the inner weight of follower *i* to the
    coefficient of a multiple regression of *j* on all followers *i*
    with *i = 1,...,I*.

Except for the path weighting scheme relatedness can come in two
flavors. If `.PLS_ignore_structural_model = TRUE` all constructs are
considered related. If `.PLS_ignore_structural_model = FALSE` (the
default) only adjacent constructs are considered. If
`.PLS_ignore_structural_model = TRUE` and
`.PLS_weight_scheme_inner = "path"` a warning is issued and
`.PLS_ignore_structural_model` is changed to `FALSE`.

## Usage

``` r
calculateInnerWeightsPLS(
  .S                           = args_default()$.S,
  .W                           = args_default()$.W,
  .csem_model                  = args_default()$.csem_model,
  .PLS_ignore_structural_model = args_default()$.PLS_ignore_structrual_model,
  .PLS_weight_scheme_inner     = args_default()$.PLS_weight_scheme_inner
)
```

## Arguments

- .S:

  The (K x K) empirical indicator correlation matrix.

- .W:

  A (J x K) matrix of weights.

- .csem_model:

  A (possibly incomplete)
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)-list.

- .PLS_ignore_structural_model:

  Logical. Should the structural model be ignored when calculating the
  inner weights of the PLS-PM algorithm? Defaults to `FALSE`. Ignored if
  `.approach_weights` is not PLS-PM.

- .PLS_weight_scheme_inner:

  Character string. The inner weighting scheme used by PLS-PM. One of:
  "*centroid*", "*factorial*", or "*path*". Defaults to "*path*".
  Ignored if `.approach_weight` is not PLS-PM.

## Value

The (J x J) matrix `E` of inner weights.
