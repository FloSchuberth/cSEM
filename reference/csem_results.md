# cSEMResults

A call to
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
results in an object with at least two class attributes. The first class
attribute is always `cSEMResults` no matter the type of data or model
provided. The second is one of `cSEMResults_default`,
`cSEMResults_multi`, or `cSEMResults_2ndorder` and depends on the
estimated model and/or the type of data provided to the `.model` and
`.data` arguments of
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md). The
third class attribute `cSEMResults_resampled` is only added if
resampling was conducted.

## Details

Depending on the type of data and/or model provided three different
output types exists.

- \_default:

  This will be the structure for the vast majority of applications. If
  the data is a single `matrix` or `data.frame` with no id-column, the
  result is a `list` with elements:

  `$Estimates`

  :   A list containing a list of estimated quantities.

  `$Information`

  :   A list containing a list of additional information.

  The resulting object has classes `cSEMResults` and
  `cSEMResults_default`.

- \_multi:

  If the data provided is a single `matrix` or `data.frame` containing
  an id-column to split the data by `G` group levels or if a list of `G`
  datasets is provided, the resulting object is a list of `G` lists,
  where `G` is equal to the number of groups or the number of datasets
  in the list of datasets provided. Each of the `G` list elements is
  itself a `cSEMResults_default` object. Hence its structure is
  identical to the structure described in `_default`.

  The resulting object has classes `cSEMResults` and
  `cSEMResults_multi`.

- \_2ndorder:

  A special output is generated if the model to estimate contains
  hierarchical constructs **and** the "2stage" or "mixed" approach is
  used to estimate the model. In this case the resulting object is a
  list containing two elements `First_stage` and ` Second_stage`.

  Each list element is itself a `cSEMResults_default` object. Hence its
  structure is identical to the structure described in `_default`.

If `.resample_method = "bootstrap"` or `.resample_method = "jackknife"`,
resamples are attached to each object. For objects of class
`cSEMResults_default` the resamples are attached to
`.object$Estimates$Estimates_resample`. For objects of class
`cSEMResults_multi` the same is done by group. For objects of class
`cSEMResults_2ndorder` the resamples are attached to the
`.object$Second_stage$Information$Resamples`. All objects containing
these elements gain the `cSEMResults_resampled` class.
