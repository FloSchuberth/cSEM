# NA

## Notes to contributors

Stick to the structure, design choices and style conventions described
below. For questions: please contact the
[author](mailto:f.schuberth@utwente.nl).

### General design choices

- The only OO system used is the S3 system. No S4 classes will be
  allowed!
- Whenever you subset a matrix using `[` use: `[..., ..., drop = FALSE]`
  to avoid accidentally dropping the `dim` attributes.
- Generally avoid using attributes (sometimes they are meaningful
  though)
- Whenever the output consists of more than 1 element use a named list!

### Style/Naming

Stick to [this styleguide](https://style.tidyverse.org/) with the
following exceptions/additions:

1.  Function and class names are always CamelCase. Function names should
    contain a verb followed by a noun like:
    [`processData()`](https://floschuberth.github.io/cSEM/reference/processData.md),
    `calculateValue()`.
2.  Verbs in function names should be consistent across the whole
    package. Avoid mixing synonyms. Example: `computeValue()`
    vs. `calculateValue()`. This package always uses `calculate` instead
    of `compute`. Similarly, `method` vs e.g. `approach`. This package
    always uses `approach` instead of `method`.
3.  Use plural in function/object names if the main output is more than
    one element, like
    [`scaleWeights()`](https://floschuberth.github.io/cSEM/reference/scaleWeights.md),
    `calculateComposites()`,
    [`handleArgs()`](https://floschuberth.github.io/cSEM/reference/handleArgs.md)
    etc. but stick to singular in other cases like
    [`parseModel()`](https://floschuberth.github.io/cSEM/reference/parseModel.md).
4.  Strive for meaningful argument names even if they are a bit longer
    than usual. People are much better at remembering arguments like
    `respect_structural_model` compared to something like `resp_sm`.
    Naming should also be consistent if possible. For example: any
    argument that describes a method or an approach should be named
    `.approach_*`.
5.  Argument names always start with a dot to distinguish them from
    other objects.
6.  Indentation: It is OK to align function arguments indented by two
    spaces below the function name instead of where the function starts
    if this help with readability.

``` r
## Both ok but second is prefered in this case
calculateInnerWeightsPLS <- function(.S                           = NULL,
                                     .W                           = NULL,
                                     .csem_model                  = NULL,
                                     .PLS_weight_scheme_inner     = c("centroid",
                                                                      "factorial", 
                                                                      "path"),
                                     .PLS_ignore_structural_model = NULL
                                     ) { }

calculateInnerWeightsPLS <- function(
  .S                            = NULL,
  .W                            = NULL,
  .csem_model                   = NULL,
  .PLS_weight_scheme_inner      = c("centroid","factorial", "path"),
  .PLS_ignore_structural_model  = NULL
  ) { }
```

### Matrices

- Whenever possible: variables belong in columns, observations belong to
  rows. This implies: whenever a matrix contains observations, variables
  are in columns, no matter if they are indicators, proxies, errors or
  anything else.
- Covariance matrices: indicators **always** belong to columns and
  proxies to rows, i.e., the matrix of weights $W$ for PLS is therefore
  $(J \times K)$ where $J$ is the number of proxies and $K$ the number
  of indicators.
- Naming: matrices within the package should be named according to the
  naming schemes of the related SEM literature.

### Arguments

All arguments used by any of the functions in the package (including
internal functions) are centrally collected in the file
`zz_arguments.R`. Whenever a new argument is introduced:

1.  Add the new argument name + a description to the `cSEMArguments`
    list in alphabetical order by writing
    `@param <argument> "Some description". "Defaults to xxx".`
2.  Add the argument to the `args` or the `args_dotdotdot_csem` list of
    the
    [`args_default()`](https://floschuberth.github.io/cSEM/reference/args_default.md)
    function and provide a default value.
    - Arguments for `args_dotdotdot_csem`: all arguments that can be
      used when calling
      [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
      or `cca()` . Practically this comprises all formal arguments of
      [`foreman()`](https://floschuberth.github.io/cSEM/reference/foreman.md)
      that are not formal arguments of
      [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).
    - Arguments for `args`: all other arguments.
3.  Add the argument to the function you want to use it in. Order
    arguments according to their importance (i.e. `.data` and `.model`
    always come first). if there is one otherwise use alphabetical
    order.
