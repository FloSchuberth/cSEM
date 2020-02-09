
## Notes to contributors

Stick to the structure, design choices and style conventions described
below. For questions: please contact the
[author](mailto:manuel.rademaker@uni-wuerzburg.de).

<!-- ### Structure -->

<!-- The package structure is best understood with reference to a (hierarchically organized)  -->

<!-- company. In analogy to company departments, two separate "function departments"  -->

<!-- exist: -->

<!-- - **Estimation functions** -->

<!-- - **Postestimation functions** -->

<!-- Each department is hierarchically structured. All functions that do not belong -->

<!-- to one of the departments are considered **utility functions**. To stay in the -->

<!-- picture, they are best understood as external consultants in charge of one -->

<!-- specific task that is not meaningfully classified as belonging to (only) one of -->

<!-- the two departments. -->

<!-- #### Estimation functions: -->

<!-- Estimation functions are all functions involved in the estimation/calculation  -->

<!-- of the quantities of the measurement and the path/structural model. There are: -->

<!-- - **2 toplevel functions** (`csem()`). These functions + the -->

<!--   postestimation function `summarize()` (see below) should be sufficient for the end -->

<!--   user most of the time. Both functions eventually call `foreman()`. -->

<!-- - **1 midlevel function**, `foreman()`, that acts like a foreman by collecting all -->

<!--   (estimation) tasks, distributing them to lower level (helper) functions, and  -->

<!--   eventually recollecting all of their results. -->

<!-- - **Lowlevel (helper) functions** that perform one specific task. A distinction -->

<!--   is made between: -->

<!--     - **Exported helper functions**: Functions that are exported to be applied -->

<!--       directly by the end user if needed (e.g. `parseModel()` or -->

<!--       `calculateWeightsPLS()`). -->

<!--     - **Internal helper functions**: Functions that are not exported. These functions  -->

<!--       can be accessed via `csem:::` but they are not generally useful to the end user -->

<!--       (e.g. `calculateConstructVCV()` or `classifyConstructs()`).  -->

<!-- #### Postestimation functions: -->

<!-- Postestimation functions are functions to be applied to an object resulting  -->

<!-- from a call to `csem()`, `cca()` or `foreman()`, namely a `cSEMResults` object.  -->

<!-- All postestimation functions are consistently named by (preferably short) verbs. -->

<!-- Postestimation functions are generic with methods for classes `cSEMResults_default`, -->

<!-- `cSEMResults_multi` and `cSEMResults_2ndorder`. -->

<!-- Currently only a subset of these functions is implemented: -->

<!-- - 1 toplevel function:  -->

<!--   - `summarize()` with class `cSEMSummarize` -->

<!-- - 2 midlevel functions: -->

<!--   - `assess()` with class `cSEMCheck` -->

<!--   - `verify` with class `cSEMVerify` -->

<!-- - 3 lowlevel functions: -->

<!--   - `testMICOM()` with class `cSEMTestMICOM` -->

<!--   - `testMGD()` with class `cSEMTestMGA` -->

<!--   - `testOMF()` with class `cSEMTestOMF` -->

<!-- Each function has (will eventually have) a distinct class and a corresponding  -->

<!-- `print` method. -->

<!-- ### Helper functions -->

<!-- #### Exported helper functions -->

<!-- Exported helper functions should be written as autonomous as possible in a sense  -->

<!-- that they can be used without having to jump to a mother function in order to allow  -->

<!-- researchers using the package to use helper function the way they need it.  -->

<!-- Flexibility will come at the price of code repetition (i.e. most exported helper  -->

<!-- functions will have to have a `parseModel()` + `processDate()` statement at the  -->

<!-- beginning) to make them autonomous. -->

<!-- #### Internal helper functions -->

<!-- Internal helper functions on the other hand do not need to be autonomous. Here -->

<!-- code compactness is preferred.  -->

<!-- Mark every internal function as internal with `@keywords internal`. -->

### General design choices

  - The only OO system used is the S3 system. No S4 classes will be
    allowed\!
  - Whenever you subset a matrix using `[` use: `[..., ..., drop =
    FALSE]` to avoid accidentally dropping the `dim` attributes.
  - Generally avoid using attributs (sometimes they are meaningful
    though)
  - Whenever the output consists of more than 1 element use a named
    list\!

### Style/Naming

Stick to [this styleguide](http://style.tidyverse.org/) with the
following exceptions/additions:

1.  Function and class names are always CamelCase. Function names should
    contain a verb followed by a noun like: `processData()`,
    `calculateValue()`.
2.  Verbs in function names should be consistent across the whole
    package. Avoid to mix synonyms. Example: `computeValue()`
    vs. `calculateValue()`. This package always uses `calculate`
    instead of `compute`. Similarly, `method` vs e.g. `approach`. This
    package always uses `approach` instead of `method`.
3.  Use plural in function/object names if the main output is more than
    one element, like `scaleWeights()`, `calculateComposites()`,
    `handleArgs()` etc. but stick to singular in other cases like
    `parseModel()`.
4.  Strive for meaningful argument names even if they are a bit longer
    than usual. People are much better at remembering arguments like
    `respect_structural_model` compared to something like `resp_sm`.
    Naming should also be consistent if possible. For example: any
    argument that describes a method or an approach should be named
    `.appraoch_*`.
5.  Argument names always start with a dot to distinguish them from
    other objects.
6.  Indentation: It is OK to align function arguments indented by two
    spaces below the function name instead of where the function starts
    if this help with clarity.

<!-- end list -->

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

  - Whenever possible: variables belong in columns, observations belong
    to rows. This implies: whenever a matrix contains observations,
    variables are in columns, no matter if they are indicators, proxies,
    errors or anything else.
  - Covariance matrices: indicators **always** belong to columns and
    proxies to rows, i.e., the matrix of weights \(W\) for PLS is
    therefore \((J \times K)\) where \(J\) is the number of proxies and
    \(K\) the number of indicators.
  - Naming: matrices within the package should be named according to the
    naming schemes of the related SEM literature.

### Arguments

All arguments used by any of the functions in the package (including
internal functions) are centrally collected in the file
`zz_arguments.R`. Whenever a new argument is introduced:

1.  Add the new argument name + a description to the `cSEMArguments`
    list in alphabetical order by writing `@param <argument> "Some
    description". "Defaults to xxx".`
2.  Add the argument to the `args` or the `args_dotdotdot_csem` list of
    the `args_default()` function and provide a default value.
      - Arguments for `args_dotdotdot_csem`: all arguments that can be
        used when calling `csem()` or `cca()` . Practically this
        comprises all formal arguments of `foreman()` that are not
        formal arguments of `csem()`.
      - Arguments for `args`: all other arguments.
3.  Add the argument to the function you want to use it in. Order
    arguments according to their importance (i.e. `.data` and `.model`
    always come first). if there is one otherwise use alphabetical
    order.
