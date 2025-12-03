# Parse lavaan model

Turns a model written in [lavaan model
syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html) into a
[cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)
list.

## Usage

``` r
parseModel(
  .model        = NULL, 
  .instruments  = NULL, 
  .check_errors = TRUE
  )
```

## Arguments

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

- .check_errors:

  Logical. Should the model to parse be checked for correctness in a
  sense that all necessary components to estimate the model are given?
  Defaults to `TRUE`.

## Value

An object of class
[cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)
is a standardized list containing the following components. J stands for
the number of constructs and K for the number of indicators.

- `$structural`:

  A matrix mimicking the structural relationship between constructs. If
  constructs are only linearly related, `structural` is of dimension (J
  x J) with row- and column names equal to the construct names. If the
  structural model contains nonlinear relationships `structural` is (J x
  (J + J\*)) where J\* is the number of nonlinear terms. Rows are
  ordered such that exogenous constructs are always first, followed by
  constructs that only depend on exogenous constructs and/or previously
  ordered constructs.

- `$measurement`:

  A (J x K) matrix mimicking the measurement/composite relationship
  between constructs and their related indicators. Rows are in the same
  order as the matrix `$structural` with row names equal to the
  construct names. The order of the columns is such that `$measurement`
  forms a block diagonal matrix.

- `$error_cor`:

  A (K x K) matrix mimicking the measurement error correlation
  relationship. The row and column order is identical to the column
  order of `$measurement`.

- `$cor_specified`:

  A matrix indicating the correlation relationships between any
  variables of the model as specified by the user. Mainly for internal
  purposes. Note that `$cor_specified` may also contain inadmissible
  correlations such as a correlation between measurement errors
  indicators and constructs.

- `$construct_type`:

  A named vector containing the names of each construct and their
  respective type ("Common factor" or "Composite").

- `$construct_order`:

  A named vector containing the names of each construct and their
  respective order ("First order" or "Second order").

- `$model_type`:

  The type of model ("Linear" or "Nonlinear").

- `$instruments`:

  Only if instruments are supplied: a list of structural equations
  relating endogenous RHS variables to instruments.

- `$indicators`:

  The names of the indicators (i.e., observed variables and/or
  first-order constructs)

- `$cons_exo`:

  The names of the exogenous constructs of the structural model (i.e.,
  variables that do not appear on the LHS of any structural equation)

- `$cons_endo`:

  The names of the endogenous constructs of the structural model (i.e.,
  variables that appear on the LHS of at least one structural equation)

- `$vars_2nd`:

  The names of the constructs modeled as second orders.

- `$vars_attached_to_2nd`:

  The names of the constructs forming or building a second order
  construct.

- `$vars_not_attached_to_2nd`:

  The names of the constructs not forming or building a second order
  construct.

It is possible to supply an incomplete list to `parseModel()`, resulting
in an incomplete
[cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)
list which can be passed to all functions that require `.csem_model` as
a mandatory argument. Currently, only the structural and the measurement
matrix are required. However, specifying an incomplete
[cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)
list may lead to unexpected behavior and errors. Use with care.

## Details

Instruments must be supplied separately as a named list of vectors of
instruments. The names of the list elements are the names of the
dependent constructs of the structural equation whose explanatory
variables are endogenous. The vectors contain the names of the
instruments corresponding to each equation. Note that exogenous
variables of a given equation **must** be supplied as instruments for
themselves.

By default `parseModel()` attempts to check if the model provided is
correct in a sense that all necessary components required to estimate
the model are specified (e.g., a construct of the structural model has
at least 1 item). To prevent checking for errors use
`.check_errors = FALSE`.

## Examples

``` r
# ===========================================================================
# Providing a model in lavaan syntax 
# ===========================================================================
model <- "
# Structural model
y1 ~ y2 + y3

# Measurement model
y1 =~ x1 + x2 + x3
y2 =~ x4 + x5
y3 =~ x6 + x7

# Error correlation
x1 ~~ x2
"

m <- parseModel(model)
m
#> $structural
#>    y2 y3 y1
#> y2  0  0  0
#> y3  0  0  0
#> y1  1  1  0
#> 
#> $measurement
#>    x4 x5 x6 x7 x1 x2 x3
#> y2  1  1  0  0  0  0  0
#> y3  0  0  1  1  0  0  0
#> y1  0  0  0  0  1  1  1
#> 
#> $error_cor
#>    x4 x5 x6 x7 x1 x2 x3
#> x4  0  0  0  0  0  0  0
#> x5  0  0  0  0  0  0  0
#> x6  0  0  0  0  0  0  0
#> x7  0  0  0  0  0  0  0
#> x1  0  0  0  0  0  1  0
#> x2  0  0  0  0  1  0  0
#> x3  0  0  0  0  0  0  0
#> 
#> $cor_specified
#>    x1 x2
#> x1  0  1
#> x2  1  0
#> 
#> $construct_type
#>              y2              y3              y1 
#> "Common factor" "Common factor" "Common factor" 
#> 
#> $construct_order
#>            y2            y3            y1 
#> "First order" "First order" "First order" 
#> 
#> $model_type
#> [1] "Linear"
#> 
#> $indicators
#> [1] "x4" "x5" "x6" "x7" "x1" "x2" "x3"
#> 
#> $cons_exo
#> [1] "y2" "y3"
#> 
#> $cons_endo
#> [1] "y1"
#> 
#> $vars_2nd
#> character(0)
#> 
#> $vars_attached_to_2nd
#> character(0)
#> 
#> $vars_not_attached_to_2nd
#> [1] "y1" "y2" "y3"
#> 
#> attr(,"class")
#> [1] "cSEMModel"

# ===========================================================================
# Providing a complete model in cSEM format (class cSEMModel)
# ===========================================================================
# If the model is already a cSEMModel object, the model is returned as is:

identical(m, parseModel(m)) # TRUE
#> [1] TRUE

# ===========================================================================
# Providing a list 
# ===========================================================================
# It is possible to provide a list that contains at least the
# elements "structural" and "measurement". This is generally discouraged
# as this may cause unexpected errors.

m_incomplete <- m[c("structural", "measurement", "construct_type")]
parseModel(m_incomplete)
#> $structural
#>    y2 y3 y1
#> y2  0  0  0
#> y3  0  0  0
#> y1  1  1  0
#> 
#> $measurement
#>    x4 x5 x6 x7 x1 x2 x3
#> y2  1  1  0  0  0  0  0
#> y3  0  0  1  1  0  0  0
#> y1  0  0  0  0  1  1  1
#> 
#> $construct_type
#>              y2              y3              y1 
#> "Common factor" "Common factor" "Common factor" 
#> 
#> attr(,"class")
#> [1] "cSEMModel"

# Providing a list containing list names that are not part of a `cSEMModel`
# causes an error:

if (FALSE) { # \dontrun{
m_incomplete[c("name_a", "name_b")] <- c("hello world", "hello universe")
parseModel(m_incomplete)
} # }

# Failing to provide "structural" or "measurement" also causes an error:

if (FALSE) { # \dontrun{
m_incomplete <- m[c("structural", "construct_type")]
parseModel(m_incomplete)
} # }
```
