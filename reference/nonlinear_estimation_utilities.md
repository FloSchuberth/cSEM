# Internal: Utility functions for the estimation of nonlinear models

Internal: Utility functions for the estimation of nonlinear models

## Usage

``` r
f1(.i, .j)

f2(.i, .j, .select_from, .Q, .H)

f3(.i, .j, .Q, .H)

f4(.i, .j, .Q, .H, .var_struc_error, .temp = NULL)

f5(.i, .j, .H, .Q, .var_struc_error)
```

## Arguments

- .i:

  Row index

- .j:

  Column index

- .select_from:

  matrix to select from

- .Q:

  A vector of composite-construct correlations with element names equal
  to the names of the J construct names used in the measurement model.
  Note Q^2 is also called the reliability coefficient.

- .H:

  The (N x J) matrix of construct scores.
