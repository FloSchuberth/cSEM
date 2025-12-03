# Internal: Estimate the structural coefficients

Estimates the coefficients of the structural model (nonlinear and
linear) using OLS, 2SLS. The latter currently work for linear models
only.

## Usage

``` r
estimatePath(
 .approach_nl      = args_default()$.approach_nl,
 .approach_paths   = args_default()$.approach_paths,
 .csem_model       = args_default()$.csem_model,
 .H                = args_default()$.H,
 .normality        = args_default()$.normality,
 .P                = args_default()$.P,
 .Q                = args_default()$.Q
 )
```

## Arguments

- .approach_nl:

  Character string. Approach used to estimate nonlinear structural
  relationships. One of: "*sequential*" or "*replace*". Defaults to
  "*sequential*".

- .approach_paths:

  Character string. Approach used to estimate the structural
  coefficients. One of: "*OLS*" or "*2SLS*". If "*2SLS*", instruments
  need to be supplied to `.instruments`. Defaults to "*OLS*".

- .csem_model:

  A (possibly incomplete)
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)-list.

- .H:

  The (N x J) matrix of construct scores.

- .normality:

  Logical. Should joint normality of \\\[\eta\_{1:p}; \zeta;
  \epsilon\]\\ be assumed in the nonlinear model? See (Dijkstra and
  Schermelleh-Engel 2014) for details. Defaults to `FALSE`. Ignored if
  the model is not nonlinear.

- .P:

  A (J x J) construct variance-covariance matrix (possibly
  disattenuated).

- .Q:

  A vector of composite-construct correlations with element names equal
  to the names of the J construct names used in the measurement model.
  Note Q^2 is also called the reliability coefficient.

## Value

A named list containing the estimated structural coefficients, the R2,
the adjusted R2, and the VIFs for each regression.
