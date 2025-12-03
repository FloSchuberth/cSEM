# Internal: Calculate construct variance-covariance matrix

Calculate the variance-covariance matrix (VCV) of the constructs, i.e.,
correlations that involve common factors/latent variables are
diattenuated.

## Usage

``` r
calculateConstructVCV(
 .C          = args_default()$.C, 
 .Q          = args_default()$.Q
 )
```

## Arguments

- .C:

  A (J x J) composite variance-covariance matrix.

- .Q:

  A vector of composite-construct correlations with element names equal
  to the names of the J construct names used in the measurement model.
  Note Q^2 is also called the reliability coefficient.

## Value

The (J x J) construct VCV matrix. Disattenuated if requested.
