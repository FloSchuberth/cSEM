# Model-implied indicator or construct variance-covariance matrix

Calculate the model-implied indicator or construct variance-covariance
(VCV) matrix. Currently only the model-implied VCV for recursive linear
models is implemented (including models containing second order
constructs).

## Usage

``` r
fit(
  .object    = NULL, 
  .saturated = args_default()$.saturated,
  .type_vcv  = args_default()$.type_vcv
  )
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .saturated:

  Logical. Should a saturated structural model be used? Defaults to
  `FALSE`.

- .type_vcv:

  Character string. Which model-implied correlation matrix should be
  calculated? One of "*indicator*" or "*construct*". Defaults to
  "*indicator*".

## Value

Either a (K x K) matrix or a (J x J) matrix depending on the `type_vcv`.

## Details

Notation is taken from Bollen (1989) . If `.saturated = TRUE` the
model-implied variance-covariance matrix is calculated for a saturated
structural model (i.e., the VCV of the constructs is replaced by their
correlation matrix). Hence: V(eta) = WSW' (possibly disattenuated).

## References

Bollen KA (1989). *Structural Equations with Latent Variables*.
Wiley-Interscience. ISBN 978-0471011712.

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[`foreman()`](https://floschuberth.github.io/cSEM/reference/foreman.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
