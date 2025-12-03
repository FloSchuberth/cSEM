# Internal: Helper for infer()

Collection of various functions that compute an inferential quantity.

## Usage

``` r
MeanResample(.first_resample)

SdResample(.first_resample, .resample_method, .n)

BiasResample(.first_resample, .resample_method, .n)

StandardCIResample(
  .first_resample,
  .bias_corrected,
  .dist = c("z", "t"),
  .df = c("type1", "type2"),
  .resample_method,
  .n,
  .probs
)

PercentilCIResample(.first_resample, .probs)

BasicCIResample(.first_resample, .bias_corrected, .probs)

TStatCIResample(
  .first_resample,
  .second_resample,
  .bias_corrected,
  .resample_method,
  .resample_method2,
  .n,
  .probs
)

BcCIResample(.first_resample, .probs)

BcaCIResample(.object, .first_resample, .probs)
```

## Arguments

- .first_resample:

  A list containing the `.R` resamples based on the original data
  obtained by resamplecSEMResults().

- .resample_method:

  Character string. The resampling method to use. One of: "*none*",
  "*bootstrap*" or "*jackknife*". Defaults to "*none*".

- .n:

  Integer. The number of observations of the original data.

- .bias_corrected:

  Logical. Should the standard and the tStat confidence interval be
  bias-corrected using the bootstrapped bias estimate? If `TRUE` the
  confidence interval for some estimated parameter `theta` is centered
  at `2*theta - theta*_hat`, where `theta*_hat` is the average over all
  `.R` bootstrap estimates of `theta`. Defaults to `TRUE`

- .dist:

  Character string. The distribution to use for the critical value. One
  of *"t"* for Student's t-distribution or *"z"* for the standard normal
  distribution. Defaults to *"z"*.

- .df:

  Character string. The method for obtaining the degrees of freedom.
  Choices are "*type1*" and "*type2*". Defaults to "*type1*" .

- .probs:

  A vector of probabilities.

- .second_resample:

  A list containing `.R2` resamples for each of the `.R` resamples of
  the first run.

- .resample_method2:

  Character string. The resampling method to use when resampling from a
  resample. One of: "*none*", "*bootstrap*" or "*jackknife*". For
  "*bootstrap*" the number of draws is provided via `.R2`. Currently,
  resampling from each resample is only required for the studentized
  confidence interval ("*CI_t_interval*") computed by the
  [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md)
  function. Defaults to "*none*".

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

## Details

Implementation and terminology of the confidence intervals is based on
Hesterberg (2015) and Davison and Hinkley (1997) .

## References

Davison AC, Hinkley DV (1997). *Bootstrap Methods and their
Application*. Cambridge University Press.
[doi:10.1017/cbo9780511802843](https://doi.org/10.1017/cbo9780511802843)
.  
  
Hesterberg TC (2015). “What Teachers Should Know About the Bootstrap:
Resampling in the Undergraduate Statistics Curriculum.” *The American
Statistician*, **69**(4), 371–386.
[doi:10.1080/00031305.2015.1089789](https://doi.org/10.1080/00031305.2015.1089789)
.
