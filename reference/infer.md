# Inference

**\[stable\]**

## Usage

``` r
infer(
 .object            = NULL,
 .quantity          = c("all", "mean", "sd", "bias", "CI_standard_z", 
                        "CI_standard_t", "CI_percentile", "CI_basic", 
                        "CI_bc", "CI_bca", "CI_t_interval"),
 .alpha             = 0.05,
 .bias_corrected    = TRUE
)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .quantity:

  Character string. Which statistic should be returned? One of "*all*",
  "*mean*", "*sd*", "*bias*", "*CI_standard_z*", "*CI_standard_t*",
  "*CI_percentile*", "*CI_basic*", "*CI_bc*", "*CI_bca*",
  "*CI_t_interval*" Defaults to "*all*" in which case all quantities
  that do not require additional resampling are returned, i.e., all
  quantities but "*CI_bca*", "*CI_t_interval*".

- .alpha:

  An integer or a numeric vector of significance levels. Defaults to
  `0.05`.

- .bias_corrected:

  Logical. Should the standard and the tStat confidence interval be
  bias-corrected using the bootstrapped bias estimate? If `TRUE` the
  confidence interval for some estimated parameter `theta` is centered
  at `2*theta - theta*_hat`, where `theta*_hat` is the average over all
  `.R` bootstrap estimates of `theta`. Defaults to `TRUE`

## Value

A list of class `cSEMInfer`.

## Details

Calculate common inferential quantities. For users interested in the
estimated standard errors, t-values, p-values and/or confidences
intervals of the path, weight or loading estimates, calling
[`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
directly will usually be more convenient as it has a much more
user-friendly print method. `infer()` is useful for comparing different
confidence interval estimates.

`infer()` is a convenience wrapper around a number of internal functions
that compute a particular inferential quantity, i.e., a value or set of
values to be used in statistical inference.

cSEM relies on resampling (bootstrap and jackknife) as the basis for the
computation of e.g., standard errors or confidence intervals.
Consequently, `infer()` requires resamples to work. Technically, the
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object used in the call to `infer()` must therefore also have class
attribute `cSEMResults_resampled`. If the object provided by the user
does not contain resamples yet, `infer()` will obtain bootstrap
resamples first. Naturally, computation will take longer in this case.

`infer()` does as much as possible in the background. Hence, every time
`infer()` is called on a
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object the quantities chosen by the user are automatically computed for
every estimated parameter contained in the object. By default all
possible quantities are computed (`.quantity = all`). The following
table list the available inferential quantities alongside a brief
description. Implementation and terminology of the confidence intervals
is based on Hesterberg (2015) and Davison and Hinkley (1997) .

- `"mean"`, `"sd"`:

  The mean or the standard deviation over all `M` resample estimates of
  a generic statistic or parameter.

- `"bias"`:

  The difference between the resample mean and the original estimate of
  a generic statistic or parameter.

- `"CI_standard_z"` and `"CI_standard_t"`:

  The standard confidence interval for a generic statistic or parameter
  with standard errors estimated by the resample standard deviation.
  While `"CI_standard_z"` assumes a standard normally distributed
  statistic, `"CI_standard_t"` assumes a t-statistic with N - 1 degrees
  of freedom.

- `"CI_percentile"`:

  The percentile confidence interval. The lower and upper bounds of the
  confidence interval are estimated as the alpha and 1-alpha quantiles
  of the distribution of the resample estimates.

- `"CI_basic"`:

  The basic confidence interval also called the reverse bootstrap
  percentile confidence interval. See Hesterberg (2015) for details.

- `"CI_bc"`:

  The bias corrected (Bc) confidence interval. See Davison and
  Hinkley (1997) for details.

- `"CI_bca"`:

  The bias-corrected and accelerated (Bca) confidence interval. Requires
  additional jackknife resampling to compute the influence values. See
  Davison and Hinkley (1997) for details.

- `"CI_t_interval"`:

  The "studentized" t-confidence interval. If based on bootstrap
  resamples the interval is also called the bootstrap t-interval
  confidence interval. See Hesterberg (2015) on page 381. Requires
  resamples of resamples. See
  [`resamplecSEMResults()`](https://floschuberth.github.io/cSEM/reference/resamplecSEMResults.md).

By default, all but the studendized t-interval confidence interval and
the bias-corrected and accelerated confidence interval are calculated.
The reason for excluding these quantities by default are that both
require an additional resampling step. The former requires jackknife
estimates to compute influence values and the latter requires double
bootstrap. Both can potentially be time consuming. Hence, computation is
triggered only if explicitly chosen.

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

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[`resamplecSEMResults()`](https://floschuberth.github.io/cSEM/reference/resamplecSEMResults.md),
[`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)

## Examples

``` r
model <- "
# Structural model
QUAL ~ EXPE
EXPE ~ IMAG
SAT  ~ IMAG + EXPE + QUAL + VAL
LOY  ~ IMAG + SAT
VAL  ~ EXPE + QUAL

# Measurement model
EXPE =~ expe1 + expe2 + expe3 + expe4 + expe5
IMAG =~ imag1 + imag2 + imag3 + imag4 + imag5
LOY  =~ loy1  + loy2  + loy3  + loy4
QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
SAT  =~ sat1  + sat2  + sat3  + sat4
VAL  =~ val1  + val2  + val3  + val4
"
  
## Estimate the model with bootstrap resampling 
a <- csem(satisfaction, model, .resample_method = "bootstrap", .R = 20,
          .handle_inadmissibles = "replace")

## Compute inferential quantities
inf <- infer(a)

inf$Path_estimates$CI_basic
#>      EXPE ~ IMAG QUAL ~ EXPE VAL ~ EXPE VAL ~ QUAL SAT ~ IMAG SAT ~ EXPE
#> 95%L   0.6167111    1.018643  -5.536211   4.632356 -0.3189832  0.9848302
#> 95%U   0.7677769    1.074531  -3.726635   6.455080  0.2976096  2.8558962
#>      SAT ~ QUAL SAT ~ VAL LOY ~ IMAG LOY ~ SAT
#> 95%L  -4.198745  1.080635  0.1970053 0.2499581
#> 95%U  -1.304647  2.670981  0.6463820 0.6286743
inf$Indirect_effect$sd
#> QUAL ~ IMAG  VAL ~ IMAG  VAL ~ EXPE  SAT ~ IMAG  SAT ~ EXPE  SAT ~ QUAL 
#>  0.04972580  0.05231434  0.58354380  0.09482609  0.41358789  1.00710014 
#>  LOY ~ IMAG  LOY ~ EXPE  LOY ~ QUAL   LOY ~ VAL 
#>  0.08212706  0.09613538  0.31373557  0.24256665 

### Compute the bias-corrected and accelerated and/or the studentized t-inverval.
## For the studentied t-interval confidence interval a double bootstrap is required.
## This is pretty time consuming.
if (FALSE) { # \dontrun{
  inf <- infer(a, .quantity = c("all", "CI_bca")) # requires jackknife estimates 
  
## Estimate the model with double bootstrap resampling:
# Notes:
#   1. The .resample_method2 arguments triggers a bootstrap of each bootstrap sample
#   2. The double bootstrap is is very time consuming, consider setting 
#      `.eval_plan = "multisession`. 
a1 <- csem(satisfaction, model, .resample_method = "bootstrap", .R = 499,
          .resample_method2 = "bootstrap", .R2 = 199, .handle_inadmissibles = "replace") 
infer(a1, .quantity = "CI_t_interval")} # }
```
