# Resample cSEMResults

Resample a
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object using bootstrap or jackknife resampling. The function is called
by [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md) if
the user sets `csem(..., .resample_method = "bootstrap")` or
`csem(..., .resample_method = "jackknife")` but may also be called
directly.

## Usage

``` r
resamplecSEMResults(
 .object                = NULL,
 .resample_method       = c("bootstrap", "jackknife"), 
 .resample_method2      = c("none", "bootstrap", "jackknife"), 
 .R                     = 499,
 .R2                    = 199,
 .handle_inadmissibles  = c("drop", "ignore", "replace"),
 .user_funs             = NULL,
 .eval_plan             = c("sequential", "multicore", "multisession"),
 .force                 = FALSE,
 .seed                  = NULL,
 .sign_change_option    = c("none","individual","individual_reestimate",
                            "construct_reestimate"),
 ...
)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .resample_method:

  Character string. The resampling method to use. One of: "*bootstrap*"
  or "*jackknife*". Defaults to "*bootstrap*".

- .resample_method2:

  Character string. The resampling method to use when resampling from a
  resample. One of: "*none*", "*bootstrap*" or "*jackknife*". For
  "*bootstrap*" the number of draws is provided via `.R2`. Currently,
  resampling from each resample is only required for the studentized
  confidence interval ("*CI_t_interval*") computed by the
  [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md)
  function. Defaults to "*none*".

- .R:

  Integer. The number of bootstrap replications. Defaults to `499`.

- .R2:

  Integer. The number of bootstrap replications to use when resampling
  from a resample. Defaults to `199`.

- .handle_inadmissibles:

  Character string. How should inadmissible results be treated? One of
  "*drop*", "*ignore*", or "*replace*". If "*drop*", all
  replications/resamples yielding an inadmissible result will be dropped
  (i.e. the number of results returned will potentially be less than
  `.R`). For "*ignore*" all results are returned even if all or some of
  the replications yielded inadmissible results (i.e. number of results
  returned is equal to `.R`). For "*replace*" resampling continues until
  there are exactly `.R` admissible solutions. Depending on the
  frequency of inadmissible solutions this may significantly increase
  computing time. Defaults to "*drop*".

- .user_funs:

  A function or a (named) list of functions to apply to every resample.
  The functions must take `.object` as its first argument (e.g.,
  `myFun <- function(.object, ...) {body-of-the-function}`). Function
  output should preferably be a (named) vector but matrices are also
  accepted. However, the output will be vectorized (columnwise) in this
  case. See the examples section for details.

- .eval_plan:

  Character string. The evaluation plan to use. One of "*sequential*",
  "*multicore*", or "*multisession*". In the two latter cases all
  available cores will be used. Defaults to "*sequential*".

- .force:

  Logical. Should .object be resampled even if it contains resamples
  already?. Defaults to `FALSE`.

- .seed:

  Integer or `NULL`. The random seed to use. Defaults to `NULL` in which
  case an arbitrary seed is chosen. Note that the scope of the seed is
  limited to the body of the function it is used in. Hence, the global
  seed will not be altered!

- .sign_change_option:

  Character string. Which sign change option should be used to handle
  flipping signs when resampling? One of "*none*","*individual*",
  "*individual_reestimate*", "*construct_reestimate*". Defaults to
  "*none*".

- ...:

  Further arguments passed to functions supplied to `.user_funs`.

## Value

The core structure is the same structure as that of `.object` with the
following elements added:

- `$Estimates_resamples`: A list containing the `.R` resamples and the
  original estimates for each of the resampled quantities
  (Path_estimates, Loading_estimates, Weight_estimates, user defined
  functions). Each list element is a list containing elements
  `$Resamples` and `$Original`. `$Resamples` is a `(.R x K)` matrix with
  each row representing one resample for each of the `K`
  parameters/statistics. `$Original` contains the original estimates
  (vectorized by column if the output of the user provided function is a
  matrix.

- `$Information_resamples`: A list containing additional information.

Use `str(<.object>, list.len = 3)` on the resulting object for an
overview.

## Details

Given `M` resamples (for bootstrap `M = .R` and for jackknife `M = N`,
where `N` is the number of observations) based on the data used to
compute the
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object provided via `.object`, `resamplecSEMResults()` essentially calls
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md) on
each resample using the arguments of the original call (ignoring any
arguments related to resampling) and returns estimates for each of a
subset of practically useful resampled parameters/statistics computed by
[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).
Currently, the following estimates are computed and returned by default
based on each resample: Path estimates, Loading estimates, Weight
estimates.

In practical application users may need to resample a specific statistic
(e.g, the heterotrait-monotrait ratio of correlations (HTMT) or
differences between path coefficients such as beta_1 - beta_2). Such
statistics may be provided by a function `fun(.object, ...)` or a list
of such functions via the `.user_funs` argument. The first argument of
these functions must always be `.object`. Internally, the function will
be applied on each resample to produce the desired statistic. Hence,
arbitrary complicated statistics may be resampled as long as the body of
the function draws on elements contained in the
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object only. Output of `fun(.object, ...)` should preferably be a
(named) vector but matrices are also accepted. However, the output will
be vectorized (columnwise) in this case. See the examples section for
details.

Both resampling the original
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object (call it "first resample") and resampling based on a resampled
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object (call it "second resample") are supported. Choices for the former
are "*bootstrap*" and "*jackknife*". Resampling based on a resample is
turned off by default (`.resample_method2 = "none"`) as this
significantly increases computation time (there are now `M * M2`
resamples to compute, where `M2` is `.R2` or `N`). Resamples of a
resample are required, e.g., for the studentized confidence interval
computed by the
[`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md)
function. Typically, bootstrap resamples are used in this case (Davison
and Hinkley 1997) .

As [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
accepts a single data set, a list of data sets as well as data sets that
contain a column name used to split the data into groups, the
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object may contain multiple data sets. In this case, resampling is done
by data set or group. Note that depending on the number of data
sets/groups, the computation may be considerably slower as resampling
will be repeated for each data set/group. However, apart from speed
considerations users don not need to worry about the type of input used
to compute the
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object as `resamplecSEMResults()` is able to deal with each case.

The number of bootstrap runs for the first and second run are given by
`.R` and `.R2`. The default is `499` for the first and `199` for the
second run but should be increased in real applications. See e.g.,
Hesterberg (2015) , p.380, Davison and Hinkley (1997) , and Efron and
Hastie (2016) for recommendations. For jackknife `.R` are `.R2` are
ignored.

Resampling may produce inadmissible results (as checked by
[`verify()`](https://floschuberth.github.io/cSEM/reference/verify.md)).
By default these results are dropped however users may choose to
`"ignore"` or `"replace"` inadmissible results in which resampling
continuous until the necessary number of admissible results is reached.

The cSEM package supports (multi)processing via the
[future](https://github.com/futureverse/future/) framework (Bengtsson
2018) . Users may simply choose an evaluation plan via `.eval_plan` and
the package takes care of all the complicated backend issues. Currently,
users may choose between standard single-core/single-session evaluation
(`"sequential"`) and multiprocessing (`"multisession"` or
`"multicore"`). The future package provides other options (e.g.,
`"cluster"` or `"remote"`), however, they probably will not be needed in
the context of the cSEM package as simulations usually do not require
high-performance clusters. Depending on the operating system, the future
package will manage to distribute tasks to multiple R sessions (Windows)
or multiple cores. Note that multiprocessing is not necessary always
faster when only a "small" number of replications is required as the
overhead of initializing new sessions or distributing tasks to different
cores will not immediately be compensated by the availability of
multiple sessions/cores.

Random number generation (RNG) uses the L'Ecuyer-CRMR RGN stream as
implemented in the [future.apply
package](https://github.com/futureverse/future.apply/) (Bengtsson 2018)
. It is independent of the evaluation plan. Hence, setting e.g.,
`.seed = 123` will generate the same random number and replicates for
both `.eval_plan = "sequential"`, `.eval_plan = "multisession"`, and
`.eval_plan = "multicore"`. See
[?future_lapply](https://future.apply.futureverse.org/reference/future_lapply.html)
for details.

## References

Bengtsson H (2018). *future: Unified Parallel and Distributed Processing
in R for Everyone*. R package version 1.10.0,
<https://CRAN.R-project.org/package=future>.  
  
Bengtsson H (2018). *future.apply: Apply Function to Elements in
Parallel using Futures*. R package version 1.0.1,
<https://CRAN.R-project.org/package=future.apply>.  
  
Davison AC, Hinkley DV (1997). *Bootstrap Methods and their
Application*. Cambridge University Press.
[doi:10.1017/cbo9780511802843](https://doi.org/10.1017/cbo9780511802843)
.  
  
Efron B, Hastie T (2016). *Computer Age Statistical Inference*.
Cambridge University Pr. ISBN 1107149894.  
  
Hesterberg TC (2015). “What Teachers Should Know About the Bootstrap:
Resampling in the Undergraduate Statistics Curriculum.” *The American
Statistician*, **69**(4), 371–386.
[doi:10.1080/00031305.2015.1089789](https://doi.org/10.1080/00031305.2015.1089789)
.

## See also

[csem](https://floschuberth.github.io/cSEM/reference/csem.md),
[`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md),
[`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Note: example not run as resampling is time consuming
# ===========================================================================
# Basic usage
# ===========================================================================
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

## Estimate the model without resampling 
a <- csem(satisfaction, model)

## Bootstrap and jackknife estimation
boot <- resamplecSEMResults(a)
jack <- resamplecSEMResults(a, .resample_method = "jackknife") 

## Alternatively use .resample_method in csem()
boot_csem <- csem(satisfaction, model, .resample_method = "bootstrap")
jack_csem <- csem(satisfaction, model, .resample_method = "jackknife")

# ===========================================================================
# Extended usage
# ===========================================================================
### Double resampling  ------------------------------------------------------
# The confidence intervals (e.g. the bias-corrected and accelearated CI) 
# require double resampling. Use .resample_method2 for this.

boot1 <- resamplecSEMResults(
  .object = a, 
  .resample_method  = "bootstrap", 
  .R                = 50,
  .resample_method2 = "bootstrap", 
  .R2               = 20,
  .seed             = 1303
  )

## Again, this is identical to using csem 
boot1_csem <- csem(
  .data             = satisfaction, 
  .model            = model, 
  .resample_method  = "bootstrap",
  .R                = 50,
  .resample_method2 = "bootstrap",
  .R2               = 20,
  .seed             = 1303
  )

identical(boot1, boot1_csem) # only true if .seed was set

### Inference ---------------------------------------------------------------
# To get inferencial quanitites such as the estimated standard error or
# the percentile confidence intervall for each resampled quantity use 
# postestimation function infer()

inference <- infer(boot1)
inference$Path_estimates$sd
inference$Path_estimates$CI_percentile

# As usual summarize() can be called directly
summarize(boot1)

# In the example above .R x .R2 = 50 x 20 = 1000. Multiprocessing will be
# faster on most systems here and is therefore recommended. Note that multiprocessing
# does not affect the random number generation

boot2 <- resamplecSEMResults(
  .object           = a, 
  .resample_method  = "bootstrap", 
  .R                = 50,
  .resample_method2 = "bootstrap", 
  .R2               = 20,
  .eval_plan        = "multisession", 
  .seed             = 1303
  )

identical(boot1, boot2)} # }
```
