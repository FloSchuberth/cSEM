# Resample data

Resample data from a data set using common resampling methods. For
bootstrap or jackknife resampling, package users usually do not need to
call this function but directly use
[`resamplecSEMResults()`](https://floschuberth.github.io/cSEM/reference/resamplecSEMResults.md)
instead.

## Usage

``` r
resampleData(
 .object          = NULL,
 .data            = NULL,
 .resample_method = c("bootstrap", "jackknife", "permutation", 
                      "cross-validation"),
 .cv_folds        = 10,
 .id              = NULL,
 .R               = 499,
 .seed            = NULL
)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .data:

  A `data.frame`, a `matrix` or a `list` of data of either type.
  Possible column types or classes of the data provided are:
  "`logical`", "`numeric`" ("`double`" or "`integer`"), "`factor`"
  (ordered and unordered) or a mix of several types. The data may also
  include **one** character column whose column name must be given to
  `.id`. This column is assumed to contain group identifiers used to
  split the data into groups. If `.data` is provided, `.object` is
  ignored. Defaults to `NULL`.

- .resample_method:

  Character string. The resampling method to use. One of: "*bootstrap*",
  "*jackknife*", "*permutation*", or "*cross-validation*". Defaults to
  "*bootstrap*".

- .cv_folds:

  Integer. The number of cross-validation folds to use. Setting
  `.cv_folds` to `N` (the number of observations) produces leave-one-out
  cross-validation samples. Defaults to `10`.

- .id:

  Character string or integer. A character string giving the name or an
  integer of the position of the column of `.data` whose levels are used
  to split `.data` into groups. Defaults to `NULL`.

- .R:

  Integer. The number of bootstrap runs, permutation runs or
  cross-validation repetitions to use. Defaults to `499`.

- .seed:

  Integer or `NULL`. The random seed to use. Defaults to `NULL` in which
  case an arbitrary seed is chosen. Note that the scope of the seed is
  limited to the body of the function it is used in. Hence, the global
  seed will not be altered!

## Value

The structure of the output depends on the type of input and the
resampling method:

- Bootstrap:

  If a `matrix` or `data.frame` without grouping variable is provided
  (i.e., `.id = NULL`), the result is a list of length `.R` (default
  `499`). Each element of that list is a bootstrap (re)sample. If a
  grouping variable is specified or a list of data is provided (where
  each list element is assumed to contain data for one group),
  resampling is done by group. Hence, the result is a list of length
  equal to the number of groups with each list element containing `.R`
  bootstrap samples based on the `N_g` observations of group `g`.

- Jackknife:

  If a `matrix` or `data.frame` without grouping variable is provided
  (`.id = NULL`), the result is a list of length equal to the number of
  observations/rows (`N`) of the data set provided. Each element of that
  list is a jackknife (re)sample. If a grouping variable is specified or
  a list of data is provided (where each list element is assumed to
  contain data for one group), resampling is done by group. Hence, the
  result is a list of length equal to the number of group levels with
  each list element containing `N` jackknife samples based on the `N_g`
  observations of group `g`.

- Permutation:

  If a `matrix` or `data.frame` without grouping variable is provided an
  error is returned as permutation will simply reorder the observations.
  If a grouping variable is specified or a list of data is provided
  (where each list element is assumed to contain data of one group),
  group membership is permuted. Hence, the result is a list of length
  `.R` where each element of that list is a permutation (re)sample.

- Cross-validation:

  If a `matrix` or `data.frame` without grouping variable is provided a
  list of length `.R` is returned. Each list element contains a list
  containing the `k` splits/folds subsequently used as test and training
  data sets. If a grouping variable is specified or a list of data is
  provided (where each list element is assumed to contain data for one
  group), cross-validation is repeated `.R` times for each group. Hence,
  the result is a list of length equal to the number of groups, each
  containing `.R` list elements (the repetitions) which in turn contain
  the `k` splits/folds.

## Details

The function `resampleData()` is general purpose. It simply resamples
data from a data set according to the resampling method provided via the
`.resample_method` argument and returns a list of resamples. Currently,
`bootstrap`, `jackknife`, `permutation`, and `cross-validation` (both
leave-one-out (LOOCV) and k-fold cross-validation) are implemented.

The user may provide the data set to resample either explicitly via the
`.data` argument or implicitly by providing a
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
objects to `.object` in which case the original data used in the call
that created the
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object is used for resampling. If both, a
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object and a data set via `.data` are provided the former is ignored.

As [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
accepts a single data set, a list of data sets as well as data sets that
contain a column name used to split the data into groups, the
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object may contain multiple data sets. In this case, resampling is done
by data set or group. Note that depending on the number of data
sets/groups provided this computation may be slower as resampling will
be repeated for each data set/group.

To split data provided via the `.data` argument into groups, the column
name or the column index of the column containing the group levels to
split the data must be given to `.id`. If data that contains grouping is
taken from a
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
object, `.id` is taken from the object information. Hence, providing
`.id` is redundant in this case and therefore ignored.

The number of bootstrap or permutation runs as well as the number of
cross-validation repetitions is given by `.R`. The default is `499` but
should be increased in real applications. See e.g., Hesterberg (2015) ,
p.380 for recommendations concerning the bootstrap. For jackknife `.R`
is ignored as it is based on the N leave-one-out data sets.

Choosing `resample_method = "permutation"` for ungrouped data causes an
error as permutation will simply reorder the observations which is
usually not meaningful. If a list of data is provided each list element
is assumed to represent the observations belonging to one group. In this
case, data is pooled and group adherence permuted.

For cross-validation the number of folds (`k`) defaults to `10`. It may
be changed via the `.cv_folds` argument. Setting `k = 2` (not 1!) splits
the data into a single training and test data set. Setting `k = N`
(where `N` is the number of observations) produces leave-one-out
cross-validation samples. Note: 1.) At least 2 folds required (`k > 1`);
2.) `k` can not be larger than `N`; 3.) If `N/k` is not not an integer
the last fold will have less observations.

Random number generation (RNG) uses the L'Ecuyer-CRMR RGN stream as
implemented in the [future.apply
package](https://github.com/futureverse/future.apply/) (Bengtsson 2018)
. See
[?future_lapply](https://future.apply.futureverse.org/reference/future_lapply.html)
for details. By default a random seed is chosen.

## References

Bengtsson H (2018). *future.apply: Apply Function to Elements in
Parallel using Futures*. R package version 1.0.1,
<https://CRAN.R-project.org/package=future.apply>.  
  
Hesterberg TC (2015). “What Teachers Should Know About the Bootstrap:
Resampling in the Undergraduate Statistics Curriculum.” *The American
Statistician*, **69**(4), 371–386.
[doi:10.1080/00031305.2015.1089789](https://doi.org/10.1080/00031305.2015.1089789)
.

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md),
[`resamplecSEMResults()`](https://floschuberth.github.io/cSEM/reference/resamplecSEMResults.md)

## Examples

``` r
# ===========================================================================
# Using the raw data 
# ===========================================================================
### Bootstrap (default) -----------------------------------------------------

res_boot1 <- resampleData(.data = satisfaction)
str(res_boot1, max.level = 3, list.len = 3)
#> List of 499
#>  $ :'data.frame':    250 obs. of  27 variables:
#>   ..$ imag1: int [1:250] 8 9 8 8 8 9 5 10 7 7 ...
#>   ..$ imag2: int [1:250] 9 9 9 9 7 9 5 10 9 7 ...
#>   ..$ imag3: int [1:250] 9 9 9 9 8 9 7 10 9 8 ...
#>   .. [list output truncated]
#>  $ :'data.frame':    250 obs. of  27 variables:
#>   ..$ imag1: int [1:250] 9 8 7 9 8 9 9 5 6 8 ...
#>   ..$ imag2: int [1:250] 10 8 10 7 5 9 8 6 6 5 ...
#>   ..$ imag3: int [1:250] 9 9 10 7 6 10 8 5 7 9 ...
#>   .. [list output truncated]
#>  $ :'data.frame':    250 obs. of  27 variables:
#>   ..$ imag1: int [1:250] 7 4 4 8 8 6 7 8 8 9 ...
#>   ..$ imag2: int [1:250] 8 6 5 9 10 7 9 7 6 7 ...
#>   ..$ imag3: int [1:250] 8 9 5 8 9 7 9 8 7 7 ...
#>   .. [list output truncated]
#>   [list output truncated]

## To replicate a bootstrap draw use .seed:
res_boot1a <- resampleData(.data = satisfaction, .seed = 2364)
res_boot1b <- resampleData(.data = satisfaction, .seed = 2364)
                           
identical(res_boot1, res_boot1a) # TRUE
#> [1] FALSE

### Jackknife ---------------------------------------------------------------

res_jack <- resampleData(.data = satisfaction, .resample_method = "jackknife")
str(res_jack, max.level = 3, list.len = 3)
#> List of 250
#>  $ :'data.frame':    249 obs. of  27 variables:
#>   ..$ imag1: int [1:249] 9 9 8 10 7 5 9 9 10 2 ...
#>   ..$ imag2: int [1:249] 9 8 9 10 8 5 9 8 10 2 ...
#>   ..$ imag3: int [1:249] 10 8 8 8 8 5 9 9 10 8 ...
#>   .. [list output truncated]
#>  $ :'data.frame':    249 obs. of  27 variables:
#>   ..$ imag1: int [1:249] 8 9 8 10 7 5 9 9 10 2 ...
#>   ..$ imag2: int [1:249] 8 8 9 10 8 5 9 8 10 2 ...
#>   ..$ imag3: int [1:249] 9 8 8 8 8 5 9 9 10 8 ...
#>   .. [list output truncated]
#>  $ :'data.frame':    249 obs. of  27 variables:
#>   ..$ imag1: int [1:249] 8 9 8 10 7 5 9 9 10 2 ...
#>   ..$ imag2: int [1:249] 8 9 9 10 8 5 9 8 10 2 ...
#>   ..$ imag3: int [1:249] 9 10 8 8 8 5 9 9 10 8 ...
#>   .. [list output truncated]
#>   [list output truncated]

### Cross-validation --------------------------------------------------------
## Create dataset for illustration:
dat <- data.frame(
  "x1" = rnorm(100),
  "x2" = rnorm(100),
  "group" = sample(c("male", "female"), size = 100, replace = TRUE),
  stringsAsFactors = FALSE)

## 10-fold cross-validation (repeated 100 times)
cv_10a <- resampleData(.data = dat, .resample_method = "cross-validation", 
                      .R = 100)
str(cv_10a, max.level = 3, list.len = 3)
#> List of 100
#>  $ :List of 10
#>   ..$ 1 :'data.frame':   10 obs. of  3 variables:
#>   .. ..$ x1   : num [1:10] -0.9456 0.0377 0.6029 -1.5955 -0.2657 ...
#>   .. ..$ x2   : num [1:10] 0.47 0.192 0.449 -0.265 0.322 ...
#>   .. ..$ group: chr [1:10] "female" "female" "female" "female" ...
#>   ..$ 2 :'data.frame':   10 obs. of  3 variables:
#>   .. ..$ x1   : num [1:10] -0.779 -0.264 -0.6259 0.8517 -0.0193 ...
#>   .. ..$ x2   : num [1:10] 0.9751 -1.1964 -0.0422 -0.9473 -0.4089 ...
#>   .. ..$ group: chr [1:10] "male" "female" "female" "male" ...
#>   ..$ 3 :'data.frame':   10 obs. of  3 variables:
#>   .. ..$ x1   : num [1:10] -1.428 -0.464 -1.568 -0.766 -0.686 ...
#>   .. ..$ x2   : num [1:10] -0.284 -0.535 2.708 -0.386 -0.982 ...
#>   .. ..$ group: chr [1:10] "male" "male" "male" "male" ...
#>   .. [list output truncated]
#>  $ :List of 10
#>   ..$ 1 :'data.frame':   10 obs. of  3 variables:
#>   .. ..$ x1   : num [1:10] 0.1516 0.0845 -0.2607 0.2064 0.7902 ...
#>   .. ..$ x2   : num [1:10] -0.284 0.786 0.286 1.901 1.227 ...
#>   .. ..$ group: chr [1:10] "male" "female" "female" "male" ...
#>   ..$ 2 :'data.frame':   10 obs. of  3 variables:
#>   .. ..$ x1   : num [1:10] -0.489 -0.628 0.496 0.486 -1.277 ...
#>   .. ..$ x2   : num [1:10] 0.709 0.736 1.373 -0.462 -1.472 ...
#>   .. ..$ group: chr [1:10] "female" "male" "male" "male" ...
#>   ..$ 3 :'data.frame':   10 obs. of  3 variables:
#>   .. ..$ x1   : num [1:10] 0.6286 -0.0193 -0.242 0.9618 1.3766 ...
#>   .. ..$ x2   : num [1:10] -0.1395 -0.4089 0.6667 -0.0325 1.0919 ...
#>   .. ..$ group: chr [1:10] "male" "female" "female" "female" ...
#>   .. [list output truncated]
#>  $ :List of 10
#>   ..$ 1 :'data.frame':   10 obs. of  3 variables:
#>   .. ..$ x1   : num [1:10] -0.752 -0.242 0.163 -0.375 0.903 ...
#>   .. ..$ x2   : num [1:10] 0.463 1.066 -0.576 -0.868 2.454 ...
#>   .. ..$ group: chr [1:10] "male" "female" "female" "male" ...
#>   ..$ 2 :'data.frame':   10 obs. of  3 variables:
#>   .. ..$ x1   : num [1:10] -0.152 -1.261 -0.393 -0.489 1.764 ...
#>   .. ..$ x2   : num [1:10] 2.011 -1.649 0.499 0.709 0.37 ...
#>   .. ..$ group: chr [1:10] "male" "male" "female" "female" ...
#>   ..$ 3 :'data.frame':   10 obs. of  3 variables:
#>   .. ..$ x1   : num [1:10] -0.2289 0.0845 -0.7555 0.611 -0.9456 ...
#>   .. ..$ x2   : num [1:10] -1.178 0.786 0.489 -1.405 0.47 ...
#>   .. ..$ group: chr [1:10] "female" "female" "female" "male" ...
#>   .. [list output truncated]
#>   [list output truncated]

# Cross-validation can be done by group if a group identifyer is provided:
cv_10 <- resampleData(.data = dat, .resample_method = "cross-validation", 
                      .id = "group", .R = 100)
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: data length is not a multiple of split variable
#> Warning: UNRELIABLE VALUE: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'future.seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced via a parallel RNG method. To disable this check, use 'future.seed = NULL', or set option 'future.rng.onMisuse' to "ignore". [future ‘future_lapply-1’ (d248ad6633a3a5683c357df86b40215d-15); on d248ad6633a3a5683c357df86b40215d@runnervmeorf1<7540>]

## Leave-one-out-cross-validation (repeated 50 times)
cv_loocv  <- resampleData(.data = dat[, -3], 
                          .resample_method = "cross-validation", 
                          .cv_folds = nrow(dat),
                          .R = 50)
str(cv_loocv, max.level = 2, list.len = 3)
#> List of 50
#>  $ :List of 100
#>   ..$ 1  :'data.frame':  1 obs. of  2 variables:
#>   ..$ 2  :'data.frame':  1 obs. of  2 variables:
#>   ..$ 3  :'data.frame':  1 obs. of  2 variables:
#>   .. [list output truncated]
#>  $ :List of 100
#>   ..$ 1  :'data.frame':  1 obs. of  2 variables:
#>   ..$ 2  :'data.frame':  1 obs. of  2 variables:
#>   ..$ 3  :'data.frame':  1 obs. of  2 variables:
#>   .. [list output truncated]
#>  $ :List of 100
#>   ..$ 1  :'data.frame':  1 obs. of  2 variables:
#>   ..$ 2  :'data.frame':  1 obs. of  2 variables:
#>   ..$ 3  :'data.frame':  1 obs. of  2 variables:
#>   .. [list output truncated]
#>   [list output truncated]

### Permuation ---------------------------------------------------------------

res_perm <- resampleData(.data = dat, .resample_method = "permutation",
                         .id = "group")
str(res_perm, max.level = 2, list.len = 3)
#> List of 499
#>  $ :'data.frame':    100 obs. of  3 variables:
#>   ..$ x1: num [1:100] -0.242 -0.468 -0.773 2.15 -1.334 ...
#>   ..$ x2: num [1:100] 1.066 -0.536 0.536 -1.829 -1.814 ...
#>   ..$ id: chr [1:100] "male" "male" "female" "female" ...
#>  $ :'data.frame':    100 obs. of  3 variables:
#>   ..$ x1: num [1:100] -0.242 -0.468 -0.773 2.15 -1.334 ...
#>   ..$ x2: num [1:100] 1.066 -0.536 0.536 -1.829 -1.814 ...
#>   ..$ id: chr [1:100] "male" "female" "male" "male" ...
#>  $ :'data.frame':    100 obs. of  3 variables:
#>   ..$ x1: num [1:100] -0.242 -0.468 -0.773 2.15 -1.334 ...
#>   ..$ x2: num [1:100] 1.066 -0.536 0.536 -1.829 -1.814 ...
#>   ..$ id: chr [1:100] "male" "female" "female" "male" ...
#>   [list output truncated]

# Forgetting to set .id causes an error
if (FALSE) { # \dontrun{
res_perm <- resampleData(.data = dat, .resample_method = "permutation")
} # }

# ===========================================================================
# Using a cSEMResults object
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
a <- csem(satisfaction, model)

# Create bootstrap and jackknife samples
res_boot <- resampleData(a, .resample_method = "bootstrap", .R = 499)
res_jack <- resampleData(a, .resample_method = "jackknife")

# Since `satisfaction` is the dataset used the following approaches yield
# identical results.
res_boot_data   <- resampleData(.data = satisfaction, .seed = 2364)
res_boot_object <- resampleData(a, .seed = 2364)

identical(res_boot_data, res_boot_object) # TRUE
#> [1] TRUE
```
