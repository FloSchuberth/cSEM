# Test measurement invariance of composites

**\[stable\]**

## Usage

``` r
testMICOM(
 .object               = NULL,
 .approach_p_adjust    = "none",
 .handle_inadmissibles = c("drop", "ignore", "replace"), 
 .R                    = 499,
 .seed                 = NULL,
 .verbose              = TRUE
 )
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .approach_p_adjust:

  Character string or a vector of character strings. Approach used to
  adjust the p-value for multiple testing. See the `methods` argument of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for a
  list of choices and their description. Defaults to "*none*".

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

- .R:

  Integer. The number of bootstrap replications. Defaults to `499`.

- .seed:

  Integer or `NULL`. The random seed to use. Defaults to `NULL` in which
  case an arbitrary seed is chosen. Note that the scope of the seed is
  limited to the body of the function it is used in. Hence, the global
  seed will not be altered!

- .verbose:

  Logical. Should information (e.g., progress bar) be printed to the
  console? Defaults to `TRUE`.

## Value

A named list of class `cSEMTestMICOM` containing the following list
element:

- `$Step2`:

  A list containing the results of the test for compositional invariance
  (Step 2).

- `$Step3`:

  A list containing the results of the test for mean and variance
  equality (Step 3).

- `$Information`:

  A list of additional information on the test.

## Details

The functions performs the permutation-based test for measurement
invariance of composites across groups proposed by Henseler et al.
(2016) . According to the authors assessing measurement invariance in
composite models can be assessed by a three-step procedure. The first
two steps involve an assessment of configural and compositional
invariance. The third steps involves mean and variance comparisons
across groups. Assessment of configural invariance is qualitative in
nature and hence not assessed by the `testMICOM()` function.

As `testMICOM()` requires at least two groups, `.object` must be of
class `cSEMResults_multi`. As of version 0.2.0 of the package,
`testMICOM()` does not support models containing second-order
constructs.

It is possible to compare more than two groups, however,
multiple-testing issues arise in this case. To adjust p-values in this
case several p-value adjustments are available via the
`approach_p_adjust` argument.

The remaining arguments set the number of permutation runs to conduct
(`.R`), the random number seed (`.seed`), instructions how inadmissible
results are to be handled (`handle_inadmissibles`), and whether the
function should be verbose in a sense that progress is printed to the
console.

The number of permutation runs defaults to `args_default()$.R` for
performance reasons. According to Henseler et al. (2016) the number of
permutations should be at least 5000 for assessment to be sufficiently
reliable.

## References

Henseler J, Ringle CM, Sarstedt M (2016). “Testing Measurement
Invariance of Composites Using Partial Least Squares.” *International
Marketing Review*, **33**(3), 405–431.
[doi:10.1108/imr-09-2014-0304](https://doi.org/10.1108/imr-09-2014-0304)
.

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md),
[`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md),
[`testMGD()`](https://floschuberth.github.io/cSEM/reference/testMGD.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# NOTE: to run the example. Download and load the newest version of cSEM.DGP
# from GitHub using devtools::install_github("M-E-Rademaker/cSEM.DGP").

# Create two data generating processes (DGPs) that only differ in how the composite
# X is build. Hence, the two groups are not compositionally invariant.
dgp1 <- "
# Structural model
Y ~ 0.6*X

# Measurement model
Y =~ 1*y1
X <~ 0.4*x1 + 0.8*x2

x1 ~~ 0.3125*x2
"

dgp2 <- "
# Structural model
Y ~ 0.6*X

# Measurement model
Y =~ 1*y1
X <~ 0.8*x1 + 0.4*x2

x1 ~~ 0.3125*x2
"

g1 <- generateData(dgp1, .N = 399, .empirical = TRUE) # requires cSEM.DGP 
g2 <- generateData(dgp2, .N = 200, .empirical = TRUE) # requires cSEM.DGP

# Model is the same for both DGPs
model <- "
# Structural model
Y ~ X

# Measurement model
Y =~ y1
X <~ x1 + x2
"

# Estimate
csem_results <- csem(.data = list("group1" = g1, "group2" = g2), model)

# Test
testMICOM(csem_results, .R = 50, .alpha = c(0.01, 0.05), .seed = 1987)
} # }
```
