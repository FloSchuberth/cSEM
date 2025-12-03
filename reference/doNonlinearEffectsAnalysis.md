# Do a nonlinear effects analysis

**\[maturing\]**

## Usage

``` r
doNonlinearEffectsAnalysis(
 .object            = NULL,
 .dependent         = NULL, 
 .independent       = NULL,
 .moderator         = NULL,
 .n_steps           = 100,
 .values_moderator  = c(-2, -1, 0, 1, 2),
 .value_independent = 0,
 .alpha             = 0.05
 )
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .dependent:

  Character string. The name of the dependent variable.

- .independent:

  Character string. The name of the independent variable.

- .moderator:

  Character string. The name of the moderator variable.

- .n_steps:

  Integer. A value giving the number of steps (the spotlights, i.e.,
  values of .moderator in surface analysis or floodlight analysis)
  between the minimum and maximum value of the moderator. Defaults to
  `100`.

- .values_moderator:

  A numeric vector. The values of the moderator in a the simple effects
  analysis. Typically these are difference from the mean (=0) measured
  in standard deviations. Defaults to `c(-2, -1, 0, 1, 2)`.

- .value_independent:

  Integer. Only required for floodlight analysis; The value of the
  independent variable in case that it appears as a higher-order term.

- .alpha:

  An integer or a numeric vector of significance levels. Defaults to
  `0.05`.

## Value

A list of class `cSEMNonlinearEffects` with a corresponding method for
[`plot()`](https://rdrr.io/r/graphics/plot.default.html). See:
[`plot.cSEMNonlinearEffects()`](https://floschuberth.github.io/cSEM/reference/plot.cSEMNonlinearEffects.md).

## Details

Calculate the expected value of the dependent variable conditional on
the values of an independent variables and a moderator variable. All
other variables in the model are assumed to be zero, i.e., they are
fixed at their mean levels. Moreover, it produces the input for the
floodlight analysis.

## See also

[`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md),
[`plot.cSEMNonlinearEffects()`](https://floschuberth.github.io/cSEM/reference/plot.cSEMNonlinearEffects.md)

## Examples

``` r
if (FALSE) { # \dontrun{
model_Int <- "
# Measurement models
INV =~ INV1 + INV2 + INV3 +INV4
SAT =~ SAT1 + SAT2 + SAT3
INT =~ INT1 + INT2

# Structrual model containing an interaction term.
INT ~ INV + SAT + INV.SAT
"
  
# Estimate model
out <- csem(.data = Switching, .model = model_Int,
            # ADANCO settings
            .PLS_weight_scheme_inner = 'factorial',
            .tolerance = 1e-06,
            .resample_method = 'bootstrap'
)
  
# Do nonlinear effects analysis
neffects <- doNonlinearEffectsAnalysis(out, 
                                       .dependent = 'INT',
                                       .moderator = 'INV',
                                       .independent = 'SAT') 

# Get an overview
neffects

# Simple effects plot
plot(neffects, .plot_type = 'simpleeffects')

# Surface plot using plotly
plot(neffects, .plot_type = 'surface', .plot_package = 'plotly')

# Surface plot using persp
plot(neffects, .plot_type = 'surface', .plot_package = 'persp')

# Floodlight analysis
plot(neffects, .plot_type = 'floodlight')
} # }
```
