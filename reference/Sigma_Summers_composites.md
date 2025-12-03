# Data: Summers

A (18 x 18) indicator correlation matrix.

## Usage

``` r
Sigma_Summers_composites
```

## Format

An object of class `matrix` (inherits from `array`) with 18 rows and 18
columns.

## Source

Own calculation based on Dijkstra and Henseler (2015) .

## Details

The indicator correlation matrix for a modified version of Summers
(1965) model. All constructs are modeled as composites.

## References

Dijkstra TK, Henseler J (2015). “Consistent and Asymptotically Normal
PLS Estimators for Linear Structural Equations.” *Computational
Statistics & Data Analysis*, **81**, 10–23.  
  
Summers R (1965). “A Capital Intensive Approach to the Small Sample
Properties of Various Simultaneous Equation Estimators.” *Econometrica*,
**33**(1), 1–41.

## Examples

``` r
require(cSEM)

model <- "
ETA1 ~ ETA2 + XI1 + XI2
ETA2 ~ ETA1 + XI3 +XI4

ETA1 ~~ ETA2

XI1  <~ x1 + x2 + x3
XI2  <~ x4 + x5 + x6
XI3  <~ x7 + x8 + x9
XI4  <~ x10 + x11 + x12
ETA1 <~ y1 + y2 + y3
ETA2 <~ y4 + y5 + y6
"

## Generate data
summers_dat <- MASS::mvrnorm(n = 300, mu = rep(0, 18), 
                             Sigma = Sigma_Summers_composites, empirical = TRUE)

## Estimate
res <- csem(.data = summers_dat, .model = model) # inconsistent

## 
# 2SLS
res_2SLS <- csem(.data = summers_dat, .model = model, .approach_paths = "2SLS",
                 .instruments = list(ETA1 = c('XI1', 'XI2', 'XI3', 'XI4'),
                                     ETA2 = c('XI1', 'XI2', 'XI3', 'XI4'))
)
```
