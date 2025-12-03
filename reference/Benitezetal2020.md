# Data: Benitezetal2020

A data frame containing 22 variables with 300 observations.

## Usage

``` r
Benitezetal2020
```

## Format

An object of class `data.frame` with 300 rows and 22 columns.

## Source

The dataset is provided as supplementary material by Benitez et al.
(2020) .

## Details

The simulated data contains variables about the social executive and
employee behavior. Moreover, it contains variables about the social
media capability and business performance. The dataset was used as an
illustrative example in Benitez et al. (2020) .

## References

Benitez J, Henseler J, Castillo A, Schuberth F (2020). “How to perform
and report an impactful analysis using partial least squares: Guidelines
for confirmatory and explanatory IS research.” *Information &
Management*, **2**(57), 103168.
[doi:10.1016/j.im.2019.05.003](https://doi.org/10.1016/j.im.2019.05.003)
.

## Examples

``` r
#============================================================================
# Example is taken from Benitez et al. (2020)
#============================================================================
model_Benitez <-"
# Reflective measurement models# Reflective measurement models
SEXB =~ SEXB1 + SEXB2 + SEXB3 +SEXB4
SEMB =~ SEMB1 + SEMB2 + SEMB3 + SEMB4

# Composite models
SMC <~ SMC1 + SMC2 + SMC3 + SMC4
BPP <~ BPP1 + BPP2 + BPP3 + BPP4 + BPP5

# Control variables
FS<~ FirmSize
Ind <~ Industry1 + Industry2 + Industry3

# Structural model
SMC ~ SEXB + SEMB 
BPP ~ SMC + Ind + FS
"

out <- csem(.data = Benitezetal2020, .model = model_Benitez,
            .PLS_weight_scheme_inner = 'factorial',
            .tolerance = 1e-06)
```
