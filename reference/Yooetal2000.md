# Data: Yooetal2000

A data frame containing 34 variables with 569 observations.

## Usage

``` r
Yooetal2000
```

## Format

An object of class `data.frame` with 569 rows and 34 columns.

## Source

Simulated data with the same correlation matrix as the data studied by
Yoo et al. (2000) .

## Details

The data is simulated and has the identical correlation matrix as the
data that was analysed by Yoo et al. (2000) to examine how five elements
of the marketing mix, namely price, store image, distribution intensity,
advertising spending, and price deals, are related to the so-called
dimensions of brand equity, i.e., perceived brand quality, brand
loyalty, and brand awareness/associations. It is also used in Henseler
(2017) and Henseler (2021) for demonstration purposes, see the
corresponding tutorial.

## References

Henseler J (2017). “Bridging Design and Behavioral Research With
Variance-Based Structural Equation Modeling.” *Journal of Advertising*,
**46**(1), 178–192.
[doi:10.1080/00913367.2017.1281780](https://doi.org/10.1080/00913367.2017.1281780)
.  
  
Henseler J (2021). *Composite-Based Structural Equation Modeling:
Analyzing Latent and Emergent Variables*. Guilford Press, New York.  
  
Yoo B, Donthu N, Lee S (2000). “An Examination of Selected Marketing Mix
Elements and Brand Equity.” *Journal of the Academy of Marketing
Science*, **28**(2), 195–211.
[doi:10.1177/0092070300282002](https://doi.org/10.1177/0092070300282002)
.

## Examples

``` r
#============================================================================
# Example is taken from Henseler (2021)
#============================================================================
model_HOC="
# Measurement models FOC
PR =~ PR1 + PR2 + PR3
IM =~ IM1 + IM2 + IM3
DI =~ DI1 + DI2 + DI3
AD =~ AD1 + AD2 + AD3
DL =~ DL1 + DL2 + DL3
AA =~ AA1 + AA2 + AA3 + AA4 + AA5 + AA6
LO =~ LO1 + LO3
QL =~ QL1 + QL2 + QL3 + QL4 + QL5 + QL6

# Composite model for SOC
BR <~ QL + LO + AA

# Structural model
BR~ PR + IM + DI + AD + DL 
"

out <- csem(.data = Yooetal2000, .model = model_HOC,
            .PLS_weight_scheme_inner = 'factorial',
            .tolerance = 1e-06)
```
