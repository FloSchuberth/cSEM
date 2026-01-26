# Data: threecommonfactors

A dataset containing 500 standardized observations on 9 indicator
generated from a population model with three concepts modeled as common
factors.

## Usage

``` r
threecommonfactors
```

## Format

A matrix with 500 rows and 9 variables:

- y11-y13:

  Indicators attached to the first common factor (`eta1`). Population
  loadings are: 0.7; 0.7; 0.7

- y21-y23:

  Indicators attached to the second common factor (`eta2`). Population
  loadings are: 0.5; 0.7; 0.8

- y31-y33:

  Indicators attached to the third common factor (`eta3`). Population
  loadings are: 0.8; 0.75; 0.7

The model is: \$\$\`eta2\` = gamma1 \* \`eta1\` + zeta1\$\$ \$\$\`eta3\`
= gamma2 \* \`eta1\` + beta \* \`eta2\` + zeta2\$\$

with population values `gamma1` = 0.6, `gamma2` = 0.4 and `beta` = 0.35.

## Examples

``` r
#============================================================================
# Correct model (the model used to generate the data)
#============================================================================
model_correct <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# Measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33 
"

a <- csem(threecommonfactors, model_correct)

## The overall model fit is evidently almost perfect:
testOMF(a, .R = 30) # .R = 30 to speed up the example
#> ________________________________________________________________________________
#> --------- Test for overall model fit based on Beran & Srivastava (1985) --------
#> 
#> Null hypothesis:
#> 
#>        ┌──────────────────────────────────────────────────────────────────┐
#>        │                                                                  │
#>        │   H0: The model-implied indicator covariance matrix equals the   │
#>        │   population indicator covariance matrix.                        │
#>        │                                                                  │
#>        └──────────────────────────────────────────────────────────────────┘
#> 
#> Test statistic and critical value: 
#> 
#>                                      Critical value
#>  Distance measure    Test statistic    95%   
#>  dG                      0.0060      0.0200  
#>  SRMR                    0.0158      0.0287  
#>  dL                      0.0112      0.0371  
#>  dML                     0.0320      0.1033  
#>  
#> 
#> Decision: 
#> 
#>                          Significance level
#>  Distance measure             95%        
#>  dG                      Do not reject  
#>  SRMR                    Do not reject  
#>  dL                      Do not reject  
#>  dML                     Do not reject  
#>  
#> Additional information:
#> 
#>  Out of 30 bootstrap replications 30 are admissible.
#>  See ?verify() for what constitutes an inadmissible result.
#> 
#>  The seed used was: -703360890
#> ________________________________________________________________________________
```
