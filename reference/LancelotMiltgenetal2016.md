# Data: LancelotMiltgenetal2016

A data frame containing 10 variables with 1090 observations.

## Usage

``` r
LancelotMiltgenetal2016
```

## Format

An object of class `data.frame` with 1090 rows and 11 columns.

## Source

This data has been collected through a cooperation with the European
Commission Joint Research Center Institute for Prospective Technological
Studies, contract “Young People and Emerging Digital Services: An
Exploratory Survey on Motivations, Perceptions, and Acceptance of Risk”
(EC JRC Contract IPTS No: 150876-2007 F1ED-FR).

## Details

The data was analysed by Lancelot-Miltgen et al. (2016) to study young
consumers’ adoption intentions of a location tracker technology in the
light of privacy concerns. It is also used in Henseler (2021) for
demonstration purposes, see the corresponding tutorial.

## References

Henseler J (2021). *Composite-Based Structural Equation Modeling:
Analyzing Latent and Emergent Variables*. Guilford Press, New York.  
  
Lancelot-Miltgen C, Henseler J, Gelhard C, Popovic A (2016).
“Introducing new products that affect consumer privacy: A mediation
model.” *Journal of Business Research*, **69**(10), 4659–4666.
[doi:10.1016/j.jbusres.2016.04.015](https://doi.org/10.1016/j.jbusres.2016.04.015)
.

## Examples

``` r
#============================================================================
# Example is taken from Henseler (2020)
#============================================================================
model_Med <- "
# Reflective measurement model
Trust =~ trust1 + trust2
PrCon =~ privcon1 + privcon2 + privcon3 + privcon4
Risk  =~ risk1 + risk2 + risk3
Int   =~ intent1 + intent2

# Structural model
Int   ~ Trust + PrCon + Risk
Risk  ~ Trust + PrCon
Trust ~ PrCon
"

out <- csem(.data = LancelotMiltgenetal2016, .model = model_Med,
            .PLS_weight_scheme_inner = 'factorial',
            .tolerance = 1e-06
)
```
