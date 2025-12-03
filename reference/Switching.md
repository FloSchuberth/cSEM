# Data: Switching

A data frame containing 26 variables with 767 observations.

## Usage

``` r
Switching
```

## Format

An object of class `data.frame` with 767 rows and 26 columns.

## Source

The dataset is provided by Jörg Henseler.

## Details

The data contains variables about the consumers’ intention to switch a
service provider. It is also used in Henseler (2021) for demonstration
purposes, see the corresponding tutorial.

## References

Henseler J (2021). *Composite-Based Structural Equation Modeling:
Analyzing Latent and Emergent Variables*. Guilford Press, New York.

## Examples

``` r
#============================================================================
# Example is taken from Henseler (2021)
#============================================================================
model_Int <-"
# Measurement models
INV =~ INV1 + INV2 + INV3 +INV4
SAT =~ SAT1 + SAT2 + SAT3
INT =~ INT1 + INT2

# Structural model containing an interaction term.
INT ~ INV + SAT + INV.SAT
"

out <- csem(.data = Switching, .model = model_Int,
            .PLS_weight_scheme_inner = 'factorial',
            .tolerance = 1e-06)
```
