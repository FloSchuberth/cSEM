# Data: Russett

A data frame containing 10 variables with 47 observations.

## Usage

``` r
Russett
```

## Format

A data frame containing the following variables for 47 countries:

- `gini`:

  The Gini index of concentration

- `farm`:

  The percentage of landholders who collectively occupy one-half of all
  the agricultural land (starting with the farmers with the smallest
  plots of land and working toward the largest)

- `rent`:

  The percentage of the total number of farms that rent all their land.
  Transformation: ln (x + 1)

- `gnpr`:

  The 1955 gross national product per capita in U.S. dollars.
  Transformation: ln (x)

- `labo`:

  The percentage of the labor force employed in agriculture.
  Transformation: ln (x)

- `inst`:

  Instability of personnel based on the term of office of the chief
  executive. Transformation: exp (x - 16.3)

- `ecks`:

  The total number of politically motivated violent incidents, from
  plots to protracted guerrilla warfare. Transformation: ln (x + 1)

- `deat`:

  The number of people killed as a result of internal group violence per
  1,000,000 people. Transformation: ln (x + 1)

- `stab`:

  One if the country has a stable democracy, and zero otherwise

- `dict`:

  One if the country experiences a dictatorship, and zero otherwise

## Source

From: Henseler (2021)

## Details

The dataset was initially compiled by Russett (1964) , discussed and
reprinted by Gifi (1990) , and partially transformed by Tenenhaus and
Tenenhaus (2011) . It is also used in Henseler (2021) for demonstration
purposes.

## References

Gifi A (1990). *Nonlinear multivariate analysis*. Wiley.  
  
Henseler J (2021). *Composite-Based Structural Equation Modeling:
Analyzing Latent and Emergent Variables*. Guilford Press, New York.  
  
Russett BM (1964). “Inequality and Instability: The Relation of Land
Tenure to Politics.” *World Politics*, **16**(3), 442–454.
[doi:10.2307/2009581](https://doi.org/10.2307/2009581) .  
  
Tenenhaus A, Tenenhaus M (2011). “Regularized generalized canonical
correlation analysis.” *Psychometrika*, **76**(2), 257–284.

## Examples

``` r
#============================================================================
# Example is taken from Henseler (2020)
#============================================================================
model_Russett="
# Composite model
AgrIneq <~ gini + farm + rent
IndDev  <~ gnpr + labo
PolInst <~ inst + ecks + deat + stab + dict

# Structural model
PolInst ~ AgrIneq + IndDev
"

out <- csem(.data = Russett, .model = model_Russett,
            .PLS_weight_scheme_inner = 'factorial',
            .tolerance = 1e-06
)
```
