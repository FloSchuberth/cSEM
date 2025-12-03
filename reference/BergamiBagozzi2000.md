# Data: BergamiBagozzi2000

A data frame containing 22 variables with 305 observations.

## Usage

``` r
BergamiBagozzi2000
```

## Format

An object of class `data.frame` with 305 rows and 22 columns.

## Source

Survey among South Korean employees conducted and reported by Bergami
and Bagozzi (2000) .

## Details

The dataset contains 22 variables and originates from a larger survey
among South Korean employees conducted and reported by Bergami and
Bagozzi (2000) . It is also used in Hwang and Takane (2004) and Henseler
(2021) for demonstration purposes, see the corresponding tutorial.

## References

Bergami M, Bagozzi RP (2000). “Self-categorization, affective commitment
and group self-esteem as distinct aspects of social identity in the
organization.” *British Journal of Social Psychology*, **39**(4),
555–577.
[doi:10.1348/014466600164633](https://doi.org/10.1348/014466600164633)
.  
  
Henseler J (2021). *Composite-Based Structural Equation Modeling:
Analyzing Latent and Emergent Variables*. Guilford Press, New York.  
  
Hwang H, Takane Y (2004). “Generalized Structured Component Analysis.”
*Psychometrika*, **69**(1), 81–99.

## Examples

``` r
#============================================================================
# Example is taken from Henseler (2021)
#============================================================================
model_Bergami_Bagozzi_Henseler="
# Measurement models
OrgPres =~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8 
OrgIden =~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
AffLove =~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
AffJoy  =~ orgcmt5 + orgcmt8
Gender  <~ gender

# Structural model 
OrgIden ~ OrgPres
AffLove ~ OrgPres + OrgIden + Gender 
AffJoy  ~ OrgPres + OrgIden + Gender 
"

out <- csem(.data = BergamiBagozzi2000, 
            .model = model_Bergami_Bagozzi_Henseler,
            .PLS_weight_scheme_inner = 'factorial',
            .tolerance = 1e-06
)

#============================================================================
# Example is taken from Hwang et al. (2004)
#============================================================================ 

model_Bergami_Bagozzi_Hwang="
# Measurement models
OrgPres =~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8 
OrgIden =~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
AffJoy =~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
AffLove  =~ orgcmt5 + orgcmt6 + orgcmt8

# Structural model 
OrgIden ~ OrgPres 
AffLove ~ OrgIden
AffJoy  ~ OrgIden"

out_Hwang <- csem(.data = BergamiBagozzi2000, 
                 .model = model_Bergami_Bagozzi_Hwang,
                 .approach_weights = "GSCA",
                 .disattenuate = FALSE,
                 .id = "gender",
                 .tolerance = 1e-06) 

```
