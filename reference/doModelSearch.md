# Automated model specification search

**\[experimental\]**

## Usage

``` r
doModelSearch(
 .object          = NULL, 
 .pop_size        = 20,
 .n_generations   = 20,
 .prob_mutation   = 0.5,
 .prob_crossover  = 0.8,
 .fbar            = -100000,
 .ms_criterion    = c('bic','aic','hq'),
 .seed            = NULL,
 .only_structural = FALSE
)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .pop_size:

  Integer. Population size used in the genetic algorithm. Defaults to
  `20`.

- .n_generations:

  Integer. Number of generations used in the genetic algorithm. Defaults
  to `20`.

- .prob_mutation:

  Scalar. Mutation probability used in the genetic algorithm. Defaults
  to `0.5`.

- .prob_crossover:

  Scalar. Crossover probability used in the genetic algorithm. Defaults
  to `0.8`.

- .fbar:

  Integer. Low fitness value that is used to penalize inadmissible
  models. Defaults to -100000.

- .ms_criterion:

  Character string. A single character string naming the model selection
  criterion to compute. One of `"bic"`, `"aic"`, or `"hq"`. Defaults to
  `"bic"`

- .seed:

  Integer or `NULL`. The random seed to use. Defaults to `NULL` in which
  case an arbitrary seed is chosen. Note that the scope of the seed is
  limited to the body of the function it is used in. Hence, the global
  seed will not be altered!

- .only_structural:

  Should the the log-likelihood be based on the structural model?
  Ignored if `.by_equation == TRUE`. Defaults to `TRUE`.

## Value

A list of class `cSEMModelSearch`. The list contains the follwing
elements:

- `$Information`:

  A list of input parameters to `doModelSearch`

- `$Results`:

  A list with elements, `original_fitness`, `best_fitness`, and
  `best_matrix`

- `$Inputcsem`:

  A list with elements, `data`, and `model`, where `data` comprises the
  dataset used to estimate the original model and `model` is a list of
  [cSEMModel](https://floschuberth.github.io/cSEM/reference/csem_model.md)
  lists.

## Details

Performs an automated model specification search using a genetic
algorithm (Holland 1992) in combination with partial least squares path
modeling (PLS-PM). Specifically, the function implements the approach
presented in Trinchera et al. (2026) . The model search is limited to
the structural model. Thereby it is ensured that the resulting
structural model is recursive, i.e., feedback loops and non-recursive
structures are prohibited. Moreover, the set of exogenous constructs is
taken from the original model provided to `doModelSearch()` and remains
exogenous during the search. If a model violates the restrictions or
does not pass
[`verify()`](https://floschuberth.github.io/cSEM/reference/verify.md)
during the model search, a penalty fitness value `.fbar` is assigned to
the model. Therefore, ensure that `.fbar` is sufficiently small. The
fitness of the model is determined by the criterion provided to
`.ms_criterion`.

To estimate the found model `Inputcsem` results can be used, see the
example.

## References

Holland JH (1992). “Genetic Algorithms.” *Scientific American*,
**267**(1), 66–72.
[doi:10.1038/scientificamerican0792-66](https://doi.org/10.1038/scientificamerican0792-66)
.  
  
Trinchera L, Pietropolli G, Castelli M, Schuberth F (2026). “Automated
specification search for composite-based structural equation modeling: a
genetic approach.” *Computational Statistics & Data Analysis*.

## See also

[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
[calculateModelSelectionCriteria](https://floschuberth.github.io/cSEM/reference/calculateModelSelectionCriteria.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Perform model search for a linear model without second-order constructs 
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

 
 outSearch <- doModelSearch(.object = out,
               .pop_size = 20,
               .n_generations = 20,
               .prob_mutation = 0.5,
               .prob_crossover = 0.8,
               .fbar = -100000,
               .ms_criterion = 'bic',
               .seed = 1234) 
 outSearch
 
#  Estimate the found model
 outNew <- csem(.data = outSearch$Inputcsem$data,
                .model = outSearch$Inputcsem$model[[1]])

 # summarize(outNew) 
} # }
```
