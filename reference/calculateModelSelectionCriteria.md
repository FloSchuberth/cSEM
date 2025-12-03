# Model selection criteria

Calculate several information or model selection criteria (MSC) such as
the Akaike information criterion (AIC), the Bayesian information
criterion (BIC) or the Hannan-Quinn criterion (HQ).

## Usage

``` r
calculateModelSelectionCriteria(
  .object          = NULL,
  .ms_criterion    = c("all", "aic", "aicc", "aicu", "bic", "fpe", "gm", "hq",
                       "hqc", "mallows_cp"),
  .by_equation     = TRUE, 
  .only_structural = TRUE 
  )
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .ms_criterion:

  Character string. Either a single character string or a vector of
  character strings naming the model selection criterion to compute.
  Defaults to `"all"`.

- .by_equation:

  Should the criteria be computed for each structural model equation
  separately? Defaults to `TRUE`.

- .only_structural:

  Should the the log-likelihood be based on the structural model?
  Ignored if `.by_equation == TRUE`. Defaults to `TRUE`.

## Value

If `.by_equation == TRUE` a named list of model selection criteria.

## Details

By default, all criteria are calculated (`.ms_criterion == "all"`). To
compute only a subset of the criteria a vector of criteria may be given.

If `.by_equation == TRUE` (the default), the criteria are computed for
each structural equation of the model separately, as suggested by Sharma
et al. (2019) in the context of PLS. The relevant formula can be found
in Table B1 of the appendix of Sharma et al. (2019) .

If `.by_equation == FALSE` the AIC, the BIC and the HQ for whole model
are calculated. All other criteria are currently ignored in this case!
The relevant formula are (see, e.g., (Akaike 1974) , Schwarz (1978) ,
Hannan and Quinn (1979) ):

\$\$AIC = - 2\*log(L) + 2\*k\$\$ \$\$BIC = - 2\*log(L) + k\*ln(n)\$\$
\$\$HQ = - 2\*log(L) + 2\*k\*ln(ln(n))\$\$

where log(L) is the log likelihood function of the multivariate normal
distribution of the observable variables, k the (total) number of
estimated parameters, and n the sample size.

If `.only_structural == TRUE`, log(L) is based on the structural model
only. The argument is ignored if `.by_equation == TRUE`.

## References

Akaike H (1974). “A New Look at the Statistical Model Identification.”
*IEEE Transactions on Automatic Control*, **19**(6), 716–723.  
  
Hannan EJ, Quinn BG (1979). “The Determination of the order of an
autoregression.” *Journal of the Royal Statistical Society: Series B
(Methodological)*, **41**(2), 190–195.  
  
Schwarz G (1978). “Estimating the Dimension of a Model.” *The Annals of
Statistics*, **6**(2), 461–464.
[doi:10.1214/aos/1176344136](https://doi.org/10.1214/aos/1176344136) .  
  
Sharma P, Sarstedt M, Shmueli G, Kim KH, Thiele KO (2019). “PLS-Based
Model Selection: The Role of Alternative Explanations in Information
Systems Research.” *Journal of the Association for Information Systems*,
**20**(4).

## See also

[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
