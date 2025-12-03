# Do a redundancy analysis

**\[stable\]**

## Usage

``` r
doRedundancyAnalysis(.object = NULL)
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

## Value

A named numeric vector of correlations. If the weighting approach used
to obtain `.object` is not `"PLS-PM"` or non of the PLS outer modes was
mode B, the function silently returns `NA`.

## Details

Perform a redundancy analysis (RA) as proposed by Hair et al. (2016)
with reference to Chin (1998) .

RA is confined to PLS-PM, specifically PLS-PM with at least one
construct whose weights are obtained by mode B. In cSEM this is the case
if the construct is modeled as a composite or if argument `.PLS_modes`
was explicitly set to mode B for at least one construct. Hence RA is
only conducted if `.approach_weights = "PLS-PM"` and if at least one
construct's mode is mode B.

The principal idea of RA is to take two different measures of the same
construct and regress the scores obtained for each measure on each
other. If they are similar they are likely to measure the same "thing"
which is then taken as evidence that both measurement models actually
measure what they are supposed to measure (validity).

There are several issues with the terminology and the reasoning behind
this logic. RA is therefore only implemented since reviewers are likely
to demand its computation, however, its actual application for validity
assessment is discouraged.

Currently, the function is not applicable to models containing
second-order constructs.

## References

Chin WW (1998). “Modern Methods for Business Research.” In Marcoulides
GA (ed.), chapter The Partial Least Squares Approach to Structural
Equation Modeling, 295–358. Mahwah, NJ: Lawrence Erlbaum.  
  
Hair JF, Hult GTM, Ringle C, Sarstedt M (2016). *A Primer on Partial
Least Squares Structural Equation Modeling (PLS-SEM)*. Sage
publications.

## See also

[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
