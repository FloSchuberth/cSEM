# cSEM 0.1.0:9000

- Fix bug in `calculateEffectSize()`/`assess()` when one of the equations
  of the structural model has only one explanatory variable. 

- Remove warning produced when printing a `cSEMAssess` object based on a
  model containing only constructs modeled as composites.

- Remove warning produced by `calculateHTMT()` when the estimated model contains
  < 2 common factors. 

- Update documentation for `assess()` and `calcaulateRhoC()`

- Update Vignettes `using-asses` and `csem-notation`
  
  
# Initial release: cSEM 0.1.0 (07.01.2020)