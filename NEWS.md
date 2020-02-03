# cSEM 0.1.0:9000

- Fix bug related to dotdotdot arguments incorrectly passed to functions supplied
  to `.user_funs` when resampling. Add additional example to `assess()` illustrating
  the use of the `.user_funs` arguments when given multiple functions. (#334) 

- Add both versions of the RMS_theta to `assess()`. `"RMS_theta"` is the RMS_theta
  based on WSW'. `"RMS_theta_mi"` uses the model-implied construct correlation matrix.
  The argument `.model_implied` is thus no longer available to assess's `...`
  arguments and therefore removed.
  
- Allow users to specify a lavaan model without a structural model. Now, users
  can specify a model with several measurement equations (via `<~` or `=~`)
  but no strucutral equations. Instead the correlations between all! constructs
  must be given. Failing to do so causes an error.

- Add CITATION file (#331)

- Add informative error message if `.data` contains missing values.

- Add the Chi_square statistic and the Chi_square statistic divided by its
  degrees of freedom to the list of fit indices
  
- Remove argument `.only_common_factors` for postestimation function `predict()`.
  Now `predict()` retruns predictions for composite models as well.
  This will break existing code that uses `predict(..., .only_common_factors = ...)`.
  You will get an `unused argument (.only_common_factors = FALSE)` error. 
  Simply remove the argument to fix it. (#330)

- Fix bug in `calculateEffectSize()`/`assess()` when one of the equations
  of the structural model has only one explanatory variable. 

- Remove warning produced when printing a `cSEMAssess` object based on a
  model containing only constructs modeled as composites.

- Remove warning produced by `calculateHTMT()` when the estimated model contains
  < 2 common factors. (#325)

- Update documentation for `assess()` and `calcaulateRhoC()`

- Update Vignettes `using-asses` and `csem-notation`
  
  
# Initial release: cSEM 0.1.0 (07.01.2020)