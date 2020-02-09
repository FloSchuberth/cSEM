# cSEM 0.1.0:9000

- Fixed error in `predict()` when the dataset used to obtain `.object` contained 
  a character column. (#345)

- When calculating the HTMT via `assess()` the geometric mean of the average 
  monotrait−heteromethod correlation construct eta_i with the average 
  monotrait−heteromethod correlation of other constructs can be negative. 
  NaNs produced are produced in this case and the HTMT was not printed. 
  Added a warning and forced printing the NaNs as well. (#346)
  
- Fix bug in `testMICOM()`. Function produced an error if the data set provided
  contained an id-column even if the id-column was correctly supplied to 
  `csem()`. (#344, ##338)
  
- Fix bug in `doFloodlightAnalysis()`. There was an internal bug. Earlier versions
  returned the wrong direct effect. If you have used `doFloodlightAnalysis()`
  from cSEM v. 0.1.0 results are likely wrong.

- Add new function `getConstructScores()`. The function returns the standardized
  or unstandardized construct scores. Requires a `cSEMResults` object as input. (#340)

- Export plot method for `cSEMFloodlight` objects.

- Update documentation for `predict()`.

- Integrate and document `cSEMPredict` method for generic function `plot()`. Now 
  users may call `plot()` on an object created by `predict()`. (#337)

- Add the density of the residuals as plot to `plot.cSEMPredict()`. (#337)

- `assess()` now also computes and prints the total and indirect effects for each
  variable as they are often used for model assessment and may thus be considered 
  a quality criteria.
  In addition, the variance accounted for (VAF) is computed and printed as well. (#335)

- Change the name of the the quality criterion "effect size (f2)" from `esize` to `f2` 
  as this is more common. (#336)
  
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