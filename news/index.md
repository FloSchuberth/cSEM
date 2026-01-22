# Changelog

## \<Development version: cSEM 0.6.1.9999

## cSEM 0.6.2

- Implement AGAS-PLS. Thanks to Gloria Pietropolli

- Update R dependency (R\>= 4.1.0)

## cSEM 0.6.1

CRAN release: 2025-05-16

- fix issue in the
  [`fit()`](https://floschuberth.github.io/cSEM/reference/fit.md)
  function to avoid partial matching. Thanks to Kss2k for this
  contribution
  ([\#558](https://github.com/FloSchuberth/cSEM/issues/558))

- Adjust description file, i.e., NeedsCompilation: no

## cSEM 0.6.0 (2025-02-24)

CRAN release: 2025-02-25

- Implement [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  function which visualizes cSEM models. Thanks to Nguyen Huu Phuc for
  this contribution
  ([\#430](https://github.com/FloSchuberth/cSEM/issues/430)).

- Bug fix in the predict function when the model contains only a single
  categorical indicator.

- Bug fix in case of a second-order composite that is formed by one
  common factor.

- Implement `calculateRelativeGoF`. Now the relative GoF can be
  obtained.

- Bug fix in `calculateGoF`. Now single-indicator constructs are
  excluded. Thanks to Mehmet Mehmetoglu and Sergio Venturini.

- Bug fix in
  [`calculateEffects()`](https://floschuberth.github.io/cSEM/reference/calculateEffects.md).
  Now it is distinguished between recursive and non-recursive models.
  For recursive models rounding is no longer necessary.

- Bug fix in
  [`calculateReliabilities()`](https://floschuberth.github.io/cSEM/reference/calculateReliabilities.md).
  Now correction for attenuation is correctly done for PLS-PM Mode B.

## cSEM 0.5.0 (2022-05-09)

CRAN release: 2022-11-24

- Bug fix in
  [`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md).
  Now the saturated argument is passed to the discrepancy/fit measures.

- Bug fix in `.resampleData()` when crossvalidation is used. Empty
  datasets are not possible anymore.

- Implemented several other prediction metrics.

- Bug fix: Revision of the predict metrics in the `predict` function.

- Update the .eval_plan argument since the multiprocess argument of the
  future package is deprecated. Now multisession or multicore need to be
  used. Note multicore does not work on Windows machines.

- Bug fix: Calculation of the R2 and adjR2 in the print function of the
  assess function

- Revise the description of the two-stage approach in the csem help file
  ([\#418](https://github.com/FloSchuberth/cSEM/issues/418))

- Bug fix: fix print method for
  [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
  when disattenuate is set to TRUE internally. Now disattenuate as
  treated in csem is reported and not the value provided by the user.
  ([\#419](https://github.com/FloSchuberth/cSEM/issues/419))

- Use singular value decomposition in GSCAm to deal with large datasets
  ([\#444](https://github.com/FloSchuberth/cSEM/issues/444))

- Bug fix: GSCAm (i.e., `.approach_weights = "GSCA"` with constructs
  modeled as common factors) no longer fails when a single indicator
  construct is supplied
  ([\#441](https://github.com/FloSchuberth/cSEM/issues/441))

- The default value for argument `.r` (the number of repetitions) of
  [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md)
  was changed from 10 to 1 since more than one repetition is hardly ever
  necessary.

- [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md)
  is now able to predict categorical indicators (a procedure known as
  OrdPLScPredict).
  [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md)
  therefore gains a number of new arguments, namely:
  `.approach_score_target`, `.sim_points`, `.treat_as_continuous`, and
  `.approach_score_benchmark`.

- Removed argument `.verbose` from
  [`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md)
  as it did not have any effect
  ([\#445](https://github.com/FloSchuberth/cSEM/issues/445)).

- Bug fix: GSCAm (i.e., `.approach_weights = "GSCA"` with constructs
  modeled as common factors) no longer fails when a single indicator
  construct is supplied
  ([\#441](https://github.com/FloSchuberth/cSEM/issues/441))

- Bug fix:
  [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md)
  no longer fails when LOOCV is used
  ([\#337](https://github.com/FloSchuberth/cSEM/issues/337))

- Bug fix: fix print method for
  [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
  when resampling with constant values (weights or loadings) is
  conducted. The standard error, t-value, p-value and CI are properly
  set to `NA` now.
  ([\#433](https://github.com/FloSchuberth/cSEM/issues/433))

- `postestimate_test_CVPAT()`: Perform a Cross-Validated Predictive
  Ability Test (CVPAT) to compare the predictive performance of two
  models ([\#455](https://github.com/FloSchuberth/cSEM/issues/455))

- [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md)
  is now able to perform predictions either based on the earliest
  antecedents, i.e., the values of the indicators associated to
  exogenous constructs or based on the direct antecedents, i.e., based
  on the values or predictions associated to the direct antecedents
  (3485)

- [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md)
  a variety of prediction metrics are added

## cSEM 0.4.0 (2021-04-20)

CRAN release: 2021-04-19

#### Major changes

- New function
  [`exportToExcel()`](https://floschuberth.github.io/cSEM/reference/exportToExcel.md).
  The function conveniently exports the results from
  [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md),
  [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md),
  [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
  and
  [`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md)
  to an .xlsx file.

#### Bug fixes

- Critical bug fix: `calculateVifModeB()` did not calculate the VIFs for
  modeB constructs correctly because of a bug in the calculation of the
  R^2. PLEASE REVIEW YOUR CALCULATIONS in cSEM version \< 0.3.1:9000!
  (thanks to [@Benjamin](https://github.com/Benjamin) Liengaard for
  pointing it out).

- Bug fix:
  [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md)
  no longer silently returns empty predictions when `.test_data` does
  not contain rownames.

- Bug fix: calculation of the MSE in `modelSelectionCriteria()` resulted
  in a vector of incorrect length. In some cases this affected the
  computation of “GM” and “Mallows_cp”.

## cSEM 0.3.1 (2021-02-14)

CRAN release: 2021-02-14

- Bug fix:
  [`summarize()`](https://floschuberth.github.io/cSEM/reference/summarize.md)
  no longer fails when `.object` is a of class `cSEMResults_2ndorder`
  and contains no indirect effects.

- Add argument `type_htmt` to
  [`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md).
  `type_htmt = "htmt2"` calculates a consistent estimator for congeneric
  measurement models.

## cSEM 0.3.0 (2020-12-10)

CRAN release: 2020-10-12

- Add lifecylce badges to postestimation
  functions.([\#376](https://github.com/FloSchuberth/cSEM/issues/376))

- Some arguments accepted by
  [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)’s
  `...` argument had not been documented properly. This has been fixed.
  See `args_assess_dotdotdot` for a complete list of available
  arguments.

- [`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md)
  now allows users to chose the type of confidence interval to use when
  computing the critical (1-alpha)% quantile of the HTMT values
  ([\#379](https://github.com/FloSchuberth/cSEM/issues/379))

- [`testMGD()`](https://floschuberth.github.io/cSEM/reference/testMGD.md)
  gains a new `.output_type` argument. By default
  (`.output_type = "structured"`), the standard output is returned. If
  `.output_type = "structured"`, however, a tibble (data frame)
  summarizing the test decisions in a user-friendly way is returned.
  ([\#398](https://github.com/FloSchuberth/cSEM/issues/398))

- Remove warning from
  [`fit()`](https://floschuberth.github.io/cSEM/reference/fit.md) when
  polychoric or polyserial indicator correlation is used during
  estimation. ([\#413](https://github.com/FloSchuberth/cSEM/issues/413))

- [`print.cSEMAssess()`](https://floschuberth.github.io/cSEM/reference/print.cSEMAssess.md)
  no longer prints zero for VIF values of constructs that are not part
  of a particular structural equation.

- [`print.cSEMAssess()`](https://floschuberth.github.io/cSEM/reference/print.cSEMAssess.md)
  now prints the results of
  [`calculateVIFModeB()`](https://floschuberth.github.io/cSEM/reference/calculateVIFModeB.md).
  This had been missing in previous releases.
  ([\#384](https://github.com/FloSchuberth/cSEM/issues/384))

- Breaking:
  [`calculateVIFModeB()`](https://floschuberth.github.io/cSEM/reference/calculateVIFModeB.md)
  now returns a matrix with the dependent construct in the rows and the
  VIFs for the coresponding weights in the columns. Previously, the
  output was a list.

- Add model selection criteria. See the
  [`calculateModelSelectionCriteria()`](https://floschuberth.github.io/cSEM/reference/calculateModelSelectionCriteria.md)
  function for details. As usual, all criteria are available via
  [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md).
  ([\#412](https://github.com/FloSchuberth/cSEM/issues/412))

- Combine functions for surface, floodlight and simple effects analysis
  in the
  [`doNonlinearEffectsAnalysis()`](https://floschuberth.github.io/cSEM/reference/doNonlinearEffectsAnalysis.md)
  function; Breaking: functions `doFloodlightAnalysis()` and
  `doSurfaceAnalysis()` have been removed!

- Progress bars are now supported for every function that does
  resampling. Progress bars are fully customizable via the `progressr`
  framework created by

  1.  Note: to suppress the progress bar use
      `progressr::handlers("void")` and then run your csem commands.
      ([\#359](https://github.com/FloSchuberth/cSEM/issues/359))

- Fix bug in the computation of the Bc and Bca interval. Computation
  failed for models that had no indirect effects.

- List element “reliability” of
  [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
  is changed to “Reliability” to be consistent with the naming scheme of
  the other list elements.

- [`infer()`](https://floschuberth.github.io/cSEM/reference/infer.md)
  automatically computes bootstrap resamples now by default if `.object`
  does not have class `cSEMResults_resampled` already.
  ([\#389](https://github.com/FloSchuberth/cSEM/issues/389))

- Remove `.alpha` argument from
  [`testMICOM()`](https://floschuberth.github.io/cSEM/reference/testMICOM.md).
  The argument is no longer required as decisions are made via (possibly
  adjusted) p-values.
  ([\#393](https://github.com/FloSchuberth/cSEM/issues/393))

- Add checks to plot methods for
  [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md),
  `doFloodlightAnalysis`, and, `doFloodlightAnalysis`.

- Several documentation updates and typo corrections.

- The Fornell-Larcker criterion is now computed by its own function
  [`calculateFLCriterion()`](https://floschuberth.github.io/cSEM/reference/calculateFLCriterion.md).
  Previously, it was only available via
  [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md).
  ([\#387](https://github.com/FloSchuberth/cSEM/issues/387))

- Implement importance-performance matrix analysis via
  [`doIPMA()`](https://floschuberth.github.io/cSEM/reference/doIPMA.md).
  A corresponding plot method is also available.

## cSEM 0.2.0 (30.03.2020)

CRAN release: 2020-03-29

### Major changes

- [`testMICOM()`](https://floschuberth.github.io/cSEM/reference/testMICOM.md)
  gains the `.approach_p_adjust` argument. The argument takes a single
  character string or a vector of character strings naming the p-value
  adjustment for multiple comparisons.
  ([\#138](https://github.com/FloSchuberth/cSEM/issues/138))

- Review
  [`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md).
  1.) Add inference; 2) fix wrong handling of single-indicator
  constructs
  ([\#351](https://github.com/FloSchuberth/cSEM/issues/351)); 3) Remove
  warning produced by
  [`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md)
  when the estimated model contains less than 2 common factors.
  ([\#325](https://github.com/FloSchuberth/cSEM/issues/325))

- Breaking: Rename argument in `doFloodlightAnalysis()`.
  ([\#343](https://github.com/FloSchuberth/cSEM/issues/343))

- New function `doSurfaceAnalysis()`. See
  `?doSurfaceAnalysis()`([\#349](https://github.com/FloSchuberth/cSEM/issues/349))

- Implement degrees of freedom calculation for second-order constructs.

- Add new function
  [`getConstructScores()`](https://floschuberth.github.io/cSEM/reference/getConstructScores.md).
  The function returns the standardized or unstandardized construct
  scores. Requires a `cSEMResults` object as input.
  ([\#340](https://github.com/FloSchuberth/cSEM/issues/340))

- Fix bug in `doFloodlightAnalysis()`. There was an internal bug.
  Earlier versions returned the wrong direct effect. If you have used
  `doFloodlightAnalysis()` from cSEM v. 0.1.0 results are likely wrong.

- Export plot method for `cSEMFloodlight` objects.

- Allow users to specify a lavaan model without a structural model. Now,
  users can specify a model with several measurement equations (via `<~`
  or `=~`) but no structural equations. Instead the correlations between
  all! constructs must be given. Failing to do so causes an error.

#### New example data

- Add indicator correlation matrix for a modified version of
  Summers (1965) model. See
  [`?Sigma_Summers_composites`](https://floschuberth.github.io/cSEM/reference/Sigma_Summers_composites.md)

- Add example data sets used in Henseler (2020). See
  [`?BergamiBagozzi2000`](https://floschuberth.github.io/cSEM/reference/BergamiBagozzi2000.md),
  [`?ITFlex`](https://floschuberth.github.io/cSEM/reference/ITFlex.md),
  [`?LancelotMiltgenetal2016`](https://floschuberth.github.io/cSEM/reference/LancelotMiltgenetal2016.md),
  [`?Russett`](https://floschuberth.github.io/cSEM/reference/Russett.md),
  [`?Switching`](https://floschuberth.github.io/cSEM/reference/Switching.md),
  and
  [`?Yooetal2000`](https://floschuberth.github.io/cSEM/reference/Yooetal2000.md).

#### assess()

- Update documentation and vignettes

- The following functions called by
  [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
  are now exported and support all of cSEM’s native classes
  ([\#357](https://github.com/FloSchuberth/cSEM/issues/357),
  [\#369](https://github.com/FloSchuberth/cSEM/issues/369)):

  - [`calculateAVE()`](https://floschuberth.github.io/cSEM/reference/calculateAVE.md)
  - [`calculateDf()`](https://floschuberth.github.io/cSEM/reference/calculateDf.md)
  - [`calculateGoF()`](https://floschuberth.github.io/cSEM/reference/calculateGoF.md)
  - [`calculateHTMT()`](https://floschuberth.github.io/cSEM/reference/calculateHTMT.md)
    (does not support models containing second-order constructs)
  - [`calculateRhoT()`](https://floschuberth.github.io/cSEM/reference/reliability.md)
  - [`calculateRhoT()`](https://floschuberth.github.io/cSEM/reference/reliability.md)
  - [`calculatef2()`](https://floschuberth.github.io/cSEM/reference/calculatef2.md)
  - [`calculateDML()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md)
  - [`calculateDG()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md)
  - [`calculateDL()`](https://floschuberth.github.io/cSEM/reference/distance_measures.md)
  - [`calculateChiSquare()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md)
  - [`calculateChiSquareDf()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md)
  - [`calculateGFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md)
  - [`calculateNFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md)
  - [`calculateNNFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md)
  - [`calculateIFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md)
  - [`calculateCFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md)
  - [`calculateSRMR()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md)
  - [`calculateRMSEA()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md)
  - [`calculateRMSTheta()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md)
  - [`calculateVIFModeB()`](https://floschuberth.github.io/cSEM/reference/calculateVIFModeB.md)

- [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
  now supports all of cSEM’s native classes.
  ([\#323](https://github.com/FloSchuberth/cSEM/issues/323))

- [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
  now also computes and prints the total and indirect effects for each
  variable as they are often used for model assessment and may thus be
  considered a quality criteria. In addition, the variance accounted for
  (VAF) is computed and printed as well.
  ([\#335](https://github.com/FloSchuberth/cSEM/issues/335))

- Breaking: change the name of the the quality criterion “effect size
  (f2)” from `esize` to `f2` and the corresponding function from
  `calculateEffectSize()` to
  [`calculatef2()`](https://floschuberth.github.io/cSEM/reference/calculatef2.md)as
  this is more common.
  ([\#336](https://github.com/FloSchuberth/cSEM/issues/336))

- Add the Chi_square statistic and the Chi_square statistic divided by
  its degrees of freedom to the list of fit indices. See:
  `?calculateChiSquare()` and `?calculateChiSquareDf()`

- Fix bug in
  [`calculatef2()`](https://floschuberth.github.io/cSEM/reference/calculatef2.md)/[`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
  when one of the equations of the structural model has only one
  explanatory variable.

- Fix bug related to dotdotdot arguments incorrectly passed to functions
  supplied to `.user_funs` when resampling. Add additional example to
  [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
  illustrating the use of the `.user_funs` arguments when given multiple
  functions. ([\#334](https://github.com/FloSchuberth/cSEM/issues/334))

- Remove warning produced when printing a `cSEMAssess` object based on a
  model containing only constructs modeled as composites.

#### predict()

- Update documentation for
  [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md).

- Integrate and document `cSEMPredict` method for generic function
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html). Now users
  may call [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on
  an object created by
  [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md).
  ([\#337](https://github.com/FloSchuberth/cSEM/issues/337))

- Add the density of the residuals as plot to
  [`plot.cSEMPredict()`](https://floschuberth.github.io/cSEM/reference/plot.cSEMPredict.md).
  ([\#337](https://github.com/FloSchuberth/cSEM/issues/337))

- Remove argument `.only_common_factors` for postestimation function
  [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md).
  Now
  [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md)
  returns predictions for composite models as well. This will break
  existing code that uses `predict(..., .only_common_factors = ...)`.
  You will get an `unused argument (.only_common_factors = FALSE)`
  error. Simply remove the argument to fix it.
  ([\#330](https://github.com/FloSchuberth/cSEM/issues/330))

- Fixed error in
  [`predict()`](https://floschuberth.github.io/cSEM/reference/predict.md)
  when the dataset used to obtain `.object` contained a character
  column. ([\#345](https://github.com/FloSchuberth/cSEM/issues/345))

#### Experimental features

- Add `.fit_measures` argument to
  [`testOMF()`](https://floschuberth.github.io/cSEM/reference/testOMF.md).
  Now other fit measures such as the RMSEA or the GFI can be used as the
  test statistic. This is a rather experimental feature and may be
  removed in future versions.

### Minor changes and bug fixes

- Using `.approach_weights = "GSCA"` for models containing nonlinear
  terms gives a more meaningful error message.
  ([\#342](https://github.com/FloSchuberth/cSEM/issues/342))

- [`print.cSEMTestMICOM()`](https://floschuberth.github.io/cSEM/reference/print.cSEMTestMICOM.md)
  no longer prints the decision but additional bootstrap information.
  (in parts: [\#339](https://github.com/FloSchuberth/cSEM/issues/339))

- If the weighting scheme is `"PLS-PM"` and `.disattenuate = TRUE`,
  dissatenuation is longer applied to constructs using modes other than
  “modeA”” or “modeB”.
  ([\#352](https://github.com/FloSchuberth/cSEM/issues/352))

- Model-implied indicator correlation matrix for non-recursive models
  should now be calculated correctly.
  ([\#264](https://github.com/FloSchuberth/cSEM/issues/264))

- [`calculatef2()`](https://floschuberth.github.io/cSEM/reference/calculatef2.md)
  gives an error when the path model estimator is not “OLS”.
  ([\#360](https://github.com/FloSchuberth/cSEM/issues/360),
  [\#370](https://github.com/FloSchuberth/cSEM/issues/370))

- Add `.type` argument to
  [`calculateGFI()`](https://floschuberth.github.io/cSEM/reference/fit_measures.md).
  Now GFI based on the ML and ULS fitting function can be computed.
  ([\#371](https://github.com/FloSchuberth/cSEM/issues/371))

- [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md)
  gives a meaningful error when the structural model contains only
  second-order constructs
  ([\#366](https://github.com/FloSchuberth/cSEM/issues/366))

- Fix bug in
  [`testMICOM()`](https://floschuberth.github.io/cSEM/reference/testMICOM.md).
  Function produced an error if the data set provided contained more
  columns than indicators used in the model used for
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).
  ([\#355](https://github.com/FloSchuberth/cSEM/issues/355))

- Fix bug in
  [`testMICOM()`](https://floschuberth.github.io/cSEM/reference/testMICOM.md).
  Function produced an error if the data set provided contained an
  id-column even if the id-column was correctly supplied to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).
  ([\#344](https://github.com/FloSchuberth/cSEM/issues/344),
  [\#338](https://github.com/FloSchuberth/cSEM/issues/338))

- When calculating the HTMT via
  [`assess()`](https://floschuberth.github.io/cSEM/reference/assess.md)
  the geometric mean of the average monotrait−heteromethod correlation
  construct eta_i with the average monotrait−heteromethod correlation of
  other constructs can be negative. NaNs produced are produced in this
  case and the HTMT was not printed. Added a warning and forced printing
  the NaNs as well.
  ([\#346](https://github.com/FloSchuberth/cSEM/issues/346))

- Add CITATION file
  ([\#331](https://github.com/FloSchuberth/cSEM/issues/331))

- Add informative error message if `.data` contains missing values.

- Update vignettes `csem-notation`

## Initial release: cSEM 0.1.0 (07.01.2020)

CRAN release: 2020-01-13
