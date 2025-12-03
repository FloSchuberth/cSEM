# Predict indicator scores

**\[maturing\]**

## Usage

``` r
predict(
 .object                   = NULL,
 .benchmark                = c("lm", "unit", "PLS-PM", "GSCA", "PCA", "MAXVAR", "NA"),
 .approach_predict         = c("earliest", "direct"),
 .cv_folds                 = 10,
 .handle_inadmissibles     = c("stop", "ignore", "set_NA"),
 .r                        = 1,
 .test_data                = NULL,
 .approach_score_target    = c("mean", "median", "mode"),
 .sim_points               = 100,
 .disattenuate             = TRUE,
 .treat_as_continuous      = TRUE,
 .approach_score_benchmark = c("mean", "median", "mode", "round"),
 .seed                     = NULL
 )
```

## Arguments

- .object:

  An R object of class
  [cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md)
  resulting from a call to
  [`csem()`](https://floschuberth.github.io/cSEM/reference/csem.md).

- .benchmark:

  Character string. The procedure to obtain benchmark predictions. One
  of "*lm*", "*unit*", "*PLS-PM*", "*GSCA*", "*PCA*", "*MAXVAR*", or
  "*NA*". Default to "*lm*".

- .approach_predict:

  Character string. Which approach should be used to perform
  predictions? One of "*earliest*" and "*direct*". If "*earliest*"
  predictions for indicators associated to endogenous constructs are
  performed using only indicators associated to exogenous constructs. If
  "*direct*", predictions for indicators associated to endogenous
  constructs are based on indicators associated to their direct
  antecedents. Defaults to "*earliest*".

- .cv_folds:

  Integer. The number of cross-validation folds to use. Setting
  `.cv_folds` to `N` (the number of observations) produces leave-one-out
  cross-validation samples. Defaults to `10`.

- .handle_inadmissibles:

  Character string. How should inadmissible results be treated? One of
  "*stop*", "*ignore*", or "*set_NA*". If "*stop*", `predict()` will
  stop immediately if estimation yields an inadmissible result. For
  "*ignore*" all results are returned even if all or some of the
  estimates yielded inadmissible results. For "*set_NA*" predictions
  based on inadmissible parameter estimates are set to `NA`. Defaults to
  "*stop*"

- .r:

  Integer. The number of repetitions to use. Defaults to `1`.

- .test_data:

  A matrix of test data with the same column names as the training data.

- .approach_score_target:

  Character string. How should the aggregation of the estimates of the
  truncated normal distribution for the predictions using OrdPLS/OrdPLSc
  be done? One of "*mean*", "*median*" or "*mode*". If "*mean*", the
  mean of the estimated endogenous indicators is calculated. If
  "*median*", the mean of the estimated endogenous indicators is
  calculated. If "*mode*", the maximum empirical density on the
  intervals defined by the thresholds is used. Defaults to "*mean*".

- .sim_points:

  Integer. How many samples from the truncated normal distribution
  should be simulated to estimate the exogenous construct scores?
  Defaults to "*100*".

- .disattenuate:

  Logical. Should the benchmark predictions be based on disattenuated
  parameter estimates? Defaults to `TRUE`.

- .treat_as_continuous:

  Logical. Should the indicators for the benchmark predictions be
  treated as continuous? If `TRUE` all indicators are treated as
  continuous and PLS-PM/PLSc is applied. If `FALSE` OrdPLS/OrdPLSc is
  applied. Defaults to `TRUE`.

- .approach_score_benchmark:

  Character string. How should the aggregation of the estimates of the
  truncated normal distribution be done for the benchmark predictions?
  Ignored if not OrdPLS or OrdPLSc is used to obtain benchmark
  predictions. One of "*mean*", "*median*", "*mode*" or "*round*". If
  "*round*", the benchmark predictions are obtained using the
  traditional prediction algorithm for PLS-PM which are rounded for
  categorical indicators. If "*mean*", the mean of the estimated
  endogenous indicators is calculated. If "*median*", the mean of the
  estimated endogenous indicators is calculated. If "*mode*", the
  maximum empirical density on the intervals defined by the thresholds
  is used. If `.treat_as_continuous = TRUE` or if all indicators are on
  a continuous scale, `.approach_score_benchmark` is ignored. Defaults
  to "*round*".

- .seed:

  Integer or `NULL`. The random seed to use. Defaults to `NULL` in which
  case an arbitrary seed is chosen. Note that the scope of the seed is
  limited to the body of the function it is used in. Hence, the global
  seed will not be altered!

## Value

An object of class `cSEMPredict` with print and plot methods.
Technically, `cSEMPredict` is a named list containing the following list
elements:

- `$Actual`:

  A matrix of the actual values/indicator scores of the endogenous
  constructs.

- `$Prediction_target`:

  A list containing matrices of the predicted indicator scores of the
  endogenous constructs based on the target model for each repetition
  .r. Target refers to procedure used to estimate the parameters in
  `.object`.

- `$Residuals_target`:

  A list of matrices of the residual indicator scores of the endogenous
  constructs based on the target model in each repetition .r.

- `$Residuals_benchmark`:

  A list of matrices of the residual indicator scores of the endogenous
  constructs based on a model estimated by the procedure given to
  `.benchmark` for each repetition .r.

- `$Prediction_metrics`:

  A data frame containing the predictions metrics MAE, RMSE, Q2_predict,
  the misclassification error rate (MER), the MAPE, the MSE2, Theil's
  forecast accuracy (U1), Theil's forecast quality (U2), Bias proportion
  of MSE (UM), Regression proportion of MSE (UR), and disturbance
  proportion of MSE (UD) (Hora and Campos 2015; Watson and
  Teelucksingh 2002) .

- `$Information`:

  A list with elements `Target`, `Benchmark`,
  `Number_of_observations_training`, `Number_of_observations_test`,
  `Number_of_folds`, `Number_of_repetitions`, and
  `Handle_inadmissibles`.

## Details

The predict function implements the procedure introduced by Shmueli et
al. (2016) in the PLS context known as "PLSPredict" (Shmueli et al.
2019) including its variants PLScPredcit, OrdPLSpredict and
OrdPLScpredict. It is used to predict the indicator scores of endogenous
constructs and to evaluate the out-of-sample predictive power of a
model. For that purpose, the predict function uses k-fold
cross-validation to randomly split the data into training and test
datasets, and subsequently predicts the values of the test data based on
the model parameter estimates obtained from the training data. The
number of cross-validation folds is 10 by default but may be changed
using the `.cv_folds` argument. By default, the procedure is not
repeated (`.r = 1`). You may choose to repeat cross-validation by
setting a higher `.r` to be sure not to have a particular (unfortunate)
split. See Shmueli et al. (2019) for details. Typically `.r = 1` should
be sufficient though.

Alternatively, users may supply a test dataset as matrix or a data frame
of `.test_data` with the same column names as those in the data used to
obtain `.object` (the training data). In this case, arguments
`.cv_folds` and `.r` are ignored and predict uses the estimated
coefficients from `.object` to predict the values in the columns of
`.test_data`.

In Shmueli et al. (2016) PLS-based predictions for indicator `i` are
compared to the predictions based on a multiple regression of indicator
`i` on all available exogenous indicators (`.benchmark = "lm"`) and a
simple mean-based prediction summarized in the Q2_predict metric.
`predict()` is more general in that is allows users to compare the
predictions based on a so-called target model/specification to
predictions based on an alternative benchmark. Available benchmarks
include predictions based on a linear model, PLS-PM weights, unit
weights (i.e. sum scores), GSCA weights, PCA weights, and MAXVAR
weights.

Each estimation run is checked for admissibility using
[`verify()`](https://floschuberth.github.io/cSEM/reference/verify.md).
If the estimation yields inadmissible results, `predict()` stops with an
error (`"stop"`). Users may choose to `"ignore"` inadmissible results or
to simply set predictions to `NA` (`"set_NA"`) for the particular run
that failed.

## References

Hora J, Campos P (2015). “A review of performance criteria to validate
simulation models.” *Expert Systems*, **32**(5), 578–595.
[doi:10.1111/exsy.12111](https://doi.org/10.1111/exsy.12111) .  
  
Shmueli G, Ray S, Estrada JMV, Chatla SB (2016). “The Elephant in the
Room: Predictive Performance of PLS Models.” *Journal of Business
Research*, **69**(10), 4552–4564.
[doi:10.1016/j.jbusres.2016.03.049](https://doi.org/10.1016/j.jbusres.2016.03.049)
.  
  
Shmueli G, Sarstedt M, Hair JF, Cheah J, Ting H, Vaithilingam S, Ringle
CM (2019). “Predictive Model Assessment in PLS-SEM: Guidelines for Using
PLSpredict.” *European Journal of Marketing*, **53**(11), 2322–2347.
[doi:10.1108/ejm-02-2019-0189](https://doi.org/10.1108/ejm-02-2019-0189)
.  
  
Watson PK, Teelucksingh SS (2002). *A practical introduction to
econometric methods: Classical and modern*. University of West Indies
Press, Mona, Jamaica.

## See also

[csem](https://floschuberth.github.io/cSEM/reference/csem.md),
[cSEMResults](https://floschuberth.github.io/cSEM/reference/csem_results.md),
[`exportToExcel()`](https://floschuberth.github.io/cSEM/reference/exportToExcel.md)

## Examples

``` r
### Anime example taken from https://github.com/ISS-Analytics/pls-predict/

# Load data
data(Anime) # data is similar to the Anime.csv found on 
            # https://github.com/ISS-Analytics/pls-predict/ but with irrelevant
            # columns removed

# Split into training and data the same way as it is done on 
# https://github.com/ISS-Analytics/pls-predict/
set.seed(123)

index     <- sample.int(dim(Anime)[1], 83, replace = FALSE)
dat_train <- Anime[-index, ]
dat_test  <- Anime[index, ]

# Specify model
model <- "
# Structural model

ApproachAvoidance ~ PerceivedVisualComplexity + Arousal

# Measurement/composite model

ApproachAvoidance         =~ AA0 + AA1 + AA2 + AA3
PerceivedVisualComplexity <~ VX0 + VX1 + VX2 + VX3 + VX4
Arousal                   <~ Aro1 + Aro2 + Aro3 + Aro4
"

# Estimate (replicating the results of the `simplePLS()` function)
res <- csem(dat_train, 
            model, 
            .disattenuate = FALSE, # original PLS
            .iter_max = 300, 
            .tolerance = 1e-07, 
            .PLS_weight_scheme_inner = "factorial"
)

# Predict using a user-supplied training data set
pp <- predict(res, .test_data = dat_test)
#> Warning: The following warning occured in the `predict()` function:
#> Disattenuation is not applicable to benchmark `lm` and ignored.
pp
#> ________________________________________________________________________________
#> ----------------------------------- Overview -----------------------------------
#> 
#>  Number of obs. training            = 100
#>  Number of obs. test                = 83
#>  Number of cv folds                 = NA
#>  Number of repetitions              = 1
#>  Handle inadmissibles               = stop
#>  Estimator target                   = 'PLS-PM'
#>  Estimator benchmark                = 'lm'
#>  Disattenuation target              = 'FALSE'
#>  Disattenuation benchmark           = 'FALSE'
#>  Approach to predict                = 'earliest'
#> 
#> ------------------------------ Prediction metrics ------------------------------
#> 
#> 
#>   Name    MAE target  MAE benchmark  RMSE target RMSE benchmark   Q2_predict
#>   AA0         1.2125         1.1621       1.5575         1.5045       0.4625
#>   AA1         1.5319         1.5711       1.9005         1.9829       0.2794
#>   AA2         0.9891         0.9804       1.3993         1.4029       0.4396
#>   AA3         1.0564         1.0472       1.4434         1.4787       0.3656
#> ________________________________________________________________________________

### Compute prediction metrics  ------------------------------------------------
res2 <- csem(Anime, # whole data set
            model, 
            .disattenuate = FALSE, # original PLS
            .iter_max = 300, 
            .tolerance = 1e-07, 
            .PLS_weight_scheme_inner = "factorial"
)

# Predict using 10-fold cross-validation
if (FALSE) { # \dontrun{
pp2 <- predict(res, .benchmark = "lm")
pp2
## There is a plot method available
plot(pp2)} # }

### Example using OrdPLScPredict -----------------------------------------------
# Transform the numerical indicators into factors
if (FALSE) { # \dontrun{
data("BergamiBagozzi2000")
data_new <- data.frame(cei1    = as.ordered(BergamiBagozzi2000$cei1),
                       cei2    = as.ordered(BergamiBagozzi2000$cei2),
                       cei3    = as.ordered(BergamiBagozzi2000$cei3),
                       cei4    = as.ordered(BergamiBagozzi2000$cei4),
                       cei5    = as.ordered(BergamiBagozzi2000$cei5),
                       cei6    = as.ordered(BergamiBagozzi2000$cei6),
                       cei7    = as.ordered(BergamiBagozzi2000$cei7),
                       cei8    = as.ordered(BergamiBagozzi2000$cei8),
                       ma1     = as.ordered(BergamiBagozzi2000$ma1),
                       ma2     = as.ordered(BergamiBagozzi2000$ma2),
                       ma3     = as.ordered(BergamiBagozzi2000$ma3),
                       ma4     = as.ordered(BergamiBagozzi2000$ma4),
                       ma5     = as.ordered(BergamiBagozzi2000$ma5),
                       ma6     = as.ordered(BergamiBagozzi2000$ma6),
                       orgcmt1 = as.ordered(BergamiBagozzi2000$orgcmt1),
                       orgcmt2 = as.ordered(BergamiBagozzi2000$orgcmt2),
                       orgcmt3 = as.ordered(BergamiBagozzi2000$orgcmt3),
                       orgcmt5 = as.ordered(BergamiBagozzi2000$orgcmt5),
                       orgcmt6 = as.ordered(BergamiBagozzi2000$orgcmt6),
                       orgcmt7 = as.ordered(BergamiBagozzi2000$orgcmt7),
                       orgcmt8 = as.ordered(BergamiBagozzi2000$orgcmt8))

model <- "
# Measurement models
OrgPres =~ cei1 + cei2 + cei3 + cei4 + cei5 + cei6 + cei7 + cei8
OrgIden =~ ma1 + ma2 + ma3 + ma4 + ma5 + ma6
AffJoy  =~ orgcmt1 + orgcmt2 + orgcmt3 + orgcmt7
AffLove =~ orgcmt5 + orgcmt 6 + orgcmt8

# Structural model
OrgIden ~ OrgPres
AffLove ~ OrgIden
AffJoy  ~ OrgIden 
"
# Estimate using cSEM; note: the fact that indicators are factors triggers OrdPLSc
res <- csem(.model = model, .data = data_new[1:250,])
summarize(res)

# Predict using OrdPLSPredict
set.seed(123)
pred <- predict(
  .object = res, 
  .benchmark = "PLS-PM",
  .test_data = data_new[(251):305,],
   .treat_as_continuous = TRUE, .approach_score_target = "median"
  )

pred 
round(pred$Prediction_metrics[, -1], 4)} # }
```
