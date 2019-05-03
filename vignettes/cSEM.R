## ---- eval = FALSE-------------------------------------------------------
#  model <- "
#  ## Structural model
#  eta2 ~ eta1
#  
#  ## Measurement model
#  eta1 <~ item1 + item2 + item3 # eta1 is modeled as a composite
#  eta2 =~ item4 + item5 + item6 # eta2 is modeled as a common factor
#  "

## ----eval=FALSE----------------------------------------------------------
#  model <- "
#  ## Structural model
#  eta2 ~ eta1
#  
#  ## Measurement model
#  eta1 <~ item1 + item2 + item3
#  eta2 =~ item4 + item5 + item6
#  "
#  
#  # Identical
#  csem(threecommonfactors, model)
#  csem(threecommonfactors, model, .disattenuate = TRUE)
#  
#  # To supress automatic disattenuation
#  csem(threecommonfactors, model, .disattenuate = FALSE)

## ----eval=FALSE----------------------------------------------------------
#  model <- "
#  ## Structural model
#  eta2 ~ eta1
#  
#  ## Measurement model (pure composite model)
#  eta1 <~ item1 + item2 + item3
#  eta2 <~ item4 + item5 + item6
#  "
#  
#  ### Currently the following weight approaches are implemented
#  # Partial least squares path modeling (PLS)
#  csem(threecommonfactors, model, .approach_weights = "PLS-PM") # default
#  
#  # Generalized canonical correlation analysis (Kettenring approaches)
#  csem(threecommonfactors, model, .approach_weights = "SUMCORR")
#  csem(threecommonfactors, model, .approach_weights = "MAXVAR")
#  csem(threecommonfactors, model, .approach_weights = "SSQCORR")
#  csem(threecommonfactors, model, .approach_weights = "MINVAR")
#  csem(threecommonfactors, model, .approach_weights = "GENVAR")
#  
#  # Generalized structured component analysis (GSCA)
#  csem(threecommonfactors, model, .approach_weights = "GSCA")
#  
#  # Principal component analysis (PCA)
#  csem(threecommonfactors, model, .approach_weights = "PCA")
#  
#  # Factor score regression (FSR) using "unit", "bartlett" or "regression" weights
#  csem(threecommonfactors, model, .approach_weights = "unit")
#  csem(threecommonfactors, model, .approach_weights = "bartlett")
#  csem(threecommonfactors, model, .approach_weights = "regression")

## ----eval=FALSE----------------------------------------------------------
#  model <- "
#  # Structural model
#  EXPE ~ IMAG
#  QUAL ~ EXPE
#  VAL  ~ EXPE + QUAL
#  SAT  ~ IMAG + EXPE + QUAL + VAL
#  LOY  ~ IMAG + SAT
#  
#  # Measurement model
#  
#  IMAG <~ imag1 + imag2 + imag3                  # composite
#  EXPE <~ expe1 + expe2 + expe3                  # composite
#  QUAL <~ qual1 + qual2 + qual3 + qual4 + qual5  # composite
#  VAL  <~ val1  + val2  + val3                   # composite
#  SAT  =~ sat1  + sat2  + sat3  + sat4           # common factor
#  LOY  =~ loy1  + loy2  + loy3  + loy4           # common factor
#  "

## ------------------------------------------------------------------------
model <- "
# Structural model
EXPE ~ IMAG + IMAG.IMAG

# Measurement model
EXPE <~ expe1 + expe2
IMAG <~ imag1 + imag2
"


## ---- eval=FALSE---------------------------------------------------------
#  model <- "
#  # Structural model
#  SAT ~ QUAL
#  VAL ~ SAT + QUAL
#  
#  # Measurement model
#  EXPE <~ expe1 + expe2
#  
#  SAT =~ sat1 + sat2
#  VAL =~ val1 + val2
#  QUAL =~ IMAG + EXPE
#  IMAG <~ imag1 + imag2
#  "

## ----warning=FALSE, message=FALSE----------------------------------------
require(cSEM)

model <- "
# Path model / Regressions 
eta2 ~ eta1
eta3 ~ eta1 + eta2

# Measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"
  
a <- csem(.data = threecommonfactors, .model = model)
a

## ---- eval=FALSE---------------------------------------------------------
#  csem(
#     .data                        = threecommonfactors,
#     .model                       = model,
#     .approach_cor_robust         = "none",
#     .approach_nl                 = "sequential",
#     .approach_paths              = "OLS",
#     .approach_weights            = "PLS-PM",
#     .conv_criterion              = "diff_absolute",
#     .disattenuate                = TRUE,
#     .dominant_indicators         = NULL,
#     .estimate_structural         = TRUE,
#     .id                          = NULL,
#     .iter_max                    = 100,
#     .normality                   = TRUE,
#     .PLS_approach_cf             = "dist_squared_euclid",
#     .PLS_ignore_structural_model = FALSE,
#     .PLS_modes                   = NULL,
#     .PLS_weight_scheme_inner     = "path",
#     .reliabilities               = NULL,
#     .starting_values             = NULL,
#     .tolerance                   = 1e-05,
#     .resample_method             = "none",
#     .resample_method2            = "none",
#     .R                           = 499,
#     .R2                          = 199,
#     .handle_inadmissibles        = "drop",
#     .user_funs                   = NULL,
#     .eval_plan                   = "sequential",
#     .seed                        = NULL,
#     .sign_change_option          = "no"
#      )

## ----eval=FALSE----------------------------------------------------------
#  # Setting `.resample_method`
#  b1 <- csem(.data = satisfaction, .model = model, .resample_method = "bootstrap")
#  b2 <- resamplecSEMResults(a)

## ----eval=FALSE----------------------------------------------------------
#  summarize(b1)
#  infer(b1, .quantity = c("CI_standard_z", "CI_percentile")) # no print method yet

## ----eval=FALSE----------------------------------------------------------
#  b <- csem(
#    .data            = satisfaction,
#    .model           = model,
#    .resample_method = "bootstrap",
#    .R               = 999,
#    .seed            = 98234,
#    .eval_plan       = "multiprocess")

