## ----eval = FALSE-------------------------------------------------------------
#  csem(.data = my_data, .model = my_model)

## -----------------------------------------------------------------------------
model <- "
# Structural model
EXPE ~ IMAG

# Reflective measurement model
EXPE =~ expe1 + expe2
IMAG =~ imag1 + imag2
"

## ----eval=FALSE---------------------------------------------------------------
#  model <- "
#  # Structural model
#  EXPE ~ IMAG
#  QUAL ~ EXPE
#  VAL  ~ EXPE + QUAL
#  SAT  ~ IMAG + EXPE + QUAL + VAL
#  LOY  ~ IMAG + SAT
#  
#  # Composite model
#  IMAG <~ imag1 + imag2 + imag3                  # composite
#  EXPE <~ expe1 + expe2 + expe3                  # composite
#  QUAL <~ qual1 + qual2 + qual3 + qual4 + qual5  # composite
#  VAL  <~ val1  + val2  + val3                   # composite
#  
#  # Reflective measurement model
#  SAT  =~ sat1  + sat2  + sat3  + sat4           # common factor
#  LOY  =~ loy1  + loy2  + loy3  + loy4           # common factor
#  
#  # Measurement error correlation
#  sat1 ~~ sat2
#  "

## -----------------------------------------------------------------------------
model <- "
# Structural model
EXPE ~ IMAG + IMAG.IMAG

# Composite model
EXPE <~ expe1 + expe2
IMAG <~ imag1 + imag2
"

## ----eval=FALSE---------------------------------------------------------------
#  model <- "
#  # Structural model
#  SAT ~ QUAL
#  VAL ~ SAT + QUAL
#  
#  # Reflective measurement model
#  SAT  =~ sat1 + sat2
#  VAL  =~ val1 + val2
#  
#  # Composite model
#  IMAG <~ imag1 + imag2
#  EXPE <~ expe1 + expe2
#  
#  # Second-order term
#  QUAL =~ IMAG + EXPE
#  "

## ----warning=FALSE, message=FALSE---------------------------------------------
require(cSEM)

model <- "
# Path model / Regressions
eta2 ~ eta1
eta3 ~ eta1 + eta2

# Reflective measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

a <- csem(.data = threecommonfactors, .model = model)
a

## ----eval=FALSE---------------------------------------------------------------
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
#     .normality                   = FALSE,
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

## ----echo=FALSE, include=FALSE------------------------------------------------
x <- runif(1) # to intialize .Random.seed

## -----------------------------------------------------------------------------
b1 <- csem(.data = threecommonfactors, .model = model, .resample_method = "bootstrap")
b2 <- resamplecSEMResults(a)

## -----------------------------------------------------------------------------
summarize(b1)

## -----------------------------------------------------------------------------
ii <- infer(b1, .quantity = c("CI_standard_z", "CI_percentile"), .alpha = c(0.01, 0.05))
ii$Path_estimates

## ----eval=FALSE---------------------------------------------------------------
#  b <- csem(
#    .data            = satisfaction,
#    .model           = model,
#    .resample_method = "bootstrap",
#    .R               = 999,
#    .seed            = 98234,
#    .eval_plan       = "multiprocess")
#  
#  # Output omitted

## ----eval=FALSE---------------------------------------------------------------
#  model <- "
#  ## Structural model
#  eta2 ~ eta1
#  
#  ## Measurement model
#  eta1 <~ y11 + y12 + y13
#  eta2 =~ y21 + y22 + y23
#  "
#  
#  # Identical
#  csem(threecommonfactors, model)
#  csem(threecommonfactors, model, .disattenuate = TRUE)
#  
#  # To supress automatic disattenuation
#  csem(threecommonfactors, model, .disattenuate = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  model <- "
#  ## Structural model
#  eta2 ~ eta1
#  
#  ## Composite model
#  eta1 <~ y11 + y12 + y13
#  eta2 <~ y21 + y22 + y23
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

