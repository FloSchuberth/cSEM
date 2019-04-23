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

