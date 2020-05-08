\dontrun{
# Note: example not run as resampling is time consuming
# ===========================================================================
# Basic usage
# ===========================================================================
model <- "
# Structural model
QUAL ~ EXPE
EXPE ~ IMAG
SAT  ~ IMAG + EXPE + QUAL + VAL
LOY  ~ IMAG + SAT
VAL  ~ EXPE + QUAL

# Measurement model
EXPE =~ expe1 + expe2 + expe3 + expe4 + expe5
IMAG =~ imag1 + imag2 + imag3 + imag4 + imag5
LOY  =~ loy1  + loy2  + loy3  + loy4
QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
SAT  =~ sat1  + sat2  + sat3  + sat4
VAL  =~ val1  + val2  + val3  + val4
"

## Estimate the model without resampling 
a <- csem(satisfaction, model)

## Bootstrap and jackknife estimation
boot <- resamplecSEMResults(a)
jack <- resamplecSEMResults(a, .resample_method = "jackknife") 

## Alternatively use .resample_method in csem()
boot_csem <- csem(satisfaction, model, .resample_method = "bootstrap")
jack_csem <- csem(satisfaction, model, .resample_method = "jackknife")

# ===========================================================================
# Extended usage
# ===========================================================================
### Double resampling  ------------------------------------------------------
# The confidence intervals (e.g. the bias-corrected and accelearated CI) 
# require double resampling. Use .resample_method2 for this.

boot1 <- resamplecSEMResults(
  .object = a, 
  .resample_method  = "bootstrap", 
  .R                = 50,
  .resample_method2 = "bootstrap", 
  .R2               = 20,
  .seed             = 1303
  )

## Again, this is identical to using csem 
boot1_csem <- csem(
  .data             = satisfaction, 
  .model            = model, 
  .resample_method  = "bootstrap",
  .R                = 50,
  .resample_method2 = "bootstrap",
  .R2               = 20,
  .seed             = 1303
  )

identical(boot1, boot1_csem) # only true if .seed was set

### Inference ---------------------------------------------------------------
# To get inferencial quanitites such as the estimated standard error or
# the percentile confidence intervall for each resampled quantity use 
# postestimation function infer()

inference <- infer(boot1)
inference$Path_estimates$sd
inference$Path_estimates$CI_percentile

# As usual summarize() can be called directly
summarize(boot1)

# In the example above .R x .R2 = 50 x 20 = 1000. Multiprocessing will be
# faster on most systems here and is therefore recommended. Note that multiprocessing
# does not affect the random number generation

boot2 <- resamplecSEMResults(
  .object           = a, 
  .resample_method  = "bootstrap", 
  .R                = 50,
  .resample_method2 = "bootstrap", 
  .R2               = 20,
  .eval_plan        = "multiprocess", 
  .seed             = 1303
  )

identical(boot1, boot2)}