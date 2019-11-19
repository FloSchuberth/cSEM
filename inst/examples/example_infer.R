\donttest{
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
  
## Estimate the model with bootstrap resampling 
a <- csem(satisfaction, model, .resample_method = "bootstrap", .R = 50)

## Compute inferential quantities
inf <- infer(a)

inf$Path_estimates$CI_basic
inf$Indirect_effect$sd

### To compute the bias-corrected and accelerated and/or the studentized t-inverval
### confidence interval:
inf <- infer(a, .quantity = c("all", "CI_bca")) # requires jackknife estimates

## For the studentied t-interval confidence interval, a double bootstrap is required:
## Estimate the model with double bootstrap resampling:
# Notes:
#   1. The .resample_method2 arguments triggers a bootstrap of each bootstrap sample
#   2. The double bootstrap is is very time consuming, consider setting `.eval_plan = "multiprocess`. 
#      To speed things up for the example, .R and .R2 are reduced. Results are
#      probably quite
a1 <- csem(satisfaction, model, .resample_method = "bootstrap", .R = 40,
          .resample_method2 = "bootstrap", .R2 = 20, .handle_inadmissibles = "replace") 
infer(a1, .quantity = "CI_t_interval") 
}