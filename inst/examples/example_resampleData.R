# ===========================================================================
# Using the raw data 
# ===========================================================================
### Bootstrap (default) -----------------------------------------------------

res_boot1 <- resampleData(.data = satisfaction)
str(res_boot1, max.level = 3, list.len = 3)

## To replicate a bootstrap draw use .seed:
res_boot1a <- resampleData(.data = satisfaction, .seed = 2364)
res_boot1b <- resampleData(.data = satisfaction, .seed = 2364)
                           
identical(res_boot1, res_boot1a) # TRUE

### Jackknife ---------------------------------------------------------------

res_jack <- resampleData(.data = satisfaction, .resample_method = "jackknife")
str(res_jack, max.level = 3, list.len = 3)

### Cross-validation --------------------------------------------------------
## Create dataset for illustration:
dat <- data.frame(
  "x1" = rnorm(100),
  "x2" = rnorm(100),
  "group" = sample(c("male", "female"), size = 100, replace = TRUE),
  stringsAsFactors = FALSE)

## 10-fold cross-validation (repeated 100 times)
cv_10a <- resampleData(.data = dat, .resample_method = "cross-validation", 
                      .R = 100)
str(cv_10a, max.level = 3, list.len = 3)

# Cross-validation can be done by group if a group identifyer is provided:
cv_10 <- resampleData(.data = dat, .resample_method = "cross-validation", 
                      .id = "group", .R = 100)

## Leave-one-out-cross-validation (repeated 50 times)
cv_loocv  <- resampleData(.data = dat[, -3], 
                          .resample_method = "cross-validation", 
                          .cv_folds = nrow(dat),
                          .R = 50)
str(cv_loocv, max.level = 2, list.len = 3)

### Permuation ---------------------------------------------------------------

res_perm <- resampleData(.data = dat, .resample_method = "permutation",
                         .id = "group")
str(res_perm, max.level = 2, list.len = 3)

# Forgetting to set .id causes an error
\dontrun{
res_perm <- resampleData(.data = dat, .resample_method = "permutation")
}

# ===========================================================================
# Using a cSEMResults object
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
a <- csem(satisfaction, model)

# Create bootstrap and jackknife samples
res_boot <- resampleData(a, .resample_method = "bootstrap", .R = 499)
res_jack <- resampleData(a, .resample_method = "jackknife")

# Since `satisfaction` is the dataset used the following approaches yield
# identical results.
res_boot_data   <- resampleData(.data = satisfaction, .seed = 2364)
res_boot_object <- resampleData(a, .seed = 2364)

identical(res_boot_data, res_boot_object) # TRUE