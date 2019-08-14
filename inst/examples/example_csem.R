# ===========================================================================
# Basic usage
# ===========================================================================
### Linear model ------------------------------------------------------------
# Most basic usage requires a dataset and a model. We use the 
#  `threecommonfactors` dataset. 

## Take a look at the dataset
#?threecommonfactors

## Specify the (correct) model
model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

## Estimate
res <- csem(threecommonfactors, model)

## Postestimation
verify(res)
summarize(res)  
assess(res)

# Notes: 
#   1. By default no inferential quantities (e.g. Std. errors, p-values, or
#      confidence intervals) are calculated. Use resampling to obtain
#      inferential quantities. See "Resampling" in the "Extended usage"
#      section below.
#   2. `summarize()` prints the full output by default. For a more condensed
#       output use:
print(summarize(res), .full_output = FALSE)

## Dealing with endogeneity -------------------------------------------------

### Nonlinear model ---------------------------------------------------------

### Second order model ------------------------------------------------------
## Take a look at the dataset
#?dgp_2ndorder_cf_of_c

model <- "
# Path model / Regressions 
c4   ~ eta1
eta2 ~ eta1 + c4

# Reflective measurement model
c1   <~ y11 + y12 
c2   <~ y21 + y22 + y23 + y24
c3   <~ y31 + y32 + y33 + y34 + y35 + y36 + y37 + y38
eta1 =~ y41 + y42 + y43
eta2 =~ y51 + y52 + y53

# Composite model
c4   =~ c1 + c2 + c3
"

res_2stage <- csem(dgp_2ndorder_cf_of_c, model, .approach_2ndorder = "2stage")
res_mixed  <- csem(dgp_2ndorder_cf_of_c, model, .approach_2ndorder = "mixed")

# The standard repeated indicators approach is done by 1.) respecifying the model
# and 2.) adding the repeated indicators to the data set

# 1.) Respecify the model
model_RI <- "
# Path model / Regressions 
c4   ~ eta1
eta2 ~ eta1 + c4
c4   ~ c1 + c2 + c3

# Reflective measurement model
c1   <~ y11 + y12 
c2   <~ y21 + y22 + y23 + y24
c3   <~ y31 + y32 + y33 + y34 + y35 + y36 + y37 + y38
eta1 =~ y41 + y42 + y43
eta2 =~ y51 + y52 + y53

# c4 is a common factor measured by composites
c4 =~ y11_temp + y12_temp + y21_temp + y22_temp + y23_temp + y24_temp +
      y31_temp + y32_temp + y33_temp + y34_temp + y35_temp + y36_temp + 
      y37_temp + y38_temp
"

# 2.) Update data set
data_RI <- dgp_2ndorder_cf_of_c
coln <- c(colnames(data_RI), paste0(colnames(data_RI), "_temp"))
data_RI <- data_RI[, c(1:ncol(data_RI), 1:ncol(data_RI))]
colnames(data_RI) <- coln

# Estimate
res_RI <- csem(data_RI, model_RI)
summarize(res_RI)

### Multigroup analysis -----------------------------------------------------

# ===========================================================================
# Extended usage
# ===========================================================================
# `csem()` provides defaults for all arguments except `.data` and `.model`.
#   Below some common options/tasks that users are likely to be interested in.
#   We use the threecommonfactors data set again:

model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# (Reflective) measurement model
eta1 =~ y11 + y12 + y13
eta2 =~ y21 + y22 + y23
eta3 =~ y31 + y32 + y33
"

### PLS vs PLSc and disattenuation
# In the model all concepts are modeled as common factors. If 
#   .approach_weights = "PLS-PM", csem() uses PLSc to disattenuate composite-indicator 
#   and composite-composite correlations.
res_plsc <- csem(threecommonfactors, model, .approach_weights = "PLS-PM")
res$Information$Model$construct_type # all common factors

# To obtain "original" (inconsistent) PLS estimates use `.disattenuate = FALSE`
res_pls <- csem(threecommonfactors, model, 
                .approach_weights = "PLS-PM",
                .disattenuate = FALSE
                )

s_plsc <- summarize(res_plsc)
s_pls  <- summarize(res_pls)

# Compare
data.frame(
  "Path"      = s_plsc$Estimates$Path_estimates$Name,
  "Pop_value" = c(0.6, 0.4, 0.35), # see ?threecommonfactors
  "PLSc"      = s_plsc$Estimates$Path_estimates$Estimate,
  "PLS"       = s_pls$Estimates$Path_estimates$Estimate
  )

### Resampling --------------------------------------------------------------
## Basic resampling
res_boot <- csem(threecommonfactors, model, .resample_method = "bootstrap", .R = 40)
res_jack <- csem(threecommonfactors, model, .resample_method = "jackknife")

# See ?resamplecSEMResults for more examples

### Choosing a different weightning scheme --------------------------------------------

res_gscam  <- csem(threecommonfactors, model, .approach_weights = "GSCA")
res_gsca   <- csem(threecommonfactors, model, 
                  .approach_weights = "GSCA",
                  .disattenuate = FALSE
                  )

s_gscam <- summarize(res_gscam)
s_gsca  <- summarize(res_gsca)

# Compare
data.frame(
  "Path"      = s_gscam$Estimates$Path_estimates$Name,
  "Pop_value" = c(0.6, 0.4, 0.35), # see ?threecommonfactors
  "GSCAm"      = s_gscam$Estimates$Path_estimates$Estimate,
  "GSCA"       = s_gsca$Estimates$Path_estimates$Estimate
)
### Fine-tuning a weighting scheme ------------------------------------------
## Setting starting values

sv <- list("eta1" = c("y12" = 10, "y13" = 4, "y11" = 1))
res <- csem(threecommonfactors, model, .starting_values = sv)

## Choosing a different inner weighting scheme 
#?args_csem_dotdotdot

res <- csem(threecommonfactors, model, .PLS_weight_scheme_inner = "factorial",
            .PLS_ignore_structural_model = TRUE)


## Choosing different modes for PLS
# By default, concepts modeled as common factors uses PLS Mode A weights.
modes <- list("eta1" = "unit", "eta2" = "modeB", "eta3" = "unit")
res   <- csem(threecommonfactors, model, .PLS_modes = modes)
summarize(res)

