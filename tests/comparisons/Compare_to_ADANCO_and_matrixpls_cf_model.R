#===============================================================================
#
#   Topic: Compare results from cSEM (0.4.0.9000), ADANCO 2.2.1., 
#          and matrixpls (v. 1.0.13)
#
#   Type of model: common factor model
#
#   Date: 19.07.2021
#
#===============================================================================
rm(list = ls())

require(cSEM)
require(xlsx)
require(matrixpls)
data("satisfaction", package = "cSEM")

path_to_results <- "tests/comparisons/results_cf_model_adanco_2_2_1_centroid.xls"


### Import ADANCO results ------------------------------------------------------
## Import weights
weights_ADANCO <- read.xlsx(path_to_results, 
                            sheetName = "Weights", 
                            startRow = 2, endRow = 29, colIndex = c(2:8),
                            colClasses = c("character", rep("numeric", 6)))
weights_ADANCO[is.na(weights_ADANCO)] <- 0
rownames(weights_ADANCO) <- weights_ADANCO$Indicator
weights_ADANCO <- weights_ADANCO[, -1]

## Import loadings
loadings_ADANCO <- read.xlsx(path_to_results, 
                             sheetName = "Loadings", 
                             startRow = 2, endRow = 29, colIndex = c(2:8),
                             colClasses = c("character", rep("numeric", 6)))
loadings_ADANCO[is.na(loadings_ADANCO)] <- 0
rownames(loadings_ADANCO) <- loadings_ADANCO$Indicator
loadings_ADANCO <- loadings_ADANCO[, -1]

## Import path
path_ADANCO <- read.xlsx(path_to_results, 
                         sheetName = "Path Coefficients", 
                         startRow = 3, endRow = 8, colIndex = c(2:7),
                         colClasses = c("character", rep("numeric", 5)))
path_ADANCO[is.na(path_ADANCO)] <- 0
rownames(path_ADANCO) <- path_ADANCO$c..IMAG....EXPE....QUAL....VAL....SAT..
path_ADANCO <- as.matrix(path_ADANCO[, -1])

## Import standardized scores
scores_ADANCO <- read.xlsx(path_to_results, 
                           sheetName = "Standardized Construct Scores", 
                           startRow = 2, endRow = 252, colIndex = c(3:8),
                           colClasses = rep("numeric", 5))

## Import Goodness of model fit (saturated)
gof_s_ADANCO <- read.xlsx(path_to_results, 
                          sheetName = "Goodness of model fit (saturate", 
                          startRow = 2, endRow = 5, colIndex = c(2:5),
                          colClasses = c("character", rep("numeric", 3)))
gof_s_ADANCO1 <- gof_s_ADANCO$Value
names(gof_s_ADANCO1) <- gof_s_ADANCO$c..SRMR....dULS....dG..
gof_s_ADANCO <- gof_s_ADANCO1

## Import Goodness of model fit (estimated)
gof_e_ADANCO <- read.xlsx(path_to_results, 
                          sheetName = "Goodness of model fit (estimate", 
                          startRow = 2, endRow = 5, colIndex = c(2:5),
                          colClasses = c("character", rep("numeric", 3)))
gof_e_ADANCO1 <- gof_e_ADANCO$Value
names(gof_e_ADANCO1) <- gof_e_ADANCO$c..SRMR....dULS....dG..
gof_e_ADANCO <- gof_e_ADANCO1

## Reliability
reliability_ADANCO <- read.xlsx(path_to_results, 
                                sheetName = "Construct Reliability", 
                                startRow = 2, endRow = 8, colIndex = c(2:5))

## AVE
ave_ADANCO <- read.xlsx(path_to_results, 
                        sheetName = "Convergent Validity", 
                        startRow = 2, endRow = 8, colIndex = c(2:3))
ave_name <- ave_ADANCO$Construct
ave_ADANCO <- ave_ADANCO$Average.variance.extracted..AVE.
names(ave_ADANCO) <- ave_name

## HTMT
htmt_ADANCO <- read.xlsx(path_to_results, 
                              sheetName = "Discriminant Validity  Heterotr", 
                              startRow = 2, endRow = 8, colIndex = c(3:8))

htmt_ADANCO <- as.matrix(htmt_ADANCO)
htmt_ADANCO <- apply(htmt_ADANCO, 2, as.numeric)
rownames(htmt_ADANCO) <- colnames(htmt_ADANCO)
htmt_ADANCO[is.na(htmt_ADANCO)] <- 0

### Fornell-Larcker Criterion
fl_ADANCO <- read.xlsx(path_to_results, 
                         sheetName = "Fornell-Larcker Criterion", 
                         startRow = 2, endRow = 8, colIndex = c(3:8))

fl_ADANCO <- as.matrix(fl_ADANCO)
fl_ADANCO <- apply(fl_ADANCO, 2, as.numeric)
fl_ADANCO[upper.tri(fl_ADANCO)] <- t(fl_ADANCO)[upper.tri(fl_ADANCO)]
rownames(fl_ADANCO) <- colnames(fl_ADANCO)

## R2
r2_ADANCO <- read.xlsx(path_to_results, 
                         sheetName = "R-Squared", 
                         startRow = 2, endRow = 7, colIndex = c(2:4))

## Cohens f2 (effect size)
f2_ADANCO <- read.xlsx(path_to_results, 
                       sheetName = "Effect Overview", 
                       startRow = 3, endRow = 17, colIndex = 6, header = FALSE,
                       colClasses = "numeric")
f2_ADANCO <- as.numeric(f2_ADANCO$X6)[!is.na(as.numeric(f2_ADANCO$X6))]

## Proxy cor
construct_cor_ADANCO <- read.xlsx(path_to_results, 
                              sheetName = "Inter-Construct Correlations", 
                              startRow = 2, endRow = 8, colIndex = c(3:8))
construct_cor_ADANCO <- as.matrix(construct_cor_ADANCO)
construct_cor_ADANCO[upper.tri(construct_cor_ADANCO)] <- t(construct_cor_ADANCO)[upper.tri(construct_cor_ADANCO)]
construct_cor_ADANCO <- apply(construct_cor_ADANCO, 2, as.numeric)
rownames(construct_cor_ADANCO) <- colnames(construct_cor_ADANCO)

## Model-implied indicator cor (saturated)
sigma_hat_saturated <- read.xlsx(path_to_results, 
                                 sheetName = "Impl_Cor Saturated Model", 
                                 startRow = 2, endRow = 29, colIndex = c(3:29))
sigma_hat_saturated <- as.matrix(sigma_hat_saturated)
rownames(sigma_hat_saturated) <- colnames(sigma_hat_saturated)

## Model-implied indicator cor (estimated)
sigma_hat_estimated <- read.xlsx(path_to_results, 
                                 sheetName = "Impl_Cor Estimated Model", 
                                 startRow = 2, endRow = 29, colIndex = c(3:29))
sigma_hat_estimated <- as.matrix(sigma_hat_estimated)
rownames(sigma_hat_estimated) <- colnames(sigma_hat_estimated)

### Load data and define model =================================================
# Model cSEM
model <- "
# Structural model
EXPE ~ IMAG
QUAL ~ EXPE
VAL  ~ EXPE + QUAL
SAT  ~ IMAG + EXPE + QUAL + VAL
LOY  ~ IMAG + SAT

# Measurement model

IMAG =~ imag1 + imag2 + imag3 + imag4 + imag5
EXPE =~ expe1 + expe2 + expe3 + expe4 + expe5
QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
VAL  =~ val1  + val2  + val3  + val4
SAT  =~ sat1  + sat2  + sat3  + sat4
LOY  =~ loy1  + loy2  + loy3  + loy4
"

### cSEM -----------------------------------------------------------------------

a1 <- csem(
  .data                        = satisfaction,
  .model                       = model,
  .approach_weights            = "PLS-PM",
  .tolerance                   = 1e-06,
  .PLS_modes                   = NULL,
  .PLS_ignore_structural_model = FALSE,
  .PLS_weight_scheme_inner     = "centroid" 
)
suma1 <- cSEM:::summarize(a1)
assa1 <- assess(a1)

### Matrixpls ------------------------------------------------------------------
a2 <- matrixpls(S           = cor(satisfaction),
                model       = model,
                standardize = TRUE,
                parametersReflective = estimator.plscLoadings,
                disattenuate = TRUE,
                innerEstim  = innerEstim.centroid,
                ignoreInnerModel = FALSE,
                outerEstim  = outerEstim.modeA,
                tol = 1e-06
)
suma2 <- summary(a2)
resid_matpls <- matrixpls:::residuals.matrixpls(a2, observed = FALSE)

W_matpls <- attr(a2, "W")
Lambda_matpls <- attr(a2, "reflective")
path_matpls <- attr(a2, "inner")
S_matpls <- attr(a2, "S")
vcv_proxy_matpls <- W_matpls %*% S_matpls %*% t(W_matpls)
vcv_construct_matpls <- attr(a2, "C")

Theta_matpls <- diag(1 - diag(Lambda_matpls %*% t(Lambda_matpls)))
sigma_hat_saturated_matpls <- Lambda_matpls %*% vcv_construct_matpls %*% t(Lambda_matpls) + Theta_matpls

sigma_hat_estimated_matpls <- matrixpls:::fitted.matrixpls(a2)
### Compare ====================================================================
# Note (19.12.2019): Apparently ADANCO has a bug. If you change the weighting scheme
# from Factor to Centroid nothing changes, i.e., the ADANCO always uses the
# factor scheme.
# Note (19.07.2021): This issue has been fixed in version 2.2.1. Now centroid and
# factorial scheme yield different results

## Compare weights
a1$Estimates$Weight_estimates - t(weights_ADANCO) # identical at the 8th sig digit
a1$Estimates$Weight_estimates - W_matpls # identical

## Compare loadings
a1$Estimates$Loading_estimates - t(loadings_ADANCO) # identical at the 8th sig digit
a1$Estimates$Loading_estimates - t(Lambda_matpls) # identical

## Compare path
a1$Estimates$Path_estimates[-1, -6] - t(path_ADANCO) # identical at the 8th sig digit
a1$Estimates$Path_estimates - path_matpls # identical

## Compare proxy/composite correlation matrix (C in cSEM)
a1$Estimates$Proxy_VCV - cor(scores_ADANCO) # identical at the 8th sig digit
a1$Estimates$Proxy_VCV - vcv_proxy_matpls # identical

## Compare construct correlation matrix (P in cSEM)
a1$Estimates$Construct_VCV - construct_cor_ADANCO # identical at the 8th sig digit
a1$Estimates$Construct_VCV - vcv_construct_matpls # identical

## Compare model_implied indicator cor (saturated)
round(fit(a1, .saturated = TRUE) - sigma_hat_saturated, 8) # identical at the 8th sig digit
fit(a1, .saturated = TRUE) - sigma_hat_saturated_matpls # identical

## Compare model_implied indicator cor (estimated)
round(fit(a1) - sigma_hat_estimated, 6) # Identical
round(fit(a1) - sigma_hat_estimated_matpls, 6) # Differences for blocks of SAT and LOY
                                               # See: fit function for explanation why.

## Compare model_implied indicator cor of constructs
# We have:
#                             Sigma = Lambda V_eta Lambda' + Theta
#   Lambda' (Sigma - Theta ) Lambda =  (Lambda'Lambda) V_eta (Lambda' Lambda)
# V_eta = (Lambda' Lambda)^-1 Lambda' (Sigma - Theta ) Lambda (Lambda' Lambda)^-1
#
# V_eta matrixpls
V_eta_matpls <- solve(t(Lambda_matpls) %*% Lambda_matpls) %*% t(Lambda_matpls) %*% 
  (sigma_hat_estimated_matpls - Theta_matpls) %*% Lambda_matpls %*% 
  solve(t(Lambda_matpls) %*% Lambda_matpls)

round(fit(a1, .type_vcv = "construct")  - V_eta_matpls, 8) 
# main diagonal elements differ, i.e., model-implied construct variances are not
# 1 (which they should be). This causes the blocks of those constructs that
# dont have a variance of 1 to be wrong in "sigma_hat_estimated_matpls".

## V_eta for ADANCO
V_eta_ADANCO <- solve(t(Lambda_matpls) %*% Lambda_matpls) %*% t(Lambda_matpls) %*% 
  (sigma_hat_estimated - Theta_matpls) %*% Lambda_matpls %*% 
  solve(t(Lambda_matpls) %*% Lambda_matpls)

round(fit(a1, .type_vcv = "construct") - V_eta_ADANCO, 8) # identical to the 8th digit

### Quality criteria -------------------------------------------------------------
## Overall model fit
# SRMR, dG dL
# Note : the fitted function in ADANCO is not correct. Reason: the model-implied
# construct correlation matrix does not always have ones on its main diagonal,
# i.e. the variance of some constructs (those deeper down the nomologial net).
# Note (19.07.2021): fixed in ADANCO version 2.2.1

## Distance measures for saturated model
cSEM:::calculateSRMR(a1, .saturated = TRUE) - gof_s_ADANCO["SRMR"] # identical at the 8th sig digit
cSEM:::calculateDL(a1, .saturated = TRUE) - gof_s_ADANCO["dULS"] # identical at the 8th sig digit
cSEM:::calculateDG(a1, .saturated = TRUE) - gof_s_ADANCO["dG"] # identical at the 8th sig digit
# no measures for matrixpls

## Distance measures for estimated model
cSEM:::calculateSRMR(a1) - gof_e_ADANCO["SRMR"] # identical at the 8th sig digit
cSEM:::calculateDL(a1) - gof_e_ADANCO["dULS"] # identical at the 8th sig digit
cSEM:::calculateDG(a1) - gof_e_ADANCO["dG"] # identical at the 8th sig digit

cSEM:::calculateSRMR(a1) - resid_matpls$indices[1] # different, since model-implied vcv is different
# for dg and DL no measure for matrixpls

## R2
assa1$R2 - r2_ADANCO$Coefficient.of.determination..R2. # identical at the 8th sig digit
assa1$R2 - suma2$r2[-1] # identical to matrixpls

## R2 adjusted
assa1$R2_adj - r2_ADANCO$Adjusted.R2 # identical at the 8th sig digit
# no measures for matrixpls

## AVE
assa1$AVE - ave_ADANCO # identical at the 8th sig digit
assa1$AVE - suma2$ave$ave # identical to matrixpls

## HTMT
assa1$HTMT - htmt_ADANCO # identical to ADANCO
assa1$HTMT - suma2$htmt # identical to matrixpls

## Fornell-Larcker
assa1$`Fornell-Larcker` - fl_ADANCO # identical to ADANCO 
# no measure for matrixpls

## Cohens f^2
f2_ADANCO
assa1$F2 # identical to ADANCO
# no measure for matrixpls

## Reliability
assa1$Reliability$`Dijkstra-Henselers_rho_A` - reliability_ADANCO$Dijkstra.Henseler.s.rho...U.03C1.A. # identical to ADNACO
assa1$Reliability$Joereskogs_rho - reliability_ADANCO$Jöreskog.s.rho...U.03C1.c. # identical to ADANCO
assa1$Reliability$Cronbachs_alpha - reliability_ADANCO$Cronbach.s.alpha.a. # identical ADANCO

assa1$Reliability$Joereskogs_rho - suma2$cr # identical to matrixpls

## GoF
assa1$GoF -  suma2$gof # identical to matrixpls
# no measure for ADANCO