### Example from Dijkstra & Hensler (2015)
## Prepartion (values are from p. 15-16 of the paper)
Lambda <- t(kronecker(diag(6), c(0.7, 0.7, 0.7)))
Phi <- matrix(c(1.0000, 0.5000, 0.5000, 0.5000, 0.0500, 0.4000, 
                0.5000, 1.0000, 0.5000, 0.5000, 0.5071, 0.6286,
                0.5000, 0.5000, 1.0000, 0.5000, 0.2929, 0.7714,
                0.5000, 0.5000, 0.5000, 1.0000, 0.2571, 0.6286,
                0.0500, 0.5071, 0.2929, 0.2571, 1.0000, sqrt(0.5),
                0.4000, 0.6286, 0.7714, 0.6286, sqrt(0.5), 1.0000), 
              ncol = 6)

## Create population indicator covariance matrix
Sigma <- t(Lambda) %*% Phi %*% Lambda
diag(Sigma) <- 1
dimnames(Sigma) <- list(paste0("x", rep(1:6, each = 3), 1:3),
                        paste0("x", rep(1:6, each = 3), 1:3))

## Generate data
dat <- MASS::mvrnorm(n = 500, mu = rep(0, 18), Sigma = Sigma, empirical = TRUE)
# empirical = TRUE to show that 2SLS is in fact able to recover the true population
# parameters.

## Model to estimate
model <- "
## Structural model (nonrecurisve)
eta5 ~ eta6 + eta1 + eta2
eta6 ~ eta5 + eta3 + eta4

## Measurement model
eta1 =~ x11 + x12 + x13
eta2 =~ x21 + x22 + x23
eta3 =~ x31 + x32 + x33
eta4 =~ x41 + x42 + x43

eta5 =~ x51 + x52 + x53
eta6 =~ x61 + x62 + x63
"

library(cSEM)

## Estimate
res_ols <- csem(dat, .model = model, .approach_paths = "OLS")
sum_res_ols <- summarize(res_ols) 

# Note: For the example the model-implied indicator correlation is irrelevant
#       the warnings can be ignored.

res_2sls <- csem(dat, .model = model, .approach_paths = "2SLS",
                 .instruments = list("eta5" = c('eta1','eta2','eta3','eta4'), 
                                     "eta6" = c('eta1','eta2','eta3','eta4')))
sum_res_2sls <- summarize(res_2sls)
# Note that exogenous constructs are supplied as instruments for themselves!

## Test for endogeneity
test_ha <- testHausman(res_2sls, .R = 200)
test_ha
