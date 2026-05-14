library(SimDesign)
library(rlang)
rlang::global_handle()
library(lobstr)
cli::pretty_print_code()
library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(2)
library(tibble)
library(future)
library(future.mirai)
plan(future.mirai::mirai_multisession)
library(future.apply)
library(boot)

# PoC --------------------------------------------------------------------

sigma <- 1
x1 <- rnorm(1000)
x2 <- rnorm(1000)
allDat <- tibble::tibble(x1, x2)
trueDat <-subset(allDat, select = c(x1, x2))
beta <- .3
y <- rnorm(1000, mean = rowSums(beta*trueDat), sd = sigma)

allDat$y <- y

lm(y ~ ., data = allDat) |> summary()




## Trying out t-test described on Page 27 of Hwang & Takane 2014 ----------

# The goal here is to examine this statement, "Similar to bootstrapping R-squared in linear regression (Ohtani 2000), we can compute the bootstrapped standard errors or confidence intervals of the difference in FIT between two models so as to examine whether there is a statistically significant difference in the FIT values. This procedure can be regarded as a nonparametric version of the paired t test for two groups of FIT values."

# It is unclear  to me whether: (1) doing a one sample t-test of the difference in two random variables against 0; versus (2) doing a paired t-test of the R^2 of two random variables  is the same or not. Here, I try to investigate this problem using the simple case of two random variables


get_one_vs_paired_comp <- function(N = 20, B = 100, mu_a = 2, mu_b = 0) {
  
  dat <- tibble::tibble(a = rnorm(N, mean = mu_a), b = rnorm(N, mean = mu_b))

  boots <- boot::boot(
    dat,
    statistic = function(dat, inds) {
      means <- colMeans(dat[inds, ])
      means['a_b'] <- means['a'] - means['b']
      return(means)
    },
    R = B
  )
  # Sanity check to make sure that there is bootstrap to bootstrap variability in statistic
  # boots$t[,1]

  # One-sample test suggested
  t_test_of_diff <- t.test(boots$t[, 3], alternative = "two.sided", mu = 0)

  # The paired t-test that was suggested
  t_test_paired <- t.test(
    boots$t[, 1],
    boots$t[, 2],
    paired = TRUE,
    alternative = "two.sided",
    mu = 0
  )

  SimDesign::nc(
    oneSampleDiff = t_test_of_diff$p.value,
    pairedDiff = t_test_paired$p.value
  )
}

notNull <- future.apply::future_replicate(
  n = 300,
  get_one_vs_paired_comp(mu_a = 2, mu_b = 0),
  simplify = TRUE,
  future.seed = TRUE
)

Null <- future.apply::future_replicate(
  300,
  get_one_vs_paired_comp(mu_a = 0, mu_b = 0),
  simplify = TRUE
)

identical(notNull['oneSampleDiff', ], notNull['pairedDiff', ])
identical(Null['oneSampleDiff', ], Null['pairedDiff', ])

# TODO: From a theoretical perspective, is a one-sample t-test of (a-b) != 0 identical to a paired t-test of a = b?

# Simulation -------------------------------------------------------------

