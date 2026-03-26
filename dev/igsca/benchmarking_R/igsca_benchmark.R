library(cSEM.DGP)
devtools::load_all()

if (FALSE) {
  Sxi1 <- matrix(c(1, 0.4, -0.3, 0.4, 1, 0.4, -0.3, 0.4, 1), 3, 3)
  w1 <- c(0.4, 0.3, 0.2)

  w1 %*% Sxi1 %*% w1

  (cmp_canonical_weights <- w1 / c(sqrt(w1 %*% Sxi1 %*% w1)))

  (single_indicator_weight <- .4 / sqrt(.4 %*% 1 %*% .4))
}

# --- Small-p DGP and Model (1 factor + 1 composite = 6 indicators) -----------

dgp_smallp <- '
  xi1 =~ 0.6*x1_1 + 0.8*x1_2 + 0.7*x1_3

  xi2 <~ 0.4*x2_1 + 0.3*x2_2 + 0.2*x2_3
  x2_1 ~~ 0.4*x2_2 + -0.3*x2_3
  x2_2 ~~ 0.4*x2_3

  xi2 ~ 0.5*xi1
'

mod_smallp <- '
  xi1 =~ x1_1 + x1_2 + x1_3
  xi2 <~ x2_1 + x2_2 + x2_3
  xi2 ~ xi1
'

# --- Big-p DGP and Model (10 factors + 10 composites = 60 indicators) --------

dgp_bigp <- '
  xi1  =~ 0.6*x1_1  + 0.8*x1_2  + 0.7*x1_3
  xi3  =~ 0.6*x3_1  + 0.8*x3_2  + 0.7*x3_3
  xi5  =~ 0.6*x5_1  + 0.8*x5_2  + 0.7*x5_3
  xi7  =~ 0.6*x7_1  + 0.8*x7_2  + 0.7*x7_3
  xi9  =~ 0.6*x9_1  + 0.8*x9_2  + 0.7*x9_3
  xi11 =~ 0.6*x11_1 + 0.8*x11_2 + 0.7*x11_3
  xi13 =~ 0.6*x13_1 + 0.8*x13_2 + 0.7*x13_3
  xi15 =~ 0.6*x15_1 + 0.8*x15_2 + 0.7*x15_3
  xi17 =~ 0.6*x17_1 + 0.8*x17_2 + 0.7*x17_3
  xi19 =~ 0.6*x19_1 + 0.8*x19_2 + 0.7*x19_3

  xi2  <~ 0.4*x2_1  + 0.3*x2_2  + 0.2*x2_3
  xi4  <~ 0.4*x4_1  + 0.3*x4_2  + 0.2*x4_3
  xi6  <~ 0.4*x6_1  + 0.3*x6_2  + 0.2*x6_3
  xi8  <~ 0.4*x8_1  + 0.3*x8_2  + 0.2*x8_3
  xi10 <~ 0.4*x10_1 + 0.3*x10_2 + 0.2*x10_3
  xi12 <~ 0.4*x12_1 + 0.3*x12_2 + 0.2*x12_3
  xi14 <~ 0.4*x14_1 + 0.3*x14_2 + 0.2*x14_3
  xi16 <~ 0.4*x16_1 + 0.3*x16_2 + 0.2*x16_3
  xi18 <~ 0.4*x18_1 + 0.3*x18_2 + 0.2*x18_3
  xi20 <~ 0.4*x20_1 + 0.3*x20_2 + 0.2*x20_3

  x2_1  ~~ 0.4*x2_2  + -0.3*x2_3
  x2_2  ~~ 0.4*x2_3
  x4_1  ~~ 0.4*x4_2  + -0.3*x4_3
  x4_2  ~~ 0.4*x4_3
  x6_1  ~~ 0.4*x6_2  + -0.3*x6_3
  x6_2  ~~ 0.4*x6_3
  x8_1  ~~ 0.4*x8_2  + -0.3*x8_3
  x8_2  ~~ 0.4*x8_3
  x10_1 ~~ 0.4*x10_2 + -0.3*x10_3
  x10_2 ~~ 0.4*x10_3
  x12_1 ~~ 0.4*x12_2 + -0.3*x12_3
  x12_2 ~~ 0.4*x12_3
  x14_1 ~~ 0.4*x14_2 + -0.3*x14_3
  x14_2 ~~ 0.4*x14_3
  x16_1 ~~ 0.4*x16_2 + -0.3*x16_3
  x16_2 ~~ 0.4*x16_3
  x18_1 ~~ 0.4*x18_2 + -0.3*x18_3
  x18_2 ~~ 0.4*x18_3
  x20_1 ~~ 0.4*x20_2 + -0.3*x20_3
  x20_2 ~~ 0.4*x20_3

  xi2  ~ 0.5*xi1
  xi4  ~ 0.5*xi3
  xi6  ~ 0.5*xi5
  xi8  ~ 0.5*xi7
  xi10 ~ 0.5*xi9
  xi12 ~ 0.5*xi11
  xi14 ~ 0.5*xi13
  xi16 ~ 0.5*xi15
  xi18 ~ 0.5*xi17
  xi20 ~ 0.5*xi19
'

mod_bigp <- '
  xi1  =~ x1_1  + x1_2  + x1_3
  xi3  =~ x3_1  + x3_2  + x3_3
  xi5  =~ x5_1  + x5_2  + x5_3
  xi7  =~ x7_1  + x7_2  + x7_3
  xi9  =~ x9_1  + x9_2  + x9_3
  xi11 =~ x11_1 + x11_2 + x11_3
  xi13 =~ x13_1 + x13_2 + x13_3
  xi15 =~ x15_1 + x15_2 + x15_3
  xi17 =~ x17_1 + x17_2 + x17_3
  xi19 =~ x19_1 + x19_2 + x19_3

  xi2  <~ x2_1  + x2_2  + x2_3
  xi4  <~ x4_1  + x4_2  + x4_3
  xi6  <~ x6_1  + x6_2  + x6_3
  xi8  <~ x8_1  + x8_2  + x8_3
  xi10 <~ x10_1 + x10_2 + x10_3
  xi12 <~ x12_1 + x12_2 + x12_3
  xi14 <~ x14_1 + x14_2 + x14_3
  xi16 <~ x16_1 + x16_2 + x16_3
  xi18 <~ x18_1 + x18_2 + x18_3
  xi20 <~ x20_1 + x20_2 + x20_3

  xi2  ~ xi1
  xi4  ~ xi3
  xi6  ~ xi5
  xi8  ~ xi7
  xi10 ~ xi9
  xi12 ~ xi11
  xi14 ~ xi13
  xi16 ~ xi15
  xi18 ~ xi17
  xi20 ~ xi19
'

# --- Benchmark Functions ------------------------------------------------------

bigN_smallp <- function(.dgp = dgp_smallp, .mod = mod_smallp) {
  set.seed(42)
  dat <- cSEM.DGP::generateData(.dgp, .empirical = TRUE, .N = 1000)
  csem(
    .data = dat,
    .model = .mod,
    .approach_weights = "GSCA",
    .disattenuate = TRUE,
    .conv_criterion = "sum_diff_absolute",
    .GSCA_modes = "CCMP",
    .tolerance = 0.0000001,
    .iter_max = 1000
  )
}

smallN_smallp <- function(.dgp = dgp_smallp, .mod = mod_smallp) {
  set.seed(42)
  dat <- cSEM.DGP::generateData(.dgp, .empirical = TRUE, .N = 50)
  csem(
    .data = dat,
    .model = .mod,
    .approach_weights = "GSCA",
    .disattenuate = TRUE,
    .conv_criterion = "sum_diff_absolute",
    .GSCA_modes = "CCMP",
    .tolerance = 0.0000001,
    .iter_max = 1000
  )
}

smallN_bigp <- function(.dgp = dgp_bigp, .mod = mod_bigp) {
  set.seed(42)
  dat <- cSEM.DGP::generateData(.dgp, .empirical = TRUE, .N = 50)
  csem(
    .data = dat,
    .model = .mod,
    .approach_weights = "GSCA",
    .disattenuate = TRUE,
    .conv_criterion = "sum_diff_absolute",
    .GSCA_modes = "CCMP",
    .tolerance = 0.0000001,
    .iter_max = 1000
  )
}

bigN_bigp <- function(.dgp = dgp_bigp, .mod = mod_bigp) {
  set.seed(42)
  dat <- cSEM.DGP::generateData(.dgp, .empirical = TRUE, .N = 1000)
  csem(
    .data = dat,
    .model = .mod,
    .approach_weights = "GSCA",
    .disattenuate = TRUE,
    .conv_criterion = "sum_diff_absolute",
    .GSCA_modes = "CCMP",
    .tolerance = 0.0000001,
    .iter_max = 1000
  )
}
