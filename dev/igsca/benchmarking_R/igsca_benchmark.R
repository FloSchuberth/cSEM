library(cSEM.DGP)
devtools::load_all()

Sxi1 <- matrix(c(1, 0.4, -0.3, 0.4, 1, 0.4, -0.3, 0.4, 1), 3, 3)
w1 <- c(0.4, 0.3, 0.2)

w1 %*% Sxi1 %*% w1

(cmp_canonical_weights <- w1 / c(sqrt(w1 %*% Sxi1 %*% w1)))

(single_indicator_weight <- .4 / sqrt(.4 %*% 1 %*% .4))

triF_triC <- 'xi1=~0.6*x11 + 0.8*x12 + 0.7*x13
               
               xi2<~0.4*x21 + 0.3*x22 + 0.2*x23
               x21~~0.4*x22 + -0.3*x23
               x22~~0.4*x23
               
               xi2 ~ 0.5*xi1'

mod <- 'xi1=~x11 + x12 + x13
        xi2<~x21 + x22 + x23
        xi2 ~ xi1'


bigN_smallp <- function(.dgp = triF_triC, .mod = mod) {
  dat <- cSEM.DGP::generateData(.dgp, .empirical = TRUE, .N = 1000)
  csem(
    .data = dat,
    .model = .mod,
    .approach_weights = 'GSCA',
    .disattenuate = TRUE,
    .conv_criterion = "sum_diff_absolute",
    .GSCA_modes = "CCMP",
    .tolerance = 0.0000001,
    .iter_max = 1000
  )
}
