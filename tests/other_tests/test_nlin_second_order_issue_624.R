devtools::load_all()

# Simulate data from two higher order variables X and Z
N <- 10000L
Sigma <- diag(2)
Sigma[1, 2] <- Sigma[2, 1] <- 0.3
Mu <- rep(0, 2L)

set.seed(2394827)
XI <- mvtnorm::rmvnorm(N, mean = Mu, sigma = Sigma)

disturbance <- \(eps) rnorm(N, mean = 0, sd = sqrt(eps))

# Second Order LVs
X <- XI[, 1]
Z <- XI[, 2]
Y <- 0.2 * X + 0.5 * Z + 0.5 * X * Z
Y <- Y + rnorm(N, mean = 0, sd = sqrt(1 - var(Y)))

# First Order LVs
X1 <- 0.9 * X + disturbance(sqrt(0.1))
X2 <- 0.8 * X + disturbance(sqrt(0.2))

Z1 <- 0.9 * Z + disturbance(sqrt(0.1))
Z2 <- 0.8 * Z + disturbance(sqrt(0.2))

# OV Indicators
x1 <- 0.9 * X1 + disturbance(sqrt(0.1))
x2 <- 0.8 * X1 + disturbance(sqrt(0.2))

x3 <- 0.9 * X2 + disturbance(sqrt(0.1))
x4 <- 0.8 * X2 + disturbance(sqrt(0.2))

z1 <- 0.9 * Z1 + disturbance(sqrt(0.1))
z2 <- 0.8 * Z1 + disturbance(sqrt(0.2))

z3 <- 0.9 * Z2 + disturbance(sqrt(0.1))
z4 <- 0.8 * Z2 + disturbance(sqrt(0.2))

Y1 <- 0.9 * Y + disturbance(sqrt(0.1))
Y2 <- 0.8 * Y + disturbance(sqrt(0.2))

y1 <- 0.9 * Y1 + disturbance(sqrt(0.1))
y2 <- 0.8 * Y1 + disturbance(sqrt(0.2))
y3 <- 0.9 * Y2 + disturbance(sqrt(0.1))
y4 <- 0.8 * Y2 + disturbance(sqrt(0.2))

data <- data.frame(x1, x2, x3, x4,
                   z1, z2, z3, z4,
                   y1, y2, y3, y4)

syntax <- '
  LX1 =~ x1 + x2
  LX2 =~ x3 + x4

  LZ1 =~ z1 + z2
  LZ2 =~ z3 + z4

  LY1 =~ y1 + y2
  LY2 =~ y3 + y4

  LX =~ LX1 + LX2
  LZ =~ LZ1 + LZ2
  LY =~ LY1 + LY2

  LY ~ LX + LZ + LX.LZ
'

testthat::expect_no_error({
  fit <- csem(data, syntax)
  summarize(fit)
})
