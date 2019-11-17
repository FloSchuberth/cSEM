\donttest{# NOTE: to run the example. Download and load the newst version of cSEM.DGP
# from GitHub using devtools::install_github("M-E-Rademaker/cSEM.DGP").

# Create two data generating processes (DGPs) that only differ in how the composite
# X is build. Hence, the two groups are not compositionally invariant.
dgp1 <- "
# Structural model
Y ~ 0.6*X

# Measurement model
Y =~ 1*y1
X <~ 0.4*x1 + 0.8*x2

x1 ~~ 0.3125*x2
"

dgp2 <- "
# Structural model
Y ~ 0.6*X

# Measurement model
Y =~ 1*y1
X <~ 0.8*x1 + 0.4*x2

x1 ~~ 0.3125*x2
"

g1 <- generateData(dgp1, .empirical = TRUE) # requires cSEM.DGP 
g2 <- generateData(dgp2, .empirical = TRUE) # requires cSEM.DGP

# Model is the same for both DGPs
model <- "
# Structural model
Y ~ X

# Measurement model
Y =~ y1
X <~ x1 + x2
"

# Estimate
csem_results <- csem(.data = list("group1" = g1, "group2" = g2), model)

# Test
testMICOM(csem_results, .R = 50, .alpha = c(0.01, 0.05))
}